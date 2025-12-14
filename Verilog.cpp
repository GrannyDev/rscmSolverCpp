//
// Created by smith on 09/12/2025.
//

#include "Verilog.h"

#include <bitset>
#include <filesystem>
#include <iostream>
#include <ranges>

VerilogGenerator::VerilogGenerator
    (
        const RSCM& solutionNode,
        const std::string& outputUri,
        const std::vector<Layer>& layers,
        const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap,
        const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap,
        const unsigned int inputBW,
        const std::vector<int>& targets,
        const std::vector<std::pair<int, std::vector<DAG>>>& scmDesigns,
        const bool overwrite
    )
    : solutionNode_(solutionNode), layers_(layers), idxToVarMap_(idxToVarMap), varToIdxMap_(varToIdxMap), inputBW_(inputBW), targets_(targets), scmDesigns_(scmDesigns), varDefsTostring_(VarDefsToString())
{
    // check if output already exists
    if (std::filesystem::exists(outputUri) && !overwrite)
    {
        throw std::runtime_error("Output file already exists.");
    }
    // overwrite
    outputFile_ = std::ofstream(outputUri, std::ofstream::trunc);
    muxCounter_ = 0;

    // Generate modules and testbench
    HandleModules();
    HandleTopModule();
    // PrintTestbench();

    // print the bitwidth of each parameters of each adder
    unsigned int adderIdx = 0;
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            std::cout << "Adder " << adderIdx << " parameters bitwidths: " << std::endl;
            for (size_t p = 0; p < adder.variables.size(); p++)
            {
                const auto varType = idxToVarMap_.at(p);
                const auto bw = solutionNode_.variableBitWidths[adderIdx * adder.variables.size() + p];
                std::cout << "\t" << varDefsTostring_(varType) << ": " << bw << " bits" << std::endl;
            }
            ++adderIdx;
        }
    }

    // print the min shit savings per parameter
    adderIdx = 0;
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            std::cout << "Adder " << adderIdx << " bit shift savings: " << std::endl;
            for (size_t p = 0; p < adder.variables.size(); p++)
            {
                const auto varType = idxToVarMap_.at(p);
                const auto bw = solutionNode_.minShiftSavings[adderIdx * adder.variables.size() + p];
                std::cout << "\t" << varDefsTostring_(varType) << ": " << bw << " bits" << std::endl;
            }
            ++adderIdx;
        }
    }
}

void VerilogGenerator::HandleModules()
{
    unsigned int nbAdder = 0;
    for (unsigned int layerIdx = 0; layerIdx < layers_.size(); layerIdx++) {
        for (unsigned int adderIdx = 0; adderIdx < layers_[layerIdx].adders.size(); adderIdx++) {
            std::stringstream adderStream = HandleAdder(layerIdx, adderIdx, nbAdder);
            outputFile_ << adderStream.str();
            nbAdder++;
        }
    }
}

void VerilogGenerator::HandleTopModule()
{
    std::stringstream topStream;
    std::vector<std::pair<unsigned int, unsigned int>> allMuxInputs; // bitPos, bitwidth
    
    // Collect all MUX inputs by processing parameters
    unsigned int bitPos = 0;
    unsigned int nbAdder = 0;
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            unsigned int paramIdx = 0;
            for (const auto& param : adder.variables) {
                VariableDefs paramType = idxToVarMap_.at(paramIdx);
                
                // Count mux inputs for this parameter
                unsigned int nbMuxInputs = 0;
                for (const int v : param.possibleValuesFusion) {
                    if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
                        ++nbMuxInputs;
                    }
                }
                
                if (nbMuxInputs > 1) {
                    unsigned int muxBW = ComputeBitWidth(static_cast<int>(nbMuxInputs - 1));
                    allMuxInputs.emplace_back(bitPos, muxBW);
                }
                
                bitPos += param.possibleValuesFusion.size();
                paramIdx++;
            }
            nbAdder++;
        }
    }
    
    // Generate RSCM module header
    topStream << "module RSCM (\n";
    topStream << "\tinput" << PrintWireBitWidth(inputBW_) << "X,\n";
    
    // Add all MUX inputs
    for (const auto& [muxBitPos, muxBW] : allMuxInputs) {
        topStream << "\tinput" << PrintWireBitWidth(muxBW) << "MUX_INPUT_" << muxBitPos << ",\n";
    }
    
    // Output is the last adder's output
    unsigned int finalOutputBW = GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, nbAdder - 1);
    topStream << "\toutput" << PrintWireBitWidth(finalOutputBW) << "OUTPUT\n";
    topStream << ");\n\n";
    
    // Generate routing and instantiation for each adder
    bitPos = 0;
    nbAdder = 0;
    for (unsigned int layerIdx = 0; layerIdx < layers_.size(); layerIdx++) {
        for (unsigned int adderIdx = 0; adderIdx < layers_[layerIdx].adders.size(); adderIdx++) {
            topStream << HandleAdderRouting(layerIdx, adderIdx, nbAdder, bitPos).str();
            nbAdder++;
        }
    }
    
    // Final output assignment
    topStream << "\t// Final Output\n";
    topStream << "\tassign OUTPUT = ADDER_" << (nbAdder - 1) << "_OUTPUT;\n\n";
    topStream << "endmodule\n";
    
    outputFile_ << topStream.str();
}

std::stringstream VerilogGenerator::HandleAdder(const unsigned int layerIndex, const unsigned int adderIdx, const unsigned int nbAdder)
{
    std::stringstream adderStream;
    std::vector<std::pair<std::string, unsigned int>> allMuxRanges;
    
    // Calculate the global parameter index offset for this adder
    // Sum up all parameters from all previous adders (adders 0 to nbAdder-1)
    unsigned int paramGlobalIdxOffset = 0;
    unsigned int addersProcessed = 0;
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            if (addersProcessed >= nbAdder) {
                break;
            }
            paramGlobalIdxOffset += adder.variables.size();
            addersProcessed++;
        }
        if (addersProcessed >= nbAdder) {
            break;
        }
    }

    // Generate module header
    unsigned int leftInputBW = GetParamBitWidth(VariableDefs::LEFT_INPUTS, nbAdder);
    unsigned int rightInputBW = GetParamBitWidth(VariableDefs::RIGHT_INPUTS, nbAdder);
    if (layerIndex == 0)
    {
        leftInputBW = inputBW_;
        rightInputBW = inputBW_;
    }
    adderStream << "module ADDER_" << nbAdder << " (\n";
    adderStream << "\tinput" << PrintWireBitWidth(leftInputBW) << "LEFT_INPUTS_ADDER_" << nbAdder << ",\n";
    adderStream << "\tinput" << PrintWireBitWidth(rightInputBW) << "RIGHT_INPUTS_ADDER_" << nbAdder << ",\n";
    
    // Add mux input ports (we'll determine these as we process parameters)
    std::stringstream muxInputPorts;
    std::stringstream insideLogic;
    
    // Process each parameter
    // Calculate the bit position offset (sum of possibleValuesFusion.size() from all previous adders)
    unsigned int bitPosOffset = 0;
    unsigned int addersProc = 0;
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            if (addersProc >= nbAdder) {
                break;
            }
            for (const auto& var : adder.variables) {
                bitPosOffset += var.possibleValuesFusion.size();
            }
            addersProc++;
        }
        if (addersProc >= nbAdder) {
            break;
        }
    }
    
    unsigned int paramIdx = 0;
    unsigned int bitPos = bitPosOffset;
    VariableDefs previousParamType; // Default starting type
    unsigned int previousParamBitWidth = 0;

    for (const auto& param : layers_[layerIndex].adders[adderIdx].variables)
    {
        unsigned int paramGlobalIdx = paramGlobalIdxOffset + paramIdx;
        VariableDefs currentParamType = idxToVarMap_.at(paramIdx);

        // Determine the correct previous parameter type based on current parameter
        if (currentParamType == VariableDefs::LEFT_SHIFTS) {
            previousParamType = VariableDefs::LEFT_INPUTS;
            previousParamBitWidth = GetParamBitWidth(VariableDefs::LEFT_INPUTS, nbAdder);
        } else if (currentParamType == VariableDefs::RIGHT_SHIFTS) {
            previousParamType = VariableDefs::RIGHT_INPUTS;
            previousParamBitWidth = GetParamBitWidth(VariableDefs::RIGHT_INPUTS, nbAdder);
        } else if (currentParamType == VariableDefs::RIGHT_MULTIPLIER) {
            previousParamType = VariableDefs::RIGHT_SHIFTS;
            previousParamBitWidth = GetParamBitWidth(VariableDefs::RIGHT_MULTIPLIER, nbAdder);
        } else if (currentParamType == VariableDefs::OUTPUTS_SHIFTS) {
            previousParamType = VariableDefs::RIGHT_MULTIPLIER;
            previousParamBitWidth = GetParamBitWidth(VariableDefs::RIGHT_MULTIPLIER, nbAdder);
        } else
        {
            paramIdx++;
            bitPos += param.possibleValuesFusion.size();
            continue;
        }

        // Handle the parameter (pass bitPos for bitset indexing, paramGlobalIdx for variableBitWidths)
        auto muxRanges = HandleParam(param, paramIdx, paramGlobalIdx, bitPos, previousParamBitWidth, nbAdder, previousParamType, insideLogic);

        // Collect mux ranges
        allMuxRanges.insert(allMuxRanges.end(), muxRanges.begin(), muxRanges.end());

        // Add mux input ports for this parameter if needed
        // BUT skip LEFT_INPUTS and RIGHT_INPUTS - those muxes are only in the top module for routing
        if (!muxRanges.empty() && currentParamType != VariableDefs::LEFT_INPUTS && currentParamType != VariableDefs::RIGHT_INPUTS) {
            for (const auto& val : muxRanges | std::views::values) {
                muxInputPorts << "\tinput" << PrintWireBitWidth(val) << "MUX_INPUT_" << bitPos << ",\n";
            }
        }
        paramIdx++;
        bitPos += param.possibleValuesFusion.size();
    }

    // Add output port
    const unsigned int outputBW = GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, nbAdder);
    adderStream << muxInputPorts.str();
    adderStream << "\toutput wire" << PrintWireBitWidth(outputBW) << "OUTPUT\n";
    adderStream << ");\n\n";

    adderStream << insideLogic.str();
    
    // Connect OUTPUT to OUTPUTS_SHIFTS
    adderStream << "\tassign OUTPUT = OUTPUTS_SHIFTS_ADDER_" << nbAdder << ";\n";
    
    adderStream << "endmodule\n\n";

    return adderStream;
}

std::vector<std::pair<std::string, unsigned int>> VerilogGenerator::HandleParam(
    const Variables& param,
    const unsigned int paramIdx,
    const unsigned int paramGlobalIdx,
    const unsigned int bitPos,
    const unsigned int inputBw,
    const unsigned int adderIdx,
    const VariableDefs& previousParamType,
    std::stringstream& inStream
    )
{
    // get param type and bitwidth
    const VariableDefs paramType = idxToVarMap_.at(paramIdx);
    const unsigned int bw = solutionNode_.variableBitWidths[paramGlobalIdx];
    std::vector<std::pair<std::string, unsigned int>> muxRanges;

    bool needsSlicing = false;
    bool needsShift = false;

    // get if is mux or not and maxShift
    unsigned int nbMuxInputs = 0;
    unsigned int maxShift = 0;
    for (const int v : param.possibleValuesFusion) {
        if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
            ++nbMuxInputs;
            if (paramType == VariableDefs::OUTPUTS_SHIFTS) {
                // Output shift uses power of 2, but only count if v > 0
                if (v > 0) {
                    maxShift = std::max(maxShift, static_cast<unsigned int>(1 << v));
                }
            } else {
                // Input shift uses direct value
                maxShift = std::max(maxShift, static_cast<unsigned int>(v));
            }
        }
    }

    // update needsSlicing and needsShift
    if (maxShift > 0)
    {
        needsShift = true;
        if (maxShift + inputBw > bw) {
            needsSlicing = true;
        }
    }

    // now we can actually print the verilog
    // Handle ADD_SUB (RIGHT_MULTIPLIER) specially
    if (paramType == VariableDefs::RIGHT_MULTIPLIER)
    {
        // Get the bitwidths of LEFT_SHIFTS and RIGHT_SHIFTS for this adder
        unsigned int leftShiftBW = GetParamBitWidth(VariableDefs::LEFT_SHIFTS, adderIdx);
        unsigned int rightShiftBW = GetParamBitWidth(VariableDefs::RIGHT_SHIFTS, adderIdx);
        unsigned int addSubBW = std::max(leftShiftBW, rightShiftBW);
        unsigned int concatBW = addSubBW + 1; // Concatenation of {cout, ADD_SUB}
        
        // For adder, we need to add or subtract LEFT_SHIFTS and RIGHT_SHIFTS
        if (nbMuxInputs < 2)
        {
            // Single operation (either + or -)
            // Determine which operation by checking the bitset
            bool isSubtract = false;
            for (const int v : param.possibleValuesFusion) {
                if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
                    isSubtract = (v == -1);
                    break;
                }
            }
            const char op = isSubtract ? '-' : '+';
            
            // Check if we need slicing
            if (concatBW != bw)
            {
                // Need to slice the concatenation to match solver's expected bitwidth
                // Declare cout and ADD_SUB_INTERNAL wires
                inStream << "\twire cout_ADDER_" << adderIdx << ";\n";
                inStream << "\twire" << PrintWireBitWidth(addSubBW) << "ADD_SUB_ADDER_" << adderIdx << "_INTERNAL;\n";
                // Cast operands to max bitwidth to avoid width expansion warnings
                inStream << "\tassign {cout_ADDER_" << adderIdx << ", ADD_SUB_ADDER_" << adderIdx << "_INTERNAL} = ";
                if (leftShiftBW < addSubBW) {
                    inStream << addSubBW << "'(LEFT_SHIFTS_ADDER_" << adderIdx << ")";
                } else {
                    inStream << "LEFT_SHIFTS_ADDER_" << adderIdx;
                }
                inStream << " " << op << " ";
                if (rightShiftBW < addSubBW) {
                    inStream << addSubBW << "'(RIGHT_SHIFTS_ADDER_" << adderIdx << ")";
                } else {
                    inStream << "RIGHT_SHIFTS_ADDER_" << adderIdx;
                }
                inStream << ";\n";
                inStream << "\twire" << PrintWireBitWidth(concatBW) << "ADD_SUB_CONCAT_ADDER_" << adderIdx << ";\n";
                inStream << "\tassign ADD_SUB_CONCAT_ADDER_" << adderIdx << " = {cout_ADDER_" << adderIdx << ", ADD_SUB_ADDER_" << adderIdx << "_INTERNAL};\n";
                inStream << "\twire" << PrintWireBitWidth(bw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << ";\n";
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = ADD_SUB_CONCAT_ADDER_" << adderIdx << "[" << (bw - 1) << ":0];\n";
            } else
            {
                // No slicing needed, concatenation matches expected bitwidth
                inStream << "\twire cout_ADDER_" << adderIdx << ";\n";
                inStream << "\twire" << PrintWireBitWidth(addSubBW) << "ADD_SUB_ADDER_" << adderIdx << ";\n";
                // Cast operands to max bitwidth to avoid width expansion warnings
                inStream << "\tassign {cout_ADDER_" << adderIdx << ", ADD_SUB_ADDER_" << adderIdx << "} = ";
                if (leftShiftBW < addSubBW) {
                    inStream << addSubBW << "'(LEFT_SHIFTS_ADDER_" << adderIdx << ")";
                } else {
                    inStream << "LEFT_SHIFTS_ADDER_" << adderIdx;
                }
                inStream << " " << op << " ";
                if (rightShiftBW < addSubBW) {
                    inStream << addSubBW << "'(RIGHT_SHIFTS_ADDER_" << adderIdx << ")";
                } else {
                    inStream << "RIGHT_SHIFTS_ADDER_" << adderIdx;
                }
                inStream << ";\n";
                inStream << "\twire" << PrintWireBitWidth(bw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << ";\n";
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = {cout_ADDER_" << adderIdx << ", ADD_SUB_ADDER_" << adderIdx << "};\n";
            }
        } else
        {
            // Multiplexer between add and subtract
            std::string muxName = "MUX_" + std::to_string(muxCounter_);
            muxRanges.emplace_back(muxName, ComputeBitWidth(static_cast<int>(nbMuxInputs - 1)));
            
            // Declare cout and ADD_SUB wires as regs (will be assigned by mux)
            inStream << "\treg cout_ADDER_" << adderIdx << ";\n";
            inStream << "\treg" << PrintWireBitWidth(addSubBW) << "ADD_SUB_ADDER_" << adderIdx << "_INTERNAL;\n";
            
            // Check if we need slicing
            if (concatBW != bw)
            {
                // Need to slice the concatenation to match solver's expected bitwidth
                inStream << "\twire" << PrintWireBitWidth(concatBW) << "ADD_SUB_CONCAT_ADDER_" << adderIdx << ";\n";
                inStream << "\tassign ADD_SUB_CONCAT_ADDER_" << adderIdx << " = {cout_ADDER_" << adderIdx << ", ADD_SUB_ADDER_" << adderIdx << "_INTERNAL};\n";
                inStream << "\twire" << PrintWireBitWidth(bw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << ";\n";
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = ADD_SUB_CONCAT_ADDER_" << adderIdx << "[" << (bw - 1) << ":0];\n";
                HandleMultiplexer(inStream, paramType, param, "LEFT_SHIFTS_ADDER_" + std::to_string(adderIdx), "ADD_SUB_ADDER_" + std::to_string(adderIdx) + "_INTERNAL", bitPos, nbMuxInputs, adderIdx, addSubBW);
            } else
            {
                // No slicing needed - concatenation matches expected bitwidth
                // The mux will assign to the internal wire, then we concatenate to get the final result
                inStream << "\twire" << PrintWireBitWidth(bw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << ";\n";
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = {cout_ADDER_" << adderIdx << ", ADD_SUB_ADDER_" << adderIdx << "_INTERNAL};\n";
                HandleMultiplexer(inStream, paramType, param, "LEFT_SHIFTS_ADDER_" + std::to_string(adderIdx), "ADD_SUB_ADDER_" + std::to_string(adderIdx) + "_INTERNAL", bitPos, nbMuxInputs, adderIdx, addSubBW);
            }
        }
    }
    // 8 cases needsShift * needsSlicing * nbMuxInputs > 1 ...
    else if (!needsShift)
    {
        inStream << "\twire" << PrintWireBitWidth(bw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << ";\n";
        inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = " << varDefsTostring_(previousParamType) << "_ADDER_" << adderIdx << ";\n";
    } else
    {
        if (nbMuxInputs < 2)
        {
            inStream << "\twire" << PrintWireBitWidth(bw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << ";\n";
            if (!needsSlicing)
            {
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = " << varDefsTostring_(previousParamType) << "_ADDER_" << adderIdx << " << " << maxShift << ";\n";
            } else
            {
                inStream << "\twire" << PrintWireBitWidth(maxShift + inputBw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << "_UNSLICED;\n";
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << "_UNSLICED = " << varDefsTostring_(previousParamType) << "_ADDER_" << adderIdx << " << " << maxShift << ";\n";
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << "_UNSLICED[" << (bw - 1) << ":0];\n";
            }
        } else
        {
            // this is an actual mux
            std::string muxName = "MUX_" + std::to_string(muxCounter_);
            muxRanges.emplace_back(muxName, ComputeBitWidth(static_cast<int>(nbMuxInputs - 1)));
            inStream << "\treg" << PrintWireBitWidth(bw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << ";\n";
            if (needsSlicing)
            {
                inStream << "\treg" << PrintWireBitWidth(maxShift + inputBw) << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << "_UNSLICED;\n";
                HandleMultiplexer(inStream, paramType, param, varDefsTostring_(previousParamType) + "_ADDER_" + std::to_string(adderIdx), varDefsTostring_(paramType) + "_ADDER_" + std::to_string(adderIdx) + "_UNSLICED", bitPos, nbMuxInputs, adderIdx, maxShift + inputBw);
                inStream << "\tassign " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << " = " << varDefsTostring_(paramType) << "_ADDER_" << adderIdx << "_UNSLICED[" << (bw - 1) << ":0];\n";
            } else
            {
                HandleMultiplexer(inStream, paramType, param, varDefsTostring_(previousParamType) + "_ADDER_" + std::to_string(adderIdx), varDefsTostring_(paramType) + "_ADDER_" + std::to_string(adderIdx), bitPos, nbMuxInputs, adderIdx, bw);
            }
        }
    }
    return muxRanges;
}

std::stringstream& VerilogGenerator::HandleMultiplexer(
    std::stringstream& muxStream,
    const VariableDefs paramName,
    const Variables& param,
    const std::string& inputWire,
    const std::string& outputWire,
    const unsigned int bitPos,
    const unsigned int nbMuxInputs,
    const unsigned int adderIdx,
    const unsigned int outputBW
    )
{
    const bool isAdder = paramName == VariableDefs::RIGHT_MULTIPLIER;
    const bool isInput = paramName == VariableDefs::LEFT_INPUTS || paramName == VariableDefs::RIGHT_INPUTS;
    const bool isOutput = paramName == VariableDefs::OUTPUTS_SHIFTS;

    // Lambda to generate assignment statement
    auto generateAssignment = [&](const int v) {
        if (isAdder) {
            const char op = v == -1 ? '-' : '+';
            const std::string leftWire = "LEFT_SHIFTS_ADDER_" + std::to_string(adderIdx);
            const std::string rightWire = "RIGHT_SHIFTS_ADDER_" + std::to_string(adderIdx);
            const std::string coutWire = "cout_ADDER_" + std::to_string(adderIdx);
            // Get bitwidths of left and right for proper extension
            unsigned int leftBW = GetParamBitWidth(VariableDefs::LEFT_SHIFTS, adderIdx);
            unsigned int rightBW = GetParamBitWidth(VariableDefs::RIGHT_SHIFTS, adderIdx);
            unsigned int maxBW = std::max(leftBW, rightBW);
            // Cast operands to match max bitwidth
            muxStream << "{" << coutWire << ", " << outputWire << "} = ";
            if (leftBW < maxBW) {
                muxStream << std::to_string(maxBW) << "'(" << leftWire << ")";
            } else {
                muxStream << leftWire;
            }
            muxStream << " " << op << " ";
            if (rightBW < maxBW) {
                muxStream << std::to_string(maxBW) << "'(" << rightWire << ")";
            } else {
                muxStream << rightWire;
            }
            muxStream << ";\n";
        } else if (isOutput) {
            // Output shift multiplexer - shift and cast to output bitwidth
            if (v == 0) {
                muxStream << outputWire << " = " << std::to_string(outputBW) << "'(" << inputWire << ");\n";
            } else {
                muxStream << outputWire << " = " << std::to_string(outputBW) << "'(" << inputWire << " << " << (1 << v) << ");\n";
            }
        } else if (!isInput) {
            // Shift multiplexers (LEFT/RIGHT_INPUT shifts) - cast to output bitwidth
            if (v == 0) {
                muxStream << outputWire << " = " << std::to_string(outputBW) << "'(" << inputWire << ");\n";
            } else {
                muxStream << outputWire << " = " << std::to_string(outputBW) << "'(" << inputWire << " << " << v << ");\n";
            }
        }
    };

    unsigned int currentSelBit = 0;
    unsigned int firstV = 0;
    muxStream << "\talways @ * begin\n\t\tcase (MUX_INPUT_" << bitPos << ")\n";
    for (const int v : param.possibleValuesFusion) {
        if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
            if (currentSelBit == 0) firstV = v;
            muxStream << "\t\t\t" << PrintBinaryMuxCode(nbMuxInputs, currentSelBit) << ": ";
            generateAssignment(v);
            ++currentSelBit;
        }
    }

    // if nbMuxInputs not a power of two, need to add a default case
    if (nbMuxInputs == 0 || (nbMuxInputs & nbMuxInputs - 1) != 0)
    {
        // make the zero case the default case
        muxStream << "\t\t\tdefault: ";
        generateAssignment(static_cast<int>(firstV));
    }

    muxStream << "\t\tendcase\n\tend\n";
    muxCounter_++;
    return muxStream;
}

std::stringstream VerilogGenerator::HandleAdderRouting(
    const unsigned int layerIndex,
    const unsigned int adderIdx,
    const unsigned int nbAdder,
    unsigned int& bitPos
) const
{
    std::stringstream routingStream;

    routingStream << "\t// Adder " << nbAdder << " Wires and Routing\n";

    // Declare output wire for this adder
    const unsigned int outputBW = GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, nbAdder);
    routingStream << "\twire" << PrintWireBitWidth(outputBW) << "ADDER_" << nbAdder << "_OUTPUT;\n";

    // Process LEFT_INPUTS and RIGHT_INPUTS parameters to generate routing
    unsigned int paramIdx = 0;
    for (const auto& param : layers_[layerIndex].adders[adderIdx].variables) {
        VariableDefs paramType = idxToVarMap_.at(paramIdx);

        if (paramType == VariableDefs::LEFT_INPUTS || paramType == VariableDefs::RIGHT_INPUTS) {
            std::string inputName = paramType == VariableDefs::LEFT_INPUTS ? "LEFT_INPUTS" : "RIGHT_INPUTS";
            const unsigned int inputBW = GetParamBitWidth(paramType, nbAdder);

            // Count options and collect values
            unsigned int nbOptions = 0;
            std::vector<int> selectedValues;
            for (const int v : param.possibleValuesFusion) {
                if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
                    selectedValues.push_back(v);
                    ++nbOptions;
                }
            }

            if (nbOptions == 1) {
                // Single input - direct assignment
                const int value = selectedValues[0];
                routingStream << "\twire" << PrintWireBitWidth(inputBW) << inputName << "_ADDER_" << nbAdder << ";\n";

                if (value < 0) {
                    // Input from X
                    if (int shiftAmount = std::abs(value + 1); shiftAmount == 0) {
                        routingStream << "\tassign " << inputName << "_ADDER_" << nbAdder << " = X;\n";
                    } else {
                        routingStream << "\tassign " << inputName << "_ADDER_" << nbAdder << " = X << " << shiftAmount << ";\n";
                    }
                } else {
                    // Input from previous adder
                    routingStream << "\tassign " << inputName << "_ADDER_" << nbAdder << " = ADDER_" << value << "_OUTPUT;\n";
                }
            } else if (nbOptions > 1) {
                // Multiple inputs - need multiplexer
                routingStream << "\treg" << PrintWireBitWidth(inputBW) << inputName << "_ADDER_" << nbAdder << ";\n";

                routingStream << "\talways @ * begin\n";
                routingStream << "\t\tcase (MUX_INPUT_" << bitPos << ")\n";

                unsigned int caseIdx = 0;
                const int firstValue = selectedValues[0];
                for (const int v : selectedValues) {
                    routingStream << "\t\t\t" << PrintBinaryMuxCode(nbOptions, caseIdx) << ": ";

                    if (v < 0) {
                        // Input from X
                        int shiftAmount = std::abs(v + 1);
                        if (shiftAmount == 0) {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = " << inputBW << "'(X);\n";
                        } else {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = " << inputBW << "'(X << " << shiftAmount << ");\n";
                        }
                    } else {
                        // Input from previous adder
                        unsigned int prevOutputBW = GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, v);
                        if (prevOutputBW != inputBW) {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = " << inputBW << "'(ADDER_" << v << "_OUTPUT);\n";
                        } else {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = ADDER_" << v << "_OUTPUT;\n";
                        }
                    }
                    caseIdx++;
                }

                // Add default case if not power of 2
                if ((nbOptions & nbOptions - 1) != 0) {
                    routingStream << "\t\t\tdefault: ";
                    if (firstValue < 0) {
                        int shiftAmount = std::abs(firstValue + 1);
                        if (shiftAmount == 0) {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = " << inputBW << "'(X);\n";
                        } else {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = " << inputBW << "'(X << " << shiftAmount << ");\n";
                        }
                    } else {
                        unsigned int prevOutputBW = GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, firstValue);
                        if (prevOutputBW != inputBW) {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = " << inputBW << "'(ADDER_" << firstValue << "_OUTPUT);\n";
                        } else {
                            routingStream << inputName << "_ADDER_" << nbAdder << " = ADDER_" << firstValue << "_OUTPUT;\n";
                        }
                    }
                }

                routingStream << "\t\tendcase\n";
                routingStream << "\tend\n";
            }
        }

        bitPos += param.possibleValuesFusion.size();
        paramIdx++;
    }

    // Instantiate the adder module
    routingStream << "\t// Adder " << nbAdder << " Instantiation\n";
    routingStream << "\tADDER_" << nbAdder << " adder_" << nbAdder << "_inst (\n";

    unsigned int leftInputBW = GetParamBitWidth(VariableDefs::LEFT_INPUTS, nbAdder);
    unsigned int rightInputBW = GetParamBitWidth(VariableDefs::RIGHT_INPUTS, nbAdder);

    routingStream << "\t\t.LEFT_INPUTS_ADDER_" << nbAdder << "(LEFT_INPUTS_ADDER_" << nbAdder << "),\n";
    routingStream << "\t\t.RIGHT_INPUTS_ADDER_" << nbAdder << "(RIGHT_INPUTS_ADDER_" << nbAdder << "),\n";

    // Add MUX inputs for this adder - need to recalculate bitPos for parameters after LEFT/RIGHT_INPUTS
    // Calculate bit position at the start of this adder
    unsigned int adderBitPosStart = 0;
    unsigned int processedAdders = 0;
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            if (processedAdders >= nbAdder) {
                break;
            }
            for (const auto& var : adder.variables) {
                adderBitPosStart += var.possibleValuesFusion.size();
            }
            processedAdders++;
        }
        if (processedAdders >= nbAdder) {
            break;
        }
    }

    paramIdx = 0;
    unsigned int currentBitPos = adderBitPosStart;
    for (const auto& param : layers_[layerIndex].adders[adderIdx].variables) {
        VariableDefs paramType = idxToVarMap_.at(paramIdx);
        
        unsigned int nbMuxInputs = 0;
        for (const int v : param.possibleValuesFusion) {
            if (solutionNode_.rscm.set.test(currentBitPos + v + param.zeroPoint)) {
                ++nbMuxInputs;
            }
        }

        // Only add MUX inputs for non-routing parameters (skip LEFT_INPUTS and RIGHT_INPUTS)
        if (nbMuxInputs > 1 && paramType != VariableDefs::LEFT_INPUTS && paramType != VariableDefs::RIGHT_INPUTS) {
            routingStream << "\t\t.MUX_INPUT_" << currentBitPos << "(MUX_INPUT_" << currentBitPos << "),\n";
        }

        currentBitPos += param.possibleValuesFusion.size();
        paramIdx++;
    }
    
    // Update bitPos to the end of this adder's parameters
    bitPos = currentBitPos;

    routingStream << "\t\t.OUTPUT(ADDER_" << nbAdder << "_OUTPUT)\n";
    routingStream << "\t);\n\n";

    return routingStream;
}

// helper methods
std::string VerilogGenerator::PrintWireBitWidth(const unsigned int bw)
{
    if (bw == 1) return " ";
    return " [" + std::to_string(bw - 1) + ":0] ";
}

std::string VerilogGenerator::PrintBinaryMuxCode(const unsigned int nbSelectBits, const unsigned int selecValue)
{
    const unsigned int bw = ComputeBitWidth(static_cast<int>(nbSelectBits - 1));
    return std::to_string(bw) + "'b" + std::bitset<32>(selecValue).to_string().substr(32 - bw);
}

unsigned int VerilogGenerator::ComputeBitWidth(const int value)
{
    return static_cast<unsigned int>(std::ceil(std::log2(value + 1)));
}

unsigned int VerilogGenerator::GetParamBitWidth(
    const VariableDefs varType,
    const unsigned int adderIdx) const
{
    const unsigned int nbVarsPerAdder = layers_.back().adders.back().variables.size();
    const auto baseIdx = varToIdxMap_.at(varType) + adderIdx * nbVarsPerAdder;
    return solutionNode_.variableBitWidths[baseIdx];
}