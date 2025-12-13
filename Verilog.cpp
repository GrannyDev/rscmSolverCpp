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
    : solutionNode_(solutionNode), layers_(layers), idxToVarMap_(idxToVarMap), varToIdxMap_(varToIdxMap), inputBW_(inputBW), targets_(targets), scmDesigns_(scmDesigns)
{
    // check if output already exists
    if (std::filesystem::exists(outputUri) && !overwrite)
    {
        throw std::runtime_error("Output file already exists.");
    }
    // overwrite
    outputFile_ = std::ofstream(outputUri, std::ofstream::trunc);

    // Generate modules and testbench
    PrintModules();
    // PrintTestbench();

    // print the bitwidth of each parameters of each adder
    unsigned int adderIdx = 0;
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            std::cout << "Adder " << adderIdx << " parameters bitwidths: " << std::endl;
            for (size_t p = 0; p < adder.variables.size(); p++)
            {
                const auto varType = idxToVarMap_.at(p);
                const auto bw = GetParamBitWidth(varType, static_cast<unsigned int>(adderIdx), layer.adders[0].variables.size());
                std::cout << "\t" << VarDefsToString()(varType) << ": " << bw << " bits" << std::endl;
            }
            ++adderIdx;
        }
    }

}

unsigned int VerilogGenerator::ComputeBitWidth(const int value)
{
    return static_cast<unsigned int>(std::ceil(std::log2(value + 1)));
}

unsigned int VerilogGenerator::GetParamBitWidth(
    const VariableDefs varType,
    const unsigned int adderIdx,
    const size_t nbVarsPerAdder) const
{
    const auto baseIdx = varToIdxMap_.at(varType) + adderIdx * nbVarsPerAdder;
    return solutionNode_.variableBitWidths[baseIdx];
}

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

std::string VerilogGenerator::GenerateMuxCode(
    const Variables& param,
    const size_t paramBitPos,
    const unsigned int paramGlobalIdx,
    const unsigned int nbPossibleValues,
    const std::string& muxName,
    const unsigned int bitwidth,
    const unsigned int adderIdx,
    const unsigned int inputBitWidth,
    const size_t nbVarsPerAdder
    ) const
{
    std::stringstream muxCode;
    const bool isAdder = muxName == "RIGHT_MULTIPLIER";
    const bool isInput = muxName == "LEFT_INPUT" || muxName == "RIGHT_INPUT";
    const bool isShifter = !isAdder && !isInput;
    
    // For shifters, find max shift amount
    unsigned int maxShift = 0;
    if (isShifter && inputBitWidth > 0) {
        for (const int v : param.possibleValuesFusion) {
            if (solutionNode_.rscm.set.test(paramBitPos + v + param.zeroPoint)) {
                if (muxName == "OUTPUT") {
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
    }
    
    // Only use unsliced/sliced pattern if maxShift > 0
    const bool needsSlicing = isShifter && inputBitWidth > 0 && maxShift > 0;
    
    // Determine wire name based on context
    std::string interWireName;
    std::string slicedWireName;
    if (isAdder) {
        interWireName = "OUTPUT";
    } else if (isInput) {
        interWireName = muxName + "_ADDER_" + std::to_string(adderIdx);
    } else if (muxName == "OUTPUT") {
        // Output shift multiplexer
        if (needsSlicing) {
            interWireName = "ADDER_" + std::to_string(adderIdx) + "_OUTPUT_UNSLICED";
            slicedWireName = "ADDER_" + std::to_string(adderIdx) + "_SHIFTED_OUTPUT";
        } else {
            interWireName = "ADDER_" + std::to_string(adderIdx) + "_SHIFTED_OUTPUT";
        }
    } else {
        // Shift multiplexers - muxName already contains full name like LEFT_INPUT_ADDER_X
        if (needsSlicing) {
            interWireName = muxName + "_SHIFTED_UNSLICED";
            slicedWireName = muxName + "_SHIFTED";
        } else {
            interWireName = muxName + "_SHIFTED";
        }
    }

    unsigned int nbMuxInputs = 0;
    for (const int v : param.possibleValuesFusion) {
        if (solutionNode_.rscm.set.test(paramBitPos + v + param.zeroPoint)) {
            nbMuxInputs++; // having the number of inputs before printing the case is needed for the binary conversion...
        }
    }

    // Declare wire (except for inputs which are ports, and output port)
    // For shifters, use the precomputed bitwidth
    // For adders, use max operand size + 1 for carry-out
    unsigned int declaredBitWidth = bitwidth;
    bool needsAdderSlicing = false;
    if (needsSlicing) {
        // For shifters with slicing, use precomputed bitwidth as the UNSLICED size
        declaredBitWidth = bitwidth;
    } else if (isAdder && inputBitWidth > 0) {
        // For adders, get the max of left and right shifted input sizes, then add 1 for carry
        const auto leftShiftedBW = GetParamBitWidth(VariableDefs::LEFT_SHIFTS, adderIdx, nbVarsPerAdder);
        const auto rightShiftedBW = GetParamBitWidth(VariableDefs::RIGHT_SHIFTS, adderIdx, nbVarsPerAdder);
        declaredBitWidth = std::max(leftShiftedBW, rightShiftedBW) + 1;
        // Always need slicing for adders to separate carry
        needsAdderSlicing = true;
    }
    
    if (!isInput) {
        // For adders with slicing and multiple values, declare OUTPUT_CARRY, OUTPUT, OUTPUT_UNSLICED, and OUTPUT_SLICED
        // Otherwise, normal declaration logic
        if (nbMuxInputs > 1) {
            if (needsAdderSlicing && isAdder) {
                const auto maxOperandBW = inputBitWidth;  // max of left and right operands
                muxCode << "\treg OUTPUT_CARRY;\n";
                muxCode << "\treg" << PrintWireBitWidth(maxOperandBW) << "OUTPUT;\n";
                muxCode << "\treg" << PrintWireBitWidth(declaredBitWidth) << "OUTPUT_UNSLICED;\n";
                muxCode << "\treg" << PrintWireBitWidth(bitwidth) << "OUTPUT_SLICED;\n";
            } else {
                muxCode << "\treg" << PrintWireBitWidth(declaredBitWidth) << interWireName << ";\n";
            }
        } else if (muxName != "OUTPUT" || (isAdder && needsAdderSlicing)) {
            // Skip declaration for adders with slicing (will be declared separately)
            if (!(isAdder && needsAdderSlicing)) {
                muxCode << "\twire" << PrintWireBitWidth(declaredBitWidth) << interWireName << ";\n";
            }
        }
    }
    
    // For shifters, declare the sliced output wire
    if (needsSlicing && !slicedWireName.empty()) {
        muxCode << "\twire" << PrintWireBitWidth(bitwidth) << slicedWireName << ";\n";
    }
    
    // For adders with slicing (single value case), declare OUTPUT_CARRY, OUTPUT, UNSLICED and SLICED wires
    if (needsAdderSlicing && isAdder && nbMuxInputs <= 1) {
        const auto maxOperandBW = inputBitWidth;  // max of left and right operands
        muxCode << "\twire OUTPUT_CARRY;\n";
        muxCode << "\twire" << PrintWireBitWidth(maxOperandBW) << "OUTPUT;\n";
        muxCode << "\twire" << PrintWireBitWidth(declaredBitWidth) << "OUTPUT_UNSLICED;\n";
        muxCode << "\twire" << PrintWireBitWidth(bitwidth) << "OUTPUT_SLICED;\n";
    }

    // Lambda to generate assignment statement
    auto generateAssignment = [&](const int v, const std::string& prefix = "") {
        if (isAdder) {
            const char op = v == -1 ? '-' : '+';
            const std::string leftWire = "LEFT_INPUT_ADDER_" + std::to_string(adderIdx) + "_SHIFTED";
            const std::string rightWire = "RIGHT_INPUT_ADDER_" + std::to_string(adderIdx) + "_SHIFTED";
            if (needsAdderSlicing) {
                muxCode << prefix << "{OUTPUT_CARRY, OUTPUT} = " << leftWire << " " << op << " " << rightWire << ";\n";
            } else {
                muxCode << prefix << interWireName << " = " << leftWire << " " << op << " " << rightWire << ";\n";
            }
        } else if (isInput) {
            if (v < 0) {
                if (v == -1) {
                    muxCode << prefix << interWireName << " = X;\n";
                } else {
                    muxCode << prefix << interWireName << " = X << " << std::abs(v + 1) << ";\n";
                }
            } else {
                muxCode << prefix << interWireName << " = ADDER_" << v << "_SHIFTED_OUTPUT;\n";
            }
        } else if (muxName == "OUTPUT") {
            // Output shift multiplexer - shift the OUTPUT_SLICED wire (or OUTPUT if no slicing)
            if (v == 0) {
                muxCode << prefix << interWireName << " = OUTPUT_SLICED;\n";
            } else {
                muxCode << prefix << interWireName << " = OUTPUT_SLICED << " << (1 << v) << ";\n";
            }
        } else {
            // Shift multiplexers (LEFT/RIGHT_INPUT shifts)
            const std::string sourceWire = muxName.substr(0, muxName.find("_ADDER_"));
            if (v == 0) {
                muxCode << prefix << interWireName << " = " << sourceWire << "_ADDER_" << adderIdx << ";\n";
            } else {
                muxCode << prefix << interWireName << " = " << sourceWire << "_ADDER_" << adderIdx << " << " << v << ";\n";
            }
        }
    };

    // Single value case - use assign
    if (nbPossibleValues == 1) {
        for (const int v : param.possibleValuesFusion) {
            if (solutionNode_.rscm.set.test(paramBitPos + v + param.zeroPoint)) {
                generateAssignment(v, "\tassign ");
                // Add slicing assignment for shifters
                if (needsSlicing && !slicedWireName.empty()) {
                    muxCode << "\tassign " << slicedWireName << " = " << interWireName 
                            << "[" << (bitwidth - 1) << ":0];\n";
                }
                // For adders with slicing, concatenate carry and output, then slice
                if (needsAdderSlicing && isAdder) {
                    muxCode << "\tassign OUTPUT_UNSLICED = {OUTPUT_CARRY, OUTPUT};\n";
                    muxCode << "\tassign OUTPUT_SLICED = OUTPUT_UNSLICED[" << (bitwidth - 1) << ":0];\n";
                }
                return muxCode.str();
            }
        }
    }

    // Multiple values case - use case statement
    muxCode << "\talways @ * begin\n\t\tcase (MUX_INPUT_" << paramGlobalIdx << ")\n";
    
    unsigned int currentSelBit = 0;

    unsigned int firstV = 0;
    for (const int v : param.possibleValuesFusion) {
        if (solutionNode_.rscm.set.test(paramBitPos + v + param.zeroPoint)) {
            if (currentSelBit == 0) firstV = v;
            muxCode << "\t\t\t" << PrintBinaryMuxCode(nbMuxInputs, currentSelBit) << ": ";
            generateAssignment(v);
            ++currentSelBit;
        }
    }

    // if nbMuxInputs not a power of two, need to add a default case
    if (nbMuxInputs == 0 || (nbMuxInputs & nbMuxInputs - 1) != 0)
    {
        // make the zero case the default case
        muxCode << "\t\t\tdefault: ";
        generateAssignment(static_cast<int>(firstV));
    }

    muxCode << "\t\tendcase\n\tend\n";
    
    // Add slicing assignment for shifters
    if (needsSlicing && !slicedWireName.empty())
    {
        muxCode << "\tassign " << slicedWireName << " = " << interWireName 
                << "[" << (bitwidth - 1) << ":0];\n";
    }
    
    // For adders with slicing, concatenate carry and output, then slice
    if (needsAdderSlicing && isAdder) {
        muxCode << "\tassign OUTPUT_UNSLICED = {OUTPUT_CARRY, OUTPUT};\n";
        muxCode << "\tassign OUTPUT_SLICED = OUTPUT_UNSLICED[" << (bitwidth - 1) << ":0];\n";
    }
    
    return muxCode.str();
}

void VerilogGenerator::ProcessParameter(
    const Variables& param,
    size_t& bitPos,
    unsigned int& paramGlobalIdx,
    const unsigned int paramInAdderIdx,
    const unsigned int adderIdx,
    const size_t nbVarsPerAdder,
    std::unordered_map<unsigned int, unsigned int>& muxInputsMap,
    std::stringstream& left_shift_mux,
    std::stringstream& right_shift_mux,
    std::stringstream& outputShiftMux,
    std::stringstream& plusMinus,
    std::stringstream& left_input,
    std::stringstream& right_input,
    const unsigned int leftInputPortBW,
    const unsigned int rightInputPortBW,
    const unsigned int adderOutputBW
) const
{
    const size_t paramBitPos = bitPos;
    
    // Count set bits
    unsigned int bitsSet = 0;
    for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
        if (solutionNode_.rscm.set.test(bitPos + j)) {
            ++bitsSet;
        }
    }
    bitPos += param.possibleValuesFusion.size();

    // Track mux inputs
    if (bitsSet > 1) {
        // Add to adder's mux map only if not input routing
        const VariableDefs varType = idxToVarMap_.at(paramInAdderIdx);
        if (varType != VariableDefs::LEFT_INPUTS && varType != VariableDefs::RIGHT_INPUTS) {
            muxInputsMap[paramGlobalIdx] = bitsSet;
        }
    }

    // Generate code based on variable type
    const VariableDefs varType = idxToVarMap_.at(paramInAdderIdx);
    const unsigned int bitWidth = GetParamBitWidth(varType, adderIdx, nbVarsPerAdder);

    struct MuxConfig {
        VariableDefs type;
        std::string name;
        std::stringstream* stream;
        bool needsInputBW;
    };

    const std::array<MuxConfig, 6> muxConfigs = {{
        {VariableDefs::LEFT_SHIFTS, "LEFT_INPUT_ADDER_" + std::to_string(adderIdx), &left_shift_mux, true},
        {VariableDefs::RIGHT_SHIFTS, "RIGHT_INPUT_ADDER_" + std::to_string(adderIdx), &right_shift_mux, true},
        {VariableDefs::OUTPUTS_SHIFTS, "OUTPUT", &outputShiftMux, true},
        {VariableDefs::RIGHT_MULTIPLIER, "RIGHT_MULTIPLIER", &plusMinus, true},
        {VariableDefs::LEFT_INPUTS, "LEFT_INPUT", &left_input, false},
        {VariableDefs::RIGHT_INPUTS, "RIGHT_INPUT", &right_input, false}
    }};

    for (const auto& [type, name, stream, needsInputBW] : muxConfigs) {
        if (varType == type) {
            unsigned int inputBW = 0;
            if (needsInputBW) {
                // Get the input bitwidth for shifters - use actual port sizes
                if (type == VariableDefs::LEFT_SHIFTS) {
                    inputBW = leftInputPortBW;
                } else if (type == VariableDefs::RIGHT_SHIFTS) {
                    inputBW = rightInputPortBW;
                } else if (type == VariableDefs::OUTPUTS_SHIFTS) {
                    inputBW = adderOutputBW;
                } else if (type == VariableDefs::RIGHT_MULTIPLIER) {
                    // For adder, pass the max of left and right shifted sizes
                    const auto leftShiftedBW = GetParamBitWidth(VariableDefs::LEFT_SHIFTS, adderIdx, nbVarsPerAdder);
                    const auto rightShiftedBW = GetParamBitWidth(VariableDefs::RIGHT_SHIFTS, adderIdx, nbVarsPerAdder);
                    inputBW = std::max(leftShiftedBW, rightShiftedBW);
                }
            }
            *stream << GenerateMuxCode(param, paramBitPos, paramGlobalIdx, bitsSet, name, bitWidth, adderIdx, inputBW, nbVarsPerAdder);
        }
    }

    ++paramGlobalIdx;
}

void VerilogGenerator::PrintModules()
{
    size_t bitPos = 0;
    unsigned int adderIdx = 0;
    unsigned int layerIdx = 0;
    unsigned int paramGlobalIdx = 0;
    const size_t nbVarsPerAdder = layers_[0].adders[0].variables.size();
    std::stringstream rscmTopModule;
    unsigned int rscmOutputBW = 0;
    std::vector<std::pair<std::stringstream, std::stringstream>> adderInputs;
    std::vector<std::unordered_map<unsigned int, unsigned int>> adderMuxMaps;
    std::unordered_map<unsigned int, unsigned int> allMuxInputs;
    std::vector<unsigned int> adderOutputSizes; // Track each adder's output size

    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            std::unordered_map<unsigned int, unsigned int> muxInputsMap;
            std::stringstream left_shift_mux, right_shift_mux, outputShiftMux, plusMinus, left_input, right_input;
            unsigned int paramInAdderIdx = 0;

            // Process all parameters
            // Calculate module port sizes first
            // For layer 0, inputs are from X
            // For other layers, we need to determine from the routing parameters
            auto leftInputBW = inputBW_;
            auto rightInputBW = inputBW_;
            auto outputBW = GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, adderIdx, nbVarsPerAdder);
            rscmOutputBW = outputBW; // store for top module
            adderOutputSizes.push_back(outputBW); // Track this adder's output size

            if (layerIdx == 0)
            {
                leftInputBW = inputBW_;
                rightInputBW = inputBW_;
            }
            else
            {
                // For non-layer-0, use conservative estimates
                // These will be recalculated after routing code generation
                leftInputBW = inputBW_;
                rightInputBW = inputBW_;
                // Check if any previous adder outputs might be used
                for (size_t i = 0; i < adderOutputSizes.size(); ++i) {
                    leftInputBW = std::max(leftInputBW, adderOutputSizes[i]);
                    rightInputBW = std::max(rightInputBW, adderOutputSizes[i]);
                }
            }
            
            // Get adder output size (input to output shifter)
            const auto adderOutputBW = GetParamBitWidth(VariableDefs::RIGHT_MULTIPLIER, adderIdx, nbVarsPerAdder);
            
            // Process all parameters
            for (const auto& param : adder.variables) {
                // Track bits set for all mux inputs (including input routing)
                const size_t currentBitPos = bitPos;
                unsigned int bitsSet = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (solutionNode_.rscm.set.test(currentBitPos + j)) {
                        ++bitsSet;
                    }
                }
                if (bitsSet > 1) {
                    allMuxInputs[paramGlobalIdx] = bitsSet;
                }
                
                ProcessParameter(param, bitPos, paramGlobalIdx, paramInAdderIdx, adderIdx,
                               nbVarsPerAdder, muxInputsMap, left_shift_mux, right_shift_mux,
                               outputShiftMux, plusMinus, left_input, right_input,
                               leftInputBW, rightInputBW, adderOutputBW);
                ++paramInAdderIdx;
            }
            adderInputs.emplace_back(std::move(left_input), std::move(right_input));
            adderMuxMaps.push_back(muxInputsMap);

            // Calculate actual port sizes from routing code
            // For layer 0, inputs are always from X
            if (layerIdx == 0) {
                leftInputBW = inputBW_;
                rightInputBW = inputBW_;
            } else {
                // Calculate based on routing code
                const bool leftHasMux = adderInputs.back().first.str().find("always") != std::string::npos;
                const bool rightHasMux = adderInputs.back().second.str().find("always") != std::string::npos;

                if (leftHasMux) {
                    leftInputBW = inputBW_;
                    const std::string& leftCode = adderInputs.back().first.str();
                    for (int shift = 0; shift <= 5; ++shift) {
                        if (leftCode.find("X << " + std::to_string(shift)) != std::string::npos) {
                            leftInputBW = std::max(leftInputBW, inputBW_ + static_cast<unsigned int>(shift));
                        }
                    }
                    for (size_t i = 0; i < adderOutputSizes.size(); ++i) {
                        if (leftCode.find("ADDER_" + std::to_string(i) + "_SHIFTED_OUTPUT") != std::string::npos) {
                            leftInputBW = std::max(leftInputBW, adderOutputSizes[i]);
                        }
                    }
                } else {
                    // Direct routing - determine from source
                    const std::string& leftCode = adderInputs.back().first.str();
                    // Check if routing from X - need to check ADDER_ after '=' to avoid matching target wire
                    size_t eqPos = leftCode.find('=');
                    size_t adderPosAfterEq = (eqPos != std::string::npos) ? leftCode.find("ADDER_", eqPos) : std::string::npos;
                    if (leftCode.find("= X") != std::string::npos && adderPosAfterEq == std::string::npos) {
                        leftInputBW = inputBW_;
                    } else {
                        // Find ADDER_ after the '=' sign
                        if (adderPosAfterEq != std::string::npos) {
                            for (size_t i = 0; i < adderOutputSizes.size(); ++i) {
                                if (leftCode.find("ADDER_" + std::to_string(i), eqPos) != std::string::npos) {
                                    leftInputBW = adderOutputSizes[i];
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (rightHasMux) {
                    rightInputBW = inputBW_;
                    const std::string& rightCode = adderInputs.back().second.str();
                    for (int shift = 0; shift <= 5; ++shift) {
                        if (rightCode.find("X << " + std::to_string(shift)) != std::string::npos) {
                            rightInputBW = std::max(rightInputBW, inputBW_ + static_cast<unsigned int>(shift));
                        }
                    }
                    for (size_t i = 0; i < adderOutputSizes.size(); ++i) {
                        if (rightCode.find("ADDER_" + std::to_string(i) + "_SHIFTED_OUTPUT") != std::string::npos) {
                            rightInputBW = std::max(rightInputBW, adderOutputSizes[i]);
                        }
                    }
                } else {
                    // Direct routing - determine from source
                    const std::string& rightCode = adderInputs.back().second.str();
                    // Check if routing from X - need to check ADDER_ after '=' to avoid matching target wire
                    size_t eqPos = rightCode.find('=');
                    size_t adderPosAfterEq = (eqPos != std::string::npos) ? rightCode.find("ADDER_", eqPos) : std::string::npos;
                    if (rightCode.find("= X") != std::string::npos && adderPosAfterEq == std::string::npos) {
                        rightInputBW = inputBW_;
                    } else {
                        // Find ADDER_ after the '=' sign
                        if (adderPosAfterEq != std::string::npos) {
                            for (size_t i = 0; i < adderOutputSizes.size(); ++i) {
                                if (rightCode.find("ADDER_" + std::to_string(i), eqPos) != std::string::npos) {
                                    rightInputBW = adderOutputSizes[i];
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // Write module header
            outputFile_ << "module ADDER_" << adderIdx << " (\n"
                       << "\tinput" << PrintWireBitWidth(leftInputBW) << "LEFT_INPUT_ADDER_" << adderIdx << ",\n"
                       << "\tinput" << PrintWireBitWidth(rightInputBW) << "RIGHT_INPUT_ADDER_" << adderIdx << ",\n";

            // Write mux select inputs
            for (const auto& [paramIdx, muxInputCount] : muxInputsMap) {
                const auto muxSelBits = ComputeBitWidth(static_cast<int>(muxInputCount) - 1);
                outputFile_ << "\tinput" << PrintWireBitWidth(muxSelBits) << "MUX_INPUT_" << paramIdx << ",\n";
            }

            outputFile_ << "\toutput" << PrintWireBitWidth(outputBW) << "ADDER_" << adderIdx << "_SHIFTED_OUTPUT\n);\n";

            // Write module body with comments
            struct MuxOutput { const char* comment; const std::stringstream* stream; };
            const std::array<MuxOutput, 4> outputs = {{
                {"Left Shift Multiplexer", &left_shift_mux},
                {"Right Shift Multiplexer", &right_shift_mux},
                {"Adder", &plusMinus},
                {"Output Shift Multiplexer", &outputShiftMux}
            }};

            for (const auto& [comment, stream] : outputs) {
                if (!stream->str().empty()) {
                    outputFile_ << "\t// " << comment << "\n" << stream->str() << "\n";
                }
            }
            
            // Assign module output if no output shift mux (otherwise it's already assigned)
            if (outputShiftMux.str().empty() && !plusMinus.str().empty()) {
                outputFile_ << "\t// Assign output\n";
                // Check if OUTPUT_UNSLICED was generated (adder needed slicing)
                if (plusMinus.str().find("OUTPUT_UNSLICED") != std::string::npos) {
                    outputFile_ << "\tassign ADDER_" << adderIdx << "_SHIFTED_OUTPUT = OUTPUT_SLICED;\n\n";
                } else {
                    outputFile_ << "\tassign ADDER_" << adderIdx << "_SHIFTED_OUTPUT = OUTPUT;\n\n";
                }
            }

            outputFile_ << "endmodule\n\n";
            ++adderIdx;
        }
        ++layerIdx;
    }

    // Build top module header with MUX inputs from collected maps
    rscmTopModule << "module RSCM (\n"
                  << "\tinput" << PrintWireBitWidth(inputBW_) <<" X,\n";
    
    // Add all MUX_INPUT ports (including input routing muxes)
    for (const auto& [paramIdx, muxInputCount] : allMuxInputs) {
        const auto muxSelBits = ComputeBitWidth(static_cast<int>(muxInputCount) - 1);
        rscmTopModule << "\tinput" << PrintWireBitWidth(muxSelBits) << "MUX_INPUT_" << paramIdx << ",\n";
    }
    
    rscmTopModule << "\toutput" << PrintWireBitWidth(rscmOutputBW) << "OUTPUT\n"
                  << ");\n\n";


    // Instantiate each adder with its routing
    unsigned int adderCounter = 0;
    for (const auto& [leftInputStream, rightInputStream] : adderInputs)
    {
        const auto outputBW = GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, adderCounter, nbVarsPerAdder);
        rscmTopModule << "\twire" << PrintWireBitWidth(outputBW) << "ADDER_" << adderCounter << "_SHIFTED_OUTPUT;\n";
        // Determine if inputs need wire declarations (for muxed inputs in top module)
        const bool leftHasMux = leftInputStream.str().find("always") != std::string::npos;
        const bool rightHasMux = rightInputStream.str().find("always") != std::string::npos;
        
        // Determine actual bitwidth for routing wires by checking what's being routed
        auto getRoutingBitWidth = [&](const std::string& streamStr) -> unsigned int {
            // Check if routing from X (input), possibly with shift
            if (streamStr.find("= X") != std::string::npos && streamStr.find("ADDER_") == std::string::npos) {
                // Check for shift
                size_t shiftPos = streamStr.find("<< ");
                if (shiftPos != std::string::npos) {
                    // Extract shift amount
                    size_t numStart = shiftPos + 3;
                    size_t numEnd = streamStr.find(';', numStart);
                    if (numEnd != std::string::npos) {
                        std::string shiftStr = streamStr.substr(numStart, numEnd - numStart);
                        // Remove any whitespace
                        shiftStr.erase(std::remove_if(shiftStr.begin(), shiftStr.end(), ::isspace), shiftStr.end());
                        int shiftAmount = std::stoi(shiftStr);
                        return inputBW_ + shiftAmount;
                    }
                }
                return inputBW_;
            }
            // Check if routing from another adder's output
            size_t adderPos = streamStr.find("ADDER_");
            // Skip if this is the LEFT_INPUT_ADDER or RIGHT_INPUT_ADDER (the target, not source)
            if (adderPos != std::string::npos) {
                // Check if this is INPUT_ADDER (target wire), not source
                if (adderPos >= 6) {
                    std::string before = streamStr.substr(adderPos - 6, 6);
                    if (before == "INPUT_") {
                        // This is the target wire, find the source after the '='
                        size_t eqPos = streamStr.find('=', adderPos);
                        if (eqPos != std::string::npos) {
                            adderPos = streamStr.find("ADDER_", eqPos);
                        }
                    }
                }
            }
            if (adderPos != std::string::npos) {
                // Extract adder index number
                size_t numStart = adderPos + 6; // length of "ADDER_"
                size_t numEnd = streamStr.find('_', numStart);
                if (numEnd != std::string::npos) {
                    int sourceAdderIdx = std::stoi(streamStr.substr(numStart, numEnd - numStart));
                    if (sourceAdderIdx >= 0 && sourceAdderIdx < static_cast<int>(adderOutputSizes.size())) {
                        return adderOutputSizes[sourceAdderIdx];
                    }
                }
            }
            // Default fallback
            return inputBW_;
        };

        if (leftHasMux) {
            // Calculate actual bitwidth needed from mux routing code
            unsigned int maxLeftBW = inputBW_;
            const std::string& leftCode = leftInputStream.str();
            
            // Check for shifts from X
            for (int shift = 0; shift <= 5; ++shift) {
                std::string shiftPattern = "X << " + std::to_string(shift);
                if (leftCode.find(shiftPattern) != std::string::npos) {
                    maxLeftBW = std::max(maxLeftBW, inputBW_ + static_cast<unsigned int>(shift));
                }
            }
            // Check for plain X
            if (leftCode.find("= X;") != std::string::npos) {
                maxLeftBW = std::max(maxLeftBW, inputBW_);
            }
            // Check for routing from previous adders
            for (size_t i = 0; i < adderOutputSizes.size(); ++i) {
                std::string adderPattern = "ADDER_" + std::to_string(i) + "_SHIFTED_OUTPUT";
                if (leftCode.find(adderPattern) != std::string::npos) {
                    maxLeftBW = std::max(maxLeftBW, adderOutputSizes[i]);
                }
            }
            rscmTopModule << "\treg" << PrintWireBitWidth(maxLeftBW) << "LEFT_INPUT_ADDER_" << adderCounter << ";\n";
        } else
        {
            const auto leftRoutingBW = getRoutingBitWidth(leftInputStream.str());
            rscmTopModule << "\twire" << PrintWireBitWidth(leftRoutingBW) << "LEFT_INPUT_ADDER_" << adderCounter << " ;\n";
        }
        if (rightHasMux) {
            // Calculate actual bitwidth needed from mux routing code
            unsigned int maxRightBW = inputBW_;
            const std::string& rightCode = rightInputStream.str();
            
            // Check for shifts from X
            for (int shift = 0; shift <= 5; ++shift) {
                std::string shiftPattern = "X << " + std::to_string(shift);
                if (rightCode.find(shiftPattern) != std::string::npos) {
                    maxRightBW = std::max(maxRightBW, inputBW_ + static_cast<unsigned int>(shift));
                }
            }
            // Check for plain X
            if (rightCode.find("= X;") != std::string::npos) {
                maxRightBW = std::max(maxRightBW, inputBW_);
            }
            // Check for routing from previous adders
            for (size_t i = 0; i < adderOutputSizes.size(); ++i) {
                std::string adderPattern = "ADDER_" + std::to_string(i) + "_SHIFTED_OUTPUT";
                if (rightCode.find(adderPattern) != std::string::npos) {
                    maxRightBW = std::max(maxRightBW, adderOutputSizes[i]);
                }
            }
            rscmTopModule << "\treg" << PrintWireBitWidth(maxRightBW) << "RIGHT_INPUT_ADDER_" << adderCounter << ";\n";
        } else
        {
            const auto rightRoutingBW = getRoutingBitWidth(rightInputStream.str());
            rscmTopModule << "\twire" << PrintWireBitWidth(rightRoutingBW) << "RIGHT_INPUT_ADDER_" << adderCounter << " ;\n";
        }
        
        rscmTopModule << "\t// Adder " << adderCounter << " Input Routing\n";
        rscmTopModule << leftInputStream.str() << rightInputStream.str();
        
        rscmTopModule << "\t// Adder " << adderCounter << " Instantiation\n";
        rscmTopModule << "\tADDER_" << adderCounter << " adder_" << adderCounter << "_inst (\n";
        rscmTopModule << "\t\t.LEFT_INPUT_ADDER_" << adderCounter << "(LEFT_INPUT_ADDER_" << adderCounter << "),\n";
        rscmTopModule << "\t\t.RIGHT_INPUT_ADDER_" << adderCounter << "(RIGHT_INPUT_ADDER_" << adderCounter << "),\n";
        
        // Connect MUX inputs for this adder
        for (const auto& paramIdx : adderMuxMaps[adderCounter] | std::views::keys) {
            rscmTopModule << "\t\t.MUX_INPUT_" << paramIdx << "(MUX_INPUT_" << paramIdx << "),\n";
        }
        
        rscmTopModule << "\t\t.ADDER_" << adderCounter << "_SHIFTED_OUTPUT(ADDER_" << adderCounter << "_SHIFTED_OUTPUT)\n";
        rscmTopModule << "\t);\n\n";
        ++adderCounter;
    }

    // Connect final output
    rscmTopModule << "\t// Final Output\n";
    rscmTopModule << "\tassign OUTPUT = ADDER_" << (adderInputs.size() - 1) << "_SHIFTED_OUTPUT;\n\n";
    rscmTopModule << "endmodule\n";
    outputFile_ << rscmTopModule.str();
}

void VerilogGenerator::PrintTestbench()
{
    if (targets_.empty() || scmDesigns_.empty()) {
        return; // No targets to test
    }

    outputFile_ << "\n\n// Testbench\n";
    outputFile_ << "module RSCM_tb;\n\n";
    
    // Declare testbench signals
    outputFile_ << "\treg [" << (inputBW_ - 1) << ":0] X;\n";
    
    // First pass: compute mux configurations per target (replicate SolveConfigToMuxMapping logic)
    std::vector<std::vector<unsigned int>> muxEncoding(targets_.size());
    std::vector<unsigned int> muxNames;
    std::unordered_map<unsigned int, unsigned int> allMuxInputs;
    
    size_t bitIndex = 0;
    unsigned int parameterIndex = 0;
    
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            for (const auto& parameter : adder.variables) {
                unsigned int nbPossibleValues = 0;
                // Get number of possible values for this parameter in merged RSCM
                for (size_t v = 0; v < parameter.possibleValuesFusion.size(); v++) {
                    if (solutionNode_.rscm.set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                        nbPossibleValues++;
                    }
                }
                
                // If more than one bit is set then we have a multiplexer
                if (nbPossibleValues > 1) {
                    unsigned int valueNumber = 0;
                    muxNames.push_back(parameterIndex);
                    allMuxInputs[parameterIndex] = nbPossibleValues;
                    
                    // For each possible value that's active in the merged RSCM
                    for (size_t v = 0; v < parameter.possibleValuesFusion.size(); v++) {
                        if (solutionNode_.rscm.set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                            // Check which targets use this specific value by looking at their individual SCMs
                            for (size_t targetIndex = 0; targetIndex < targets_.size(); targetIndex++) {
                                const auto scmIndex = solutionNode_.scmIndexes[targetIndex];
                                const auto& targetSCM = scmDesigns_[targetIndex].second[scmIndex];
                                
                                if (targetSCM.set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                                    muxEncoding[targetIndex].push_back(valueNumber);
                                }
                            }
                            valueNumber++;
                        }
                    }
                }
                bitIndex += parameter.possibleValuesFusion.size();
                parameterIndex++;
            }
        }
    }
    
    // Declare MUX_INPUT signals
    for (const auto& [paramIdx, muxInputCount] : allMuxInputs) {
        const auto muxSelBits = ComputeBitWidth(static_cast<int>(muxInputCount) - 1);
        outputFile_ << "\treg [" << (muxSelBits - 1) << ":0] MUX_INPUT_" << paramIdx << ";\n";
    }
    
    const auto outputBW =
        GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, static_cast<unsigned int>(layers_.back().adders.size()) - 1,
                         layers_[0].adders[0].variables.size());
    outputFile_ << "\twire" << PrintWireBitWidth(outputBW) << "OUTPUT;\n\n";
    
    // Instantiate DUT
    outputFile_ << "\t// Instantiate RSCM module\n";
    outputFile_ << "\tRSCM dut (\n";
    outputFile_ << "\t\t.X(X),\n";
    for (const auto& paramIdx : allMuxInputs | std::views::keys) {
        outputFile_ << "\t\t.MUX_INPUT_" << paramIdx << "(MUX_INPUT_" << paramIdx << "),\n";
    }
    outputFile_ << "\t\t.OUTPUT(OUTPUT)\n";
    outputFile_ << "\t);\n\n";
    
    // Test stimulus
    outputFile_ << "\tinitial begin\n";
    outputFile_ << "\t\t$display(\"Starting RSCM Testbench\");\n";
    outputFile_ << "\t\t$display(\"Testing %0d target configurations\", " << targets_.size() << ");\n\n";
    
    // Generate test cases for each target constant
    for (size_t targetIdx = 0; targetIdx < targets_.size(); ++targetIdx) {
        const int target = targets_[targetIdx];
        outputFile_ << "\t\t// ========================================\n";
        outputFile_ << "\t\t// Test configuration for target constant: " << target << "\n";
        outputFile_ << "\t\t// ========================================\n";
        
        // Set MUX inputs based on the configuration for this target
        // Use muxEncoding[targetIdx] which contains the mux select values
        size_t muxIdx = 0;
        for (const auto& muxParamIdx : muxNames) {
            if (muxIdx < muxEncoding[targetIdx].size()) {
                outputFile_ << "\t\tMUX_INPUT_" << muxParamIdx << " = " << muxEncoding[targetIdx][muxIdx] << ";\n";
            }
            muxIdx++;
        }
        
        outputFile_ << "\t\t#1; // Allow mux changes to propagate\n\n";
        
        // Test with various input values
        const std::vector<int> testInputs = {1, 2, 3, 5, 7, 11, 13};
        for (const int testInput : testInputs) {
            if (testInput >= 1 << inputBW_) continue; // Skip if input too large
            
            const long long expected = static_cast<long long>(testInput) * target;
            
            outputFile_ << "\t\tX = " << testInput << ";\n";
            outputFile_ << "\t\t#10;\n";
            outputFile_ << "\t\tif (OUTPUT !== " << expected << ") begin\n";
            outputFile_ << "\t\t\t$display(\"ERROR: Target=%0d, X=%0d, Expected=%0d, Got=%0d\", "
                       << target << ", X, " << expected << ", OUTPUT);\n";
            outputFile_ << "\t\tend else begin\n";
            outputFile_ << "\t\t\t$display(\"PASS: Target=%0d, X=%0d, Output=%0d\", "
                       << target << ", X, OUTPUT);\n";
            outputFile_ << "\t\tend\n\n";
        }
    }
    
    outputFile_ << "\t\t$display(\"========================================\");\n";
    outputFile_ << "\t\t$display(\"Testbench completed\");\n";
    outputFile_ << "\t\t$finish;\n";
    outputFile_ << "\tend\n\n";
    outputFile_ << "endmodule\n";
}
