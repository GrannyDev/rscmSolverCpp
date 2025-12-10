//
// Created by smith on 09/12/2025.
//

#include "Verilog.h"

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
    PrintTestbench();
}

unsigned int VerilogGenerator::ComputeBitWidth(const int maxVal, const bool is_signed, const bool is_mux) const
{
    unsigned int bw;
    if (!is_signed) bw = static_cast<unsigned int>(std::ceil(std::log2(maxVal + 1)));
    else bw = static_cast<unsigned int>(std::ceil(std::log2(std::max(std::abs(maxVal + 1), std::abs(maxVal))))) + 1;

    if (!is_mux) bw += inputBW_;
    return bw;
}

unsigned int VerilogGenerator::GetMaxOutputValue(
    const VariableDefs varType,
    const unsigned int adderIdx,
    const size_t nbVarsPerAdder) const
{
    const auto baseIdx = varToIdxMap_.at(varType) + adderIdx * nbVarsPerAdder;
    return solutionNode_.rscm.maxOutputValue[baseIdx];
}

std::string VerilogGenerator::GenerateMuxCode(
    const Variables& param,
    const size_t paramBitPos,
    const unsigned int paramGlobalIdx,
    const unsigned int nbPossibleValues,
    const std::string& muxName,
    const unsigned int bitwidth,
    const unsigned int adderIdx
    ) const
{
    std::stringstream muxCode;
    const bool isAdder = muxName == "RIGHT_MULTIPLIER";
    const bool isInput = muxName == "LEFT_INPUT" || muxName == "RIGHT_INPUT";
    
    // Determine wire name based on context
    std::string interWireName;
    std::string outputWireName;
    if (isAdder) {
        interWireName = "OUTPUT";
        outputWireName = "ADDER_" + std::to_string(adderIdx) + "_OUTPUT";
    } else if (isInput) {
        interWireName = muxName + "_ADDER_" + std::to_string(adderIdx);
    } else if (muxName == "OUTPUT") {
        // Output shift multiplexer
        interWireName = "ADDER_" + std::to_string(adderIdx) + "_OUTPUT";
    } else {
        // Shift multiplexers - muxName already contains full name like LEFT_INPUT_ADDER_X
        interWireName = muxName + "_SHIFTED";
    }

    // Declare wire (except for inputs which are ports, and output port)
    const bool isOutputPort = (muxName == "OUTPUT");
    if (!isInput && !isOutputPort) {
        muxCode << "\twire [" << (bitwidth - 1) << ":0] " << interWireName << ";\n";
    }

    // Lambda to generate assignment statement
    auto generateAssignment = [&](const int v, const std::string& prefix = "") {
        if (isAdder) {
            const char op = v == -1 ? '-' : '+';
            const std::string leftWire = "LEFT_INPUT_ADDER_" + std::to_string(adderIdx) + "_SHIFTED";
            const std::string rightWire = "RIGHT_INPUT_ADDER_" + std::to_string(adderIdx) + "_SHIFTED";
            muxCode << prefix << interWireName << " = " << leftWire << " " << op << " " << rightWire << ";\n";
        } else if (isInput) {
            if (v < 0) {
                if (v == -1) {
                    muxCode << prefix << interWireName << " = X;\n";
                } else {
                    muxCode << prefix << interWireName << " = X << " << std::abs(v + 1) << ";\n";
                }
            } else {
                muxCode << prefix << interWireName << " = ADDER_" << v << "_OUTPUT;\n";
            }
        } else if (muxName == "OUTPUT") {
            // Output shift multiplexer - shift the OUTPUT wire
            if (v == 0) {
                muxCode << prefix << interWireName << " = OUTPUT;\n";
            } else {
                muxCode << prefix << interWireName << " = OUTPUT << " << (1 << v) << ";\n";
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
                return muxCode.str();
            }
        }
    }

    // Multiple values case - use case statement
    muxCode << "\talways @ *\n\t\tbegin\n\t\t\tcase (MUX_INPUT_" << paramGlobalIdx << ")\n";
    
    unsigned int currentSelBit = 0;
    for (const int v : param.possibleValuesFusion) {
        if (solutionNode_.rscm.set.test(paramBitPos + v + param.zeroPoint)) {
            muxCode << "\t\t\t\t" << currentSelBit << ": ";
            generateAssignment(v);
            ++currentSelBit;
        }
    }
    
    muxCode << "\t\t\tendcase\n\t\tend\n";
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
    std::stringstream& right_input
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
    const unsigned int bitWidth = ComputeBitWidth(
        GetMaxOutputValue(varType, adderIdx, nbVarsPerAdder), true
    );

    struct MuxConfig {
        VariableDefs type;
        std::string name;
        std::stringstream* stream;
    };

    const std::array<MuxConfig, 6> muxConfigs = {{
        {VariableDefs::LEFT_SHIFTS, "LEFT_INPUT_ADDER_" + std::to_string(adderIdx), &left_shift_mux},
        {VariableDefs::RIGHT_SHIFTS, "RIGHT_INPUT_ADDER_" + std::to_string(adderIdx), &right_shift_mux},
        {VariableDefs::OUTPUTS_SHIFTS, "OUTPUT", &outputShiftMux},
        {VariableDefs::RIGHT_MULTIPLIER, "RIGHT_MULTIPLIER", &plusMinus},
        {VariableDefs::LEFT_INPUTS, "LEFT_INPUT", &left_input},
        {VariableDefs::RIGHT_INPUTS, "RIGHT_INPUT", &right_input}
    }};

    for (const auto& [type, name, stream] : muxConfigs) {
        if (varType == type) {
            *stream << GenerateMuxCode(param, paramBitPos, paramGlobalIdx, bitsSet, name, bitWidth, adderIdx);
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

    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            std::unordered_map<unsigned int, unsigned int> muxInputsMap;
            std::stringstream left_shift_mux, right_shift_mux, outputShiftMux, plusMinus, left_input, right_input;
            unsigned int paramInAdderIdx = 0;

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
                               outputShiftMux, plusMinus, left_input, right_input);
                ++paramInAdderIdx;
            }
            adderInputs.emplace_back(std::move(left_input), std::move(right_input));
            adderMuxMaps.push_back(muxInputsMap);

            // Write module header
            auto leftInputBW = ComputeBitWidth(GetMaxOutputValue(VariableDefs::LEFT_INPUTS, adderIdx, nbVarsPerAdder), true);
            auto rightInputBW = ComputeBitWidth(GetMaxOutputValue(VariableDefs::RIGHT_INPUTS, adderIdx, nbVarsPerAdder), true);
            const auto outputBW = ComputeBitWidth(GetMaxOutputValue(VariableDefs::OUTPUTS_SHIFTS, adderIdx, nbVarsPerAdder), true);
            rscmOutputBW = outputBW; // store for top module

            // inputs of layer 0 adders are always -1, needs to be handled
            if (layerIdx == 0)
            {
                leftInputBW = inputBW_;
                rightInputBW = inputBW_;
            }

            outputFile_ << "module ADDER_" << adderIdx << " (\n"
                       << "\tinput [" << (leftInputBW - 1) << ":0] LEFT_INPUT_ADDER_" << adderIdx << ",\n"
                       << "\tinput [" << (rightInputBW - 1) << ":0] RIGHT_INPUT_ADDER_" << adderIdx << ",\n";

            // Write mux select inputs
            for (const auto& [paramIdx, muxInputCount] : muxInputsMap) {
                const auto muxSelBits = ComputeBitWidth(static_cast<int>(muxInputCount) - 1, false, true);
                outputFile_ << "\tinput [" << (muxSelBits - 1) << ":0] MUX_INPUT_" << paramIdx << ",\n";
            }

            outputFile_ << "\toutput [" << (outputBW - 1) << ":0] ADDER_" << adderIdx << "_OUTPUT\n);\n";

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
                outputFile_ << "\t// Assign output\n";\
                outputFile_ << "\tassign ADDER_" << adderIdx << "_OUTPUT = OUTPUT;\n\n";
            }

            outputFile_ << "endmodule\n\n";
            ++adderIdx;
        }
        ++layerIdx;
    }

    // Build top module header with MUX inputs from collected maps
    rscmTopModule << "module RSCM (\n"
                  << "\tinput [" << (inputBW_ - 1) << ":0] X,\n";
    
    // Add all MUX_INPUT ports (including input routing muxes)
    for (const auto& [paramIdx, muxInputCount] : allMuxInputs) {
        const auto muxSelBits = ComputeBitWidth(static_cast<int>(muxInputCount) - 1, false, true);
        rscmTopModule << "\tinput [" << (muxSelBits - 1) << ":0] MUX_INPUT_" << paramIdx << ",\n";
    }
    
    rscmTopModule << "\toutput [" << (rscmOutputBW - 1) << ":0] OUTPUT\n"
                  << ");\n\n";

    // Declare output wires for each adder
    for (unsigned int i = 0; i < adderInputs.size(); ++i) {
        const auto outputBW = ComputeBitWidth(GetMaxOutputValue(VariableDefs::OUTPUTS_SHIFTS, i, nbVarsPerAdder), true);
        rscmTopModule << "\twire [" << (outputBW - 1) << ":0] ADDER_" << i << "_OUTPUT;\n";
    }
    rscmTopModule << "\n";

    // Instantiate each adder with its routing
    unsigned int adderCounter = 0;
    for (const auto& [leftInputStream, rightInputStream] : adderInputs)
    {
        // Determine if inputs need wire declarations (for muxed inputs in top module)
        const bool leftHasMux = leftInputStream.str().find("always") != std::string::npos;
        const bool rightHasMux = rightInputStream.str().find("always") != std::string::npos;
        
        if (leftHasMux || rightHasMux) {
            rscmTopModule << "\t// Adder " << adderCounter << " Input Wires\n";
            if (leftHasMux) {
                const auto leftBW = ComputeBitWidth(GetMaxOutputValue(VariableDefs::LEFT_INPUTS, adderCounter, nbVarsPerAdder), true);
                rscmTopModule << "\twire [" << (leftBW - 1) << ":0] LEFT_INPUT_ADDER_" << adderCounter << ";\n";
            }
            if (rightHasMux) {
                const auto rightBW = ComputeBitWidth(GetMaxOutputValue(VariableDefs::RIGHT_INPUTS, adderCounter, nbVarsPerAdder), true);
                rscmTopModule << "\twire [" << (rightBW - 1) << ":0] RIGHT_INPUT_ADDER_" << adderCounter << ";\n";
            }
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
        
        rscmTopModule << "\t\t.ADDER_" << adderCounter << "_OUTPUT(ADDER_" << adderCounter << "_OUTPUT)\n";
        rscmTopModule << "\t);\n\n";
        ++adderCounter;
    }

    // Connect final output
    rscmTopModule << "\t// Final Output\n";
    rscmTopModule << "\tassign OUTPUT = ADDER_" << (adderInputs.size() - 1) << "_OUTPUT;\n\n";
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
        const auto muxSelBits = ComputeBitWidth(static_cast<int>(muxInputCount) - 1, false, true);
        outputFile_ << "\treg [" << (muxSelBits - 1) << ":0] MUX_INPUT_" << paramIdx << ";\n";
    }
    
    const auto outputBW = ComputeBitWidth(
        GetMaxOutputValue(VariableDefs::OUTPUTS_SHIFTS, static_cast<unsigned int>(layers_.back().adders.size()) - 1,
                         layers_[0].adders[0].variables.size()), true);
    outputFile_ << "\twire [" << (outputBW - 1) << ":0] OUTPUT;\n\n";
    
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
