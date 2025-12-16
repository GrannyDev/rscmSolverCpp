//
// Created by smith on 16/12/2025.
//

#include "JSONDumper.h"
# include <filesystem>

JSONDumper::JSONDumper(
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

    Dump();
}

void JSONDumper::Dump() {
    unsigned int nbTabs = 0;
    outputFile_ << "{\n";
    nbTabs++;

    // adders
    unsigned int adderNb = 0;
    unsigned int bitPos = 0;
    
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            outputFile_ << std::string(nbTabs, '\t') << "\"Adder_" << adderNb << "\": {\n";
            nbTabs++;

            // Get parameter indices for this adder
            unsigned int paramGlobalIdxOffset = 0;
            unsigned int addersProcessed = 0;
            for (const auto& l : layers_) {
                for (const auto& a : l.adders) {
                    if (addersProcessed >= adderNb) {
                        break;
                    }
                    paramGlobalIdxOffset += a.variables.size();
                    addersProcessed++;
                }
                if (addersProcessed >= adderNb) {
                    break;
                }
            }

            // Process each parameter
            unsigned int paramIdx = 0;
            for (const auto& param : adder.variables) {
                unsigned int paramGlobalIdx = paramGlobalIdxOffset + paramIdx;
                VariableDefs currentParamType = idxToVarMap_.at(paramIdx);
                const unsigned int computedBW = solutionNode_.variableBitWidths[paramGlobalIdx];

                // Skip non-shift and non-adder parameters
                if (currentParamType == VariableDefs::LEFT_INPUTS || 
                    currentParamType == VariableDefs::RIGHT_INPUTS) {
                    paramIdx++;
                    bitPos += param.possibleValuesFusion.size();
                    continue;
                }

                // Collect selected shifts and find max shift
                std::vector<int> selectedShifts;
                unsigned int maxShift = 0;
                for (const int v : param.possibleValuesFusion) {
                    if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
                        selectedShifts.push_back(v);
                        if (currentParamType == VariableDefs::OUTPUTS_SHIFTS) {
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

                // Check if this is a mux (more than one option selected)
                bool isMux = selectedShifts.size() > 1;
                if (currentParamType == VariableDefs::RIGHT_MULTIPLIER) {
                    // For ADD_SUB, it's a mux if it can do both add and subtract
                    int addCount = 0, subCount = 0;
                    for (const int v : selectedShifts) {
                        if (v == -1) subCount++;
                        else if (v == 1) addCount++;
                    }
                    isMux = (addCount > 0 && subCount > 0);
                }

                // Get theoretical bitwidth
                unsigned int theoreticalBW = 0;
                unsigned int inputBW = 0;
                if (currentParamType == VariableDefs::LEFT_SHIFTS) {
                    inputBW = GetParamBitWidth(VariableDefs::LEFT_INPUTS, adderNb);
                    theoreticalBW = inputBW + maxShift;
                } else if (currentParamType == VariableDefs::RIGHT_SHIFTS) {
                    inputBW = GetParamBitWidth(VariableDefs::RIGHT_INPUTS, adderNb);
                    theoreticalBW = inputBW + maxShift;
                } else if (currentParamType == VariableDefs::OUTPUTS_SHIFTS) {
                    inputBW = GetParamBitWidth(VariableDefs::RIGHT_MULTIPLIER, adderNb);
                    theoreticalBW = inputBW + maxShift;
                } else if (currentParamType == VariableDefs::RIGHT_MULTIPLIER) {
                    // For ADD_SUB, theoretical is max of left and right shift bitwidths + 1
                    unsigned int leftShiftBW = GetParamBitWidth(VariableDefs::LEFT_SHIFTS, adderNb);
                    unsigned int rightShiftBW = GetParamBitWidth(VariableDefs::RIGHT_SHIFTS, adderNb);
                    theoreticalBW = std::max(leftShiftBW, rightShiftBW) + 1;
                }

                // Output parameter information
                outputFile_ << std::string(nbTabs, '\t') << "\"" << varDefsTostring_(currentParamType) << "\": {\n";
                nbTabs++;

                if (currentParamType == VariableDefs::RIGHT_MULTIPLIER) {
                    // Handle ADD_SUB specially
                    outputFile_ << std::string(nbTabs, '\t') << "\"type\": \"";
                    
                    // Determine if it's add, subtract, or both
                    int addCount = 0, subCount = 0;
                    for (const int v : selectedShifts) {
                        if (v == -1) subCount++;
                        else if (v == 1) addCount++;
                    }
                    
                    if (addCount > 0 && subCount > 0) {
                        outputFile_ << "both";
                    } else if (subCount > 0) {
                        outputFile_ << "subtractor";
                    } else {
                        outputFile_ << "adder";
                    }
                    outputFile_ << "\",\n";
                    
                    outputFile_ << std::string(nbTabs, '\t') << "\"theoretical_bitwidth\": " << theoreticalBW << ",\n";
                    outputFile_ << std::string(nbTabs, '\t') << "\"computed_bitwidth\": " << computedBW;
                    if (isMux) {
                        outputFile_ << ",\n" << std::string(nbTabs, '\t') << "\"mux_nb\": " << paramIdx;
                    }
                    outputFile_ << "\n";
                } else {
                    // Handle shift parameters
                    outputFile_ << std::string(nbTabs, '\t') << "\"shifts\": [";
                    
                    bool first = true;
                    for (const int v : selectedShifts) {
                        if (!first) outputFile_ << ", ";
                        outputFile_ << v;
                        first = false;
                    }
                    outputFile_ << "],\n";
                    
                    outputFile_ << std::string(nbTabs, '\t') << "\"theoretical_wire_length\": " << theoreticalBW << ",\n";
                    outputFile_ << std::string(nbTabs, '\t') << "\"computed_wire_length\": " << computedBW;
                    if (isMux) {
                        outputFile_ << ",\n" << std::string(nbTabs, '\t') << "\"mux_nb\": " << paramIdx;
                    }
                    outputFile_ << "\n";
                }

                nbTabs--;
                outputFile_ << std::string(nbTabs, '\t') << "}";

                // Add comma if not the last parameter
                paramIdx++;
                bitPos += param.possibleValuesFusion.size();
                if (paramIdx < adder.variables.size()) {
                    outputFile_ << ",";
                }
                outputFile_ << "\n";
            }

            nbTabs--;
            outputFile_ << std::string(nbTabs, '\t') << "},\n";
            adderNb++;
        }
    }

    outputFile_ << std::string(nbTabs, '\t') << "\"RSCM\": {\n";
    nbTabs++;
    
    // Reset for second pass to generate RSCM top module
    adderNb = 0;
    bitPos = 0;
    
    // First, collect mux encoding information and track mux assignments
    // This logic is copied from Solver::SolveConfigToMuxMapping()
    std::vector<std::vector<unsigned int>> muxEncoding(targets_.size());
    unsigned int bitIndex = 0;
    unsigned int currentAdder = 0;
    
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            unsigned int paramIdx = 0;
            for (const auto& parameter : adder.variables) {
                // Count possible values for this parameter
                unsigned int nbPossibleValues = 0;
                for (size_t v = 0; v < parameter.possibleValuesFusion.size(); v++) {
                    if (solutionNode_.rscm.set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                        nbPossibleValues++;
                    }
                }
                
                // If more than one value, we have a mux
                if (nbPossibleValues > 1) {
                    unsigned int valueNumber = 0;
                    for (size_t v = 0; v < parameter.possibleValuesFusion.size(); v++) {
                        if (solutionNode_.rscm.set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                            for (size_t targetIndex = 0; targetIndex < targets_.size(); targetIndex++) {
                                if (scmDesigns_[targetIndex].second[solutionNode_.scmIndexes[targetIndex]].set.test(
                                    bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                                    muxEncoding[targetIndex].push_back(valueNumber);
                                }
                            }
                            valueNumber++;
                        }
                    }
                }
                bitIndex += parameter.possibleValuesFusion.size();
                paramIdx++;
            }
            currentAdder++;
        }
    }
    
    // Now generate the adder information
    adderNb = 0;
    bitPos = 0;

    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            outputFile_ << std::string(nbTabs, '\t') << "\"ADDER_" << adderNb << "\": {\n";
            nbTabs++;
            
            // Get the left and right input information
            std::vector<std::string> leftInputs;
            std::vector<unsigned int> leftShifts;
            std::vector<std::string> rightInputs;
            std::vector<unsigned int> rightShifts;
            unsigned int maxShiftLeft = 0;
            unsigned int maxShiftRight = 0;
            int leftInputMuxNb = -1;
            int rightInputMuxNb = -1;
            
            unsigned int paramIdx = 0;
            for (const auto& param : adder.variables) {
                VariableDefs paramType = idxToVarMap_.at(paramIdx);
                
                if (paramType == VariableDefs::LEFT_INPUTS || paramType == VariableDefs::RIGHT_INPUTS) {
                    std::vector<std::string>& inputs = (paramType == VariableDefs::LEFT_INPUTS) ? leftInputs : rightInputs;
                    std::vector<unsigned int>& shifts = (paramType == VariableDefs::LEFT_INPUTS) ? leftShifts : rightShifts;
                    unsigned int& maxShift = (paramType == VariableDefs::LEFT_INPUTS) ? maxShiftLeft : maxShiftRight;
                    // Check if this parameter is a mux
                    unsigned int bitsSet = 0;
                    int& inputMuxNb = (paramType == VariableDefs::LEFT_INPUTS) ? leftInputMuxNb : rightInputMuxNb;
                    
                    // Collect inputs with embedded shifts
                    for (const int v : param.possibleValuesFusion) {
                        if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
                            bitsSet++;
                            if (v < 0) {
                                if (v == -1) {
                                    inputs.emplace_back("X");
                                    shifts.push_back(0);
                                } else {
                                    unsigned int shiftAmount = std::abs(v + 1);
                                    inputs.emplace_back("X");
                                    shifts.push_back(shiftAmount);
                                    maxShift = std::max(maxShift, shiftAmount);
                                }
                            } else {
                                inputs.push_back("ADDER_" + std::to_string(v));
                                shifts.push_back(0);
                            }
                        }
                    }
                    maxShift += inputBW_;
                    if (bitsSet > 1) {
                        inputMuxNb = static_cast<int>(paramIdx);
                    }
                }
                
                bitPos += param.possibleValuesFusion.size();
                paramIdx++;
            }
            
            // Output left inputs
            outputFile_ << std::string(nbTabs, '\t') << "\"left_inputs\": [";
            for (size_t i = 0; i < leftInputs.size(); ++i) {
                if (i > 0) outputFile_ << ", ";
                outputFile_ << "\"" << leftInputs[i] << "\"";
            }
            outputFile_ << "],\n";
            
            outputFile_ << std::string(nbTabs, '\t') << "\"left_shifts\": [";
            for (size_t i = 0; i < leftShifts.size(); ++i) {
                if (i > 0) outputFile_ << ", ";
                outputFile_ << leftShifts[i];
            }
            outputFile_ << "],\n";
            
            // Add mux_nb if there are multiple left inputs
            if (leftInputMuxNb >= 0) {
                outputFile_ << std::string(nbTabs, '\t') << "\"left_mux_nb\": " << leftInputMuxNb << ",\n";
            }
            
            // Output right inputs
            outputFile_ << std::string(nbTabs, '\t') << "\"right_inputs\": [";
            for (size_t i = 0; i < rightInputs.size(); ++i) {
                if (i > 0) outputFile_ << ", ";
                outputFile_ << "\"" << rightInputs[i] << "\"";
            }
            outputFile_ << "],\n";
            
            outputFile_ << std::string(nbTabs, '\t') << "\"right_shifts\": [";
            for (size_t i = 0; i < rightShifts.size(); ++i) {
                if (i > 0) outputFile_ << ", ";
                outputFile_ << rightShifts[i];
            }
            outputFile_ << "],\n";
            
            // Add mux_nb if there are multiple right inputs
            if (rightInputMuxNb >= 0) {
                outputFile_ << std::string(nbTabs, '\t') << "\"right_mux_nb\": " << rightInputMuxNb << ",\n";
            }
            
            // Compute theoretical and computed bitwidths for left and right
            unsigned int leftComputedBW = GetParamBitWidth(VariableDefs::LEFT_INPUTS, adderNb);
            unsigned int rightComputedBW = GetParamBitWidth(VariableDefs::RIGHT_INPUTS, adderNb);
            maxShiftLeft = maxShiftLeft == inputBW_ ? leftComputedBW : maxShiftLeft;
            maxShiftRight = maxShiftRight == inputBW_ ? rightComputedBW : maxShiftRight;
            
            outputFile_ << std::string(nbTabs, '\t') << "\"left_theoretical_bitwidth\": " << maxShiftLeft << ",\n";
            outputFile_ << std::string(nbTabs, '\t') << "\"left_computed_bitwidth\": " << leftComputedBW << ",\n";
            outputFile_ << std::string(nbTabs, '\t') << "\"right_theoretical_bitwidth\": " << maxShiftRight << ",\n";
            outputFile_ << std::string(nbTabs, '\t') << "\"right_computed_bitwidth\": " << rightComputedBW << "\n";
            
            nbTabs--;
            outputFile_ << std::string(nbTabs, '\t') << "},\n";
            adderNb++;
        }
    }
    
    outputFile_ << std::string(nbTabs, '\t') << "\"output\": \"ADDER_" << (GetTotalAdderCount() - 1) << "\",\n";
    
    // Output encoding for each target
    outputFile_ << std::string(nbTabs, '\t') << "\"encoding\": {\n";
    nbTabs++;
    for (size_t i = 0; i < targets_.size(); i++) {
        outputFile_ << std::string(nbTabs, '\t') << "\"" << scmDesigns_[i].first << "\": [";
        for (size_t j = 0; j < muxEncoding[i].size(); j++) {
            if (j > 0) outputFile_ << ", ";
            outputFile_ << muxEncoding[i][j];
        }
        outputFile_ << "]";
        if (i < targets_.size() - 1) {
            outputFile_ << ",";
        }
        outputFile_ << "\n";
    }
    nbTabs--;
    outputFile_ << std::string(nbTabs, '\t') << "}\n";
    
    nbTabs--;
    outputFile_ << std::string(nbTabs, '\t') << "}\n";

    outputFile_ << "}\n";
}

unsigned int JSONDumper::GetParamBitWidth(
    const VariableDefs varType,
    const unsigned int adderIdx) const
{
    // Get the number of variables per adder (should be consistent)
    const unsigned int nbVarsPerAdder = layers_.back().adders.back().variables.size();
    const auto baseIdx = varToIdxMap_.at(varType) + adderIdx * nbVarsPerAdder;
    return solutionNode_.variableBitWidths[baseIdx];
}

unsigned int JSONDumper::GetTotalAdderCount() const
{
    unsigned int totalAdders = 0;
    for (const auto& layer : layers_) {
        totalAdders += layer.adders.size();
    }
    return totalAdders;
}