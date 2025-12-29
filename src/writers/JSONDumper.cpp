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
    const std::unordered_map<std::string, std::optional<unsigned int>>& costs,
    const bool isSymmetric,
    const bool overwrite
)
    : solutionNode_(solutionNode), layers_(layers), idxToVarMap_(idxToVarMap), varToIdxMap_(varToIdxMap), inputBW_(inputBW), targets_(targets), scmDesigns_(scmDesigns), costs_(costs), isSymmetric_(isSymmetric), varDefsTostring_(VarDefsToString())
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
    unsigned int totalNbMuxes = 0;
    std::vector<unsigned int> muxNbToMuxBw;
    
    for (const auto& layer : layers_) {
        for (const auto& adder : layer.adders) {
            outputFile_ << std::string(nbTabs, '\t') << "\"ADDER_" << adderNb << "\": {\n";
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
                const unsigned int paramGlobalIdx = paramGlobalIdxOffset + paramIdx;
                const VariableDefs currentParamType = idxToVarMap_.at(paramIdx);
                const unsigned int computedBW = solutionNode_.variableBitWidths[paramGlobalIdx];

                // Collect selected shifts and find max shift
                std::vector<int> selectedShifts;
                std::vector<std::string> inputs;
                unsigned int maxShift = 0;
                for (const int v : param.possibleValuesFusion) {
                    if (solutionNode_.rscm.set.test(bitPos + v + param.zeroPoint)) {
                        if (currentParamType == VariableDefs::OUTPUTS_SHIFTS)
                        {
                            selectedShifts.push_back(v);
                            // Output shift uses power of 2, but only count if v > 0
                            if (v > 0) {
                                maxShift = std::max(maxShift, static_cast<unsigned int>(1 << v));
                            }
                        } else if (currentParamType == VariableDefs::LEFT_INPUTS || currentParamType == VariableDefs::RIGHT_INPUTS) {
                            if (v < 0) {
                                if (v == -1) {
                                    inputs.emplace_back("X");
                                    selectedShifts.push_back(0);
                                } else {
                                    unsigned int shiftAmount = std::abs(v + 1);
                                    inputs.emplace_back("X");
                                    selectedShifts.push_back(static_cast<int>(shiftAmount));
                                    maxShift = std::max(maxShift, shiftAmount);
                                }
                            } else {
                                inputs.push_back("ADDER_" + std::to_string(v));
                                selectedShifts.push_back(0);
                            }
                        } else {
                            selectedShifts.push_back(v);
                            // Input shift uses direct value
                            maxShift = std::max(maxShift, static_cast<unsigned int>(v));
                        }
                    }
                }

                // Check if this is a mux (more than one option selected)
                bool isMux = selectedShifts.size() > 1;
                if (isMux)
                {
                    muxNbToMuxBw.push_back(std::ceil(std::log2(selectedShifts.size())));
                    totalNbMuxes++;
                }
                if (currentParamType == VariableDefs::RIGHT_MULTIPLIER) {
                    // For ADD_SUB, it's a mux if it can do both add and subtract
                    int addCount = 0, subCount = 0;
                    for (const int v : selectedShifts) {
                        if (v == -1) subCount++;
                        else if (v == 1) addCount++;
                    }
                    isMux = addCount > 0 && subCount > 0;
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
                } else
                {
                    theoreticalBW = computedBW; // For inputs, theoretical equals computed
                }

                // Output parameter information
                outputFile_ << std::string(nbTabs, '\t') << "\"" << varDefsTostring_(currentParamType) << "\": {\n";
                nbTabs++;

                if (currentParamType == VariableDefs::LEFT_INPUTS || currentParamType == VariableDefs::RIGHT_INPUTS)
                {
                    // add the inputs values
                    outputFile_ << std::string(nbTabs, '\t') << "\"inputs\": [";
                    bool first = true;
                    for (const std::string& input : inputs) {
                        if (!first) outputFile_ << ", ";
                        outputFile_ << "\"" << input << "\"";
                        first = false;
                    }
                    outputFile_ << "],\n";
                }
                if (currentParamType == VariableDefs::RIGHT_MULTIPLIER) {
                    // Handle ADD_SUB specially
                    outputFile_ << std::string(nbTabs, '\t') << R"("type": ")";
                    
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
                }
                outputFile_ << std::string(nbTabs, '\t') << "\"theoretical_bitwidth\": " << theoreticalBW << ",\n";
                outputFile_ << std::string(nbTabs, '\t') << "\"computed_bitwidth\": " << computedBW;
                if (isMux) {
                    outputFile_ << ",\n" << std::string(nbTabs, '\t') << "\"mux_nb\": " << totalNbMuxes - 1;
                }
                outputFile_ << "\n";

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

    // metadata
    outputFile_ << std::string(nbTabs, '\t') << "\"meta\": {\n";
    nbTabs++;
    outputFile_ << std::string(nbTabs, '\t') << "\"wIn\": " << inputBW_ << ",\n";
    unsigned int wConf = std::floor(std::log2(targets_.size() + 1));
    unsigned int wMuxTotal = 0;
    for (const auto& bw : muxNbToMuxBw) {
        wMuxTotal += bw;
    }
    outputFile_ << std::string(nbTabs, '\t') << "\"wConf\": " << wConf << ",\n";
    outputFile_ << std::string(nbTabs, '\t') << "\"wSelec\": " << wMuxTotal << ",\n";
    outputFile_ << std::string(nbTabs, '\t') << "\"needs_table\": " << (wConf != wMuxTotal ? "true" : "false") << ",\n";
    outputFile_ << std::string(nbTabs, '\t') << "\"wOut\": " << GetParamBitWidth(VariableDefs::OUTPUTS_SHIFTS, GetTotalAdderCount() - 1) << ",\n";
    outputFile_ << std::string(nbTabs, '\t') << R"("output": "ADDER_)" << GetTotalAdderCount() - 1 << "\",\n";
    outputFile_ << std::string(nbTabs, '\t') << "\"is_symmetric\": " << (isSymmetric_ ? "true" : "false") << ",\n";
    outputFile_ << std::string(nbTabs, '\t') << "\"muxes_bw\": [";
    for (size_t i = 0; i < muxNbToMuxBw.size(); i++) {
        if (i > 0) outputFile_ << ", ";
        outputFile_ << muxNbToMuxBw[i];
    }
    outputFile_ << "],\n";

    outputFile_ << std::string(nbTabs, '\t') << "\"costs\": {\n";
    nbTabs++;
    size_t costIdx = 0;
    for (const auto& [name, value] : costs_) {
        outputFile_ << std::string(nbTabs, '\t') << "\"" << name << "\": ";
        if (value.has_value()) {
            outputFile_ << *value;
        } else {
            outputFile_ << "\"not implemented\"";
        }
        if (costIdx < costs_.size() - 1) {
            outputFile_ << ",";
        }
        outputFile_ << "\n";
        ++costIdx;
    }
    nbTabs--;
    outputFile_ << std::string(nbTabs, '\t') << "},\n";
    outputFile_ << std::string(nbTabs, '\t') << "\"constants\": [";
    for (int i = 0; i < targets_.size(); i++) {
        if (i > 0) outputFile_ << ", ";
        outputFile_ << scmDesigns_[i].first;
    }
    outputFile_ << "],\n";

    // Output encoding for each target
    outputFile_ << std::string(nbTabs, '\t') << R"("encoding": [)" << "\n";
    nbTabs++;
    for (size_t targetIndex = 0; targetIndex < targets_.size(); targetIndex++) {
        outputFile_ << std::string(nbTabs, '\t') << "[";
        for (size_t i = 0; i < muxEncoding[targetIndex].size(); i++) {
            if (i > 0) outputFile_ << ", ";
            outputFile_ << muxEncoding[targetIndex][i];
        }
        outputFile_ << "]";
        if (targetIndex < targets_.size() - 1) {
            outputFile_ << ",";
        }
        outputFile_ << "\n";
    }
    nbTabs--;
    outputFile_ << std::string(nbTabs, '\t') << "]\n";
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
