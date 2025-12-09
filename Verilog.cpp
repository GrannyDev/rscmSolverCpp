//
// Created by smith on 09/12/2025.
//

#include "Verilog.h"

#include <filesystem>
#include <iostream>

VerilogGenerator::VerilogGenerator
    (
        const RSCM& solutionNode,
        const std::string& outputUri,
        const std::vector<Layer>& layers,
        const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap,
        const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap,
        const bool overwrite
    )
    : solutionNode_(solutionNode), layers_(layers), idxToVarMap_(idxToVarMap), varToIdxMap_(varToIdxMap)
{
    // check if output already exists
    if (std::filesystem::exists(outputUri) && !overwrite)
    {
        throw std::runtime_error("Output file already exists.");
    }
    // overwrite
    outputFile_ = std::ofstream(outputUri, std::ofstream::trunc);

    // test for now
    PrintAdderModules();
}

void VerilogGenerator::PrintAdderModules()
{
    // Prepare to count
    size_t bitPos = 0;
    unsigned int adderIdx = 0;
    unsigned int paramGlobalIdx = 0;

    const auto leftShiftBase  = varToIdxMap_.at(VariableDefs::LEFT_SHIFTS);
    const auto rightShiftBase = varToIdxMap_.at(VariableDefs::RIGHT_SHIFTS);
    const auto outputBase = varToIdxMap_.at(VariableDefs::OUTPUTS_SHIFTS);
    const size_t nbVarsPerAdder = layers_[0].adders[0].variables.size();

    // Lambda to compute bit width from max value
    auto computeBitWidth = [](const int maxVal, const bool is_signed) -> unsigned int {
        if (maxVal == 0) return 1;
        if (!is_signed) return static_cast<unsigned int>(std::ceil(std::log2(std::abs(maxVal))));
        return static_cast<unsigned int>(std::ceil(std::log2(std::max(std::abs(maxVal + 1), std::abs(maxVal))))) + 1;
    };

    // Iterate layers → adders → variables
    // This will determinate our module inputs
    for (const auto& layer : layers_)
    {
        for (const auto& adder : layer.adders)
        {
            std::unordered_map<unsigned int, unsigned int> muxInputsMap;
            unsigned int leftInputBitWidth = 0;
            unsigned int rightInputBitWidth = 0;
            unsigned int outputBitWidth = 0;
            unsigned paramInAdderIdx = 0;

            for (const auto& param : adder.variables) {
                // Count how many bits in this parameter are set
                unsigned int bitsSet = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (solutionNode_.rscm.set.test(bitPos + j)) {
                        ++bitsSet;
                    }
                }
                bitPos += param.possibleValuesFusion.size();

                // If more than one bit, we need a multiplexer
                if (bitsSet > 1)
                {
                    // save the number of inputs in the muxInputsMap;
                    muxInputsMap[paramGlobalIdx] = bitsSet;
                }
                // the adder input
                if (idxToVarMap_.at(paramInAdderIdx) == VariableDefs::RIGHT_MULTIPLIER)
                {
                    const auto leftIdx  = leftShiftBase  + adderIdx * nbVarsPerAdder;
                    const auto rightIdx = rightShiftBase + adderIdx * nbVarsPerAdder;
                    const auto outputIdx = outputBase + adderIdx * nbVarsPerAdder;
                    const int left_max = solutionNode_.rscm.maxOutputValue[leftIdx];
                    const int right_max = solutionNode_.rscm.maxOutputValue[rightIdx];
                    const int output_max = solutionNode_.rscm.maxOutputValue[outputIdx];
                    std::cout << "output max: " << output_max << std::endl;
                    // compute bitwidth to handle max value
                    leftInputBitWidth = computeBitWidth(left_max, true);
                    rightInputBitWidth = computeBitWidth(right_max, true);
                    outputBitWidth = computeBitWidth(output_max, true);
                }
                ++paramGlobalIdx;
                ++paramInAdderIdx;
            }

            // we now have all the infos to print our add module skeleton
            // let's build the chain
            outputFile_ << "module ADDER_" << adderIdx << " (\n";
            outputFile_ << "\tinput [" << (leftInputBitWidth - 1) << ":0] LEFT_INPUT,\n";
            outputFile_ << "\tinput [" << (rightInputBitWidth - 1) << ":0] RIGHT_INPUT,\n";
            // mux inputs
            for (const auto& [paramIdx, muxInputCount] : muxInputsMap)
            {
                const unsigned int muxSelBits = computeBitWidth(static_cast<int>(muxInputCount) - 1, false);
                outputFile_ << "\tinput [" << (muxSelBits - 1) << ":0] MUX_INPUT_" << paramIdx << ",\n";
            }
            outputFile_ << "\toutput [" << (outputBitWidth - 1) << ":0] OUTPUT\n";
            outputFile_ << ");\n\n";
            outputFile_ << "\t// Adder logic here\n\n";
            outputFile_ << "endmodule\n\n";
        }
        ++adderIdx;
    }
}
