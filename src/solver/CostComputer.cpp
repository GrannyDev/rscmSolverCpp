//
// Created by smith on 28/12/2025.
//

#include "CostComputer.h"

#include <iostream>

#include "Solver.h"

unsigned int CostComputer::MuxCountComputer::merge(RSCM& node, DAG const& scm) const
{
    // 1) Merge resource sets and update bounds
    node.rscm.set |= scm.set;

    // Prepare to count
    size_t bitPos = 0;
    unsigned int totalParams = 0;
    unsigned int muxCount = 0;

    // Iterate layers → adders → variables
    for (const auto& layer : solver->layers) {
        for (const auto& adder : layer.adders) {
            unsigned paramInAdderIdx = 0;
            for (const auto& param : adder.variables) {
                // Count how many bits in this parameter are set
                unsigned int bitsSet = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (node.rscm.set.test(bitPos + j)) {
                        ++bitsSet;
                    }
                }
                bitPos += param.possibleValuesFusion.size();

                // If more than one bit, we need (bitsSet - 1) multiplexers
                if (bitsSet > 1 && solver->idxToVarMap.at(paramInAdderIdx) != VariableDefs::RIGHT_MULTIPLIER)
                {
                    // compute the number of number of multiplexers needed
                    muxCount += bitsSet - 1;
                }

                ++totalParams;
                ++paramInAdderIdx;
            }
        }
    }

    // Store result
    node.cost = muxCount;
    return muxCount;
}

unsigned int CostComputer::MuxBitsComputer::merge(RSCM& node, DAG const& scm) const
{
    // 1) Merge resource sets and update bounds
    node.rscm.set |= scm.set;

    // Prepare to count
    size_t bitPos = 0;
    unsigned int totalParams = 0;
    unsigned int bitCount = 0;

    // Iterate layers → adders → variables
    for (const auto& layer : solver->layers) {
        for (const auto& adder : layer.adders) {
            for (const auto& param : adder.variables) {
                // Count how many bits in this parameter are set
                unsigned int bitsSet = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (node.rscm.set.test(bitPos + j)) {
                        ++bitsSet;
                    }
                }
                bitPos += param.possibleValuesFusion.size();

                // If more than one bit, we need (bitsSet - 1) multiplexers
                if (bitsSet > 1)
                {
                    // compute the number of bits to select instead
                    bitCount += std::ceil(std::log2(bitsSet));
                }

                ++totalParams;
            }
        }
    }

    // Store result
    node.cost = bitCount;
    return bitCount;
}

unsigned int CostComputer::LutsCostComputer::merge(RSCM& node, DAG const& scm) const
{
    unsigned int lutCost = 0;

    // 1) Merge resource sets and update bounds
    node.rscm.set |= scm.set;
    std::transform(
        node.rscm.maxOutputValue.begin(),
        node.rscm.maxOutputValue.end(),
        scm.maxOutputValue.begin(),
        node.rscm.maxOutputValue.begin(),
        [](const int a, const int b) { return std::max(a, b); }
    ); // max output values
    std::transform(
        node.rscm.minOutputValue.begin(),
        node.rscm.minOutputValue.end(),
        scm.minOutputValue.begin(),
        node.rscm.minOutputValue.begin(),
        [](const int a, const int b) { return std::min(a, b); }
    ); // min output values

    // update bitwidths
    for (size_t i = 0; i < node.variableBitWidths.size(); ++i) {
        node.variableBitWidths[i] = std::max(BitLength(node.rscm.maxOutputValue[i]), BitLength(node.rscm.minOutputValue[i]));
    }

    auto updateMinShift = [&](const unsigned idx) {
        // Use coefficient's trailing zeros, not the multiplied values'
        // This is correct because any input in the range could have 0 trailing zeros
        node.minShiftSavings[idx] = std::min(node.minShiftSavings[idx], scm.coefficientTrailingZeros[idx]);
    };

    auto muxTreeLuts = [&](const unsigned inputs, const unsigned bw) -> unsigned int {
        auto baseCost = [](const unsigned k) -> double {
            switch (k) {
            case 0:
            case 1:
                return 0.0;
            case 2:
                return 0.5;
            case 3:
            case 4:
                return 1.0;
            case 5:
                return 1.5;
            case 6:
            case 7:
            case 8:
                return 2.0;
            case 16:
                return 4.0;
            default:
                return 0.0;
            }
        };

        double perBitCost = 0.0;
        unsigned int remaining = inputs;

        while (remaining > 16) {
            perBitCost += baseCost(16);
            remaining = remaining - 16 + 1;
        }

        if (remaining == 16) {
            perBitCost += baseCost(16);
        } else if (remaining > 8) {
            perBitCost += baseCost(8);
            remaining = remaining - 8 + 1;
            perBitCost += baseCost(remaining);
        } else {
            perBitCost += baseCost(remaining);
        }

        const double total = std::ceil(perBitCost * static_cast<double>(bw));
        return static_cast<unsigned int>(total);
    };

    auto muxTreeLutsOverhead = [&](const unsigned inputs, const unsigned bw) {
        if (inputs <= 1 || bw == 0) return 0u;
        unsigned tokens = inputs;
        unsigned luts6 = 0;
        unsigned total = 0;
        while (tokens > solver->lutWidth_) {
            ++luts6;
            tokens = tokens - (solver->lutWidth_ - 1);
        }
        if (tokens >= 3) ++luts6; // one last LUT
        else total += (bw + 2 - 1) / 2; // decompose in luts5
        total += luts6 * bw;
        return total; // final LUT handles the rest
    };

    size_t bitPos = 0;
    unsigned int adderIdx = 0;
    unsigned int paramGlobalIdx = 0;
    unsigned int nbEncodingBits = 0;
    std::array<std::pair<unsigned int, unsigned int>, 6> muxTracker; // indexed by static_cast<size_t>(VariableDefs)
    std::array<std::vector<int>, 6> selectedValues; // selected possibleValuesFusion per VariableDefs

    // Iterate layers → adders → variables
    for (const auto& layer : solver->layers) {
        for (const auto& adder : layer.adders) {
            // zero the tracker cheaply
            for (auto& entry : muxTracker) entry = {0u, 0u};
            for (auto& entry : selectedValues) entry.clear();
            // Track plus-minus flag
            node.isPlusMinus[adderIdx] = node.isPlusMinus[adderIdx] || node.rscm.isMinus[adderIdx] != scm.isMinus[adderIdx];

            unsigned paramInAdderIdx = 0;
            for (const auto& param : adder.variables) {
                // 2) Update shift savings early
                updateMinShift(paramGlobalIdx);

                // 3) Count the number of bits at one for this variable
                unsigned bitCount = 0;
                std::vector<int> selected;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (node.rscm.set.test(bitPos + j)) {
                        ++bitCount;
                        selected.push_back(param.possibleValuesFusion[j]);
                    }
                }
                bitPos += param.possibleValuesFusion.size();

                if (bitCount > 1) nbEncodingBits += std::ceil(std::log2(bitCount));

                // 5) Handle RIGHT_MULTIPLIER adder cost (an adder)
                const auto varType = solver->idxToVarMap.at(paramInAdderIdx);
                if (varType == VariableDefs::RIGHT_MULTIPLIER) {
                    const size_t mapSize = solver->varToIdxMap.size();
                    const auto rightMultIdx =  solver->varToIdxMap.at(VariableDefs::RIGHT_MULTIPLIER) + adderIdx * mapSize;
                    const auto leftIdx  = solver->varToIdxMap.at(VariableDefs::LEFT_SHIFTS) + adderIdx * mapSize;
                    const auto rightIdx = solver->varToIdxMap.at(VariableDefs::RIGHT_SHIFTS) + adderIdx * mapSize;
                    lutCost += node.variableBitWidths[rightMultIdx];
                    if (node.variableBitWidths[rightMultIdx] != node.variableBitWidths[leftIdx] && node.variableBitWidths[rightMultIdx] != node.variableBitWidths[rightIdx]) lutCost--; // no need for a carry lut
                } else {
                    const unsigned normalizedBw = node.variableBitWidths[paramGlobalIdx] - node.minShiftSavings[paramGlobalIdx];
                    muxTracker[static_cast<size_t>(varType)] = {bitCount > 1 ? bitCount : 0, normalizedBw};
                    selectedValues[static_cast<size_t>(varType)] = std::move(selected);
                }

                ++paramGlobalIdx;
                ++paramInAdderIdx;
            }

            // output never merged
            const auto outputPathInputs = muxTracker[static_cast<size_t>(VariableDefs::OUTPUTS_SHIFTS)].first;
            const auto outputPathBw = muxTracker[static_cast<size_t>(VariableDefs::OUTPUTS_SHIFTS)].second;
            lutCost += muxTreeLuts(outputPathInputs, outputPathBw);

            // prepare for merge logic
            size_t leftMergeCandiate = static_cast<size_t>(VariableDefs::LEFT_SHIFTS);
            size_t rightMergeCandiate = static_cast<size_t>(VariableDefs::RIGHT_SHIFTS);
            if (muxTracker[leftMergeCandiate].first == 0) {
                leftMergeCandiate = static_cast<size_t>(VariableDefs::LEFT_INPUTS);
            } else {
                lutCost += muxTreeLuts(
                    muxTracker[static_cast<size_t>(VariableDefs::LEFT_INPUTS)].first,
                    muxTracker[static_cast<size_t>(VariableDefs::LEFT_INPUTS)].second
                );
            }
            if (muxTracker[rightMergeCandiate].first == 0) {
                rightMergeCandiate = static_cast<size_t>(VariableDefs::RIGHT_INPUTS);
            } else {
                lutCost += muxTreeLuts(
                    muxTracker[static_cast<size_t>(VariableDefs::RIGHT_INPUTS)].first,
                    muxTracker[static_cast<size_t>(VariableDefs::RIGHT_INPUTS)].second
                );
            }
            // look for mux2
            bool bothMux2 = true;
            if (muxTracker[leftMergeCandiate].first > 2) {
                lutCost += muxTreeLuts(
                    muxTracker[leftMergeCandiate].first,
                    muxTracker[leftMergeCandiate].second
                );
                bothMux2 = false;
            }
            if (muxTracker[rightMergeCandiate].first > 2) {
                lutCost += muxTreeLuts(
                    muxTracker[rightMergeCandiate].first,
                    muxTracker[rightMergeCandiate].second
                );
                bothMux2 = false;
            }
            // if both are mux2, determine if both can be merged or which one
            if (bothMux2) {
                const auto leftInputs = muxTracker[leftMergeCandiate].first;
                const auto rightInputs = muxTracker[rightMergeCandiate].first;
                if (leftInputs <= 1 || rightInputs <= 1) {
                    // nothing to add
                } else {
                    auto extractShifts = [&](const size_t varIdx) {
                        std::vector<int> shifts;
                        const auto varType = static_cast<VariableDefs>(varIdx);
                        for (const int v : selectedValues[varIdx]) {
                            if (varType == VariableDefs::LEFT_SHIFTS ||
                                varType == VariableDefs::RIGHT_SHIFTS ||
                                varType == VariableDefs::OUTPUTS_SHIFTS) {
                                shifts.push_back(v);
                            } else if (varType == VariableDefs::LEFT_INPUTS ||
                                       varType == VariableDefs::RIGHT_INPUTS) {
                                if (v < 0) {
                                    const int shift = (v == -1) ? 0 : std::abs(v + 1);
                                    shifts.push_back(shift);
                                } else {
                                    // Encode adder inputs as negative values to avoid clashing with shifts.
                                    shifts.push_back(-static_cast<int>(v) - 1);
                                }
                            }
                        }
                        return shifts;
                    };

                    bool hasCommonShift = false;
                    if (!node.isPlusMinus[adderIdx]) {
                        const auto leftShifts = extractShifts(leftMergeCandiate);
                        const auto rightShifts = extractShifts(rightMergeCandiate);
                        for (const int l : leftShifts) {
                            for (const int r : rightShifts) {
                                if (l == r) {
                                    hasCommonShift = true;
                                    break;
                                }
                            }
                            if (hasCommonShift) break;
                        }
                    }

                    if (!hasCommonShift) {
                        const auto leftBw = muxTracker[leftMergeCandiate].second;
                        const auto rightBw = muxTracker[rightMergeCandiate].second;
                        if (leftBw <= rightBw) {
                            lutCost += muxTreeLuts(leftInputs, leftBw);
                        } else {
                            lutCost += muxTreeLuts(rightInputs, rightBw);
                        }
                    }
                }
            }

            ++adderIdx;
        }
    }

    // add decoding overhead if nbEncodingBits > nbMinEncodingBits_
    if (nbEncodingBits > solver->nbMinEncodingBits_)
    {
        lutCost +=  muxTreeLutsOverhead(solver->nbMinEncodingBits_, nbEncodingBits);
    }

    // Store result
    node.cost = lutCost;
    return lutCost;
}

unsigned int CostComputer::FPGADelayComputer::merge(RSCM& /*node*/, DAG const& /*scm*/) const
{
    throw std::logic_error("ASICDelayComputer::merge not implemented");
}

unsigned int CostComputer::ASICDelayComputer::merge(RSCM& /*node*/, DAG const& /*scm*/) const
{
    throw std::logic_error("ASICDelayComputer::merge not implemented");
}

unsigned int CostComputer::AreaCostComputer::merge(RSCM& node, DAG const& scm) const
{
    unsigned int fineGrainCost = 0;

    // 1) Merge resource sets and update bounds
    node.rscm.set |= scm.set;
    std::transform(
        node.rscm.maxOutputValue.begin(),
        node.rscm.maxOutputValue.end(),
        scm.maxOutputValue.begin(),
        node.rscm.maxOutputValue.begin(),
        [](const int a, const int b) { return std::max(a, b); }
    ); // max output values
    std::transform(
        node.rscm.minOutputValue.begin(),
        node.rscm.minOutputValue.end(),
        scm.minOutputValue.begin(),
        node.rscm.minOutputValue.begin(),
        [](const int a, const int b) { return std::min(a, b); }
    ); // min output values

    // update bitwidths
    for (size_t i = 0; i < node.variableBitWidths.size(); ++i) {
        node.variableBitWidths[i] = std::max(BitLength(node.rscm.maxOutputValue[i]), BitLength(node.rscm.minOutputValue[i]));
    }

    const auto leftShiftBase  = solver->varToIdxMap.at(VariableDefs::LEFT_SHIFTS);
    const auto rightShiftBase = solver->varToIdxMap.at(VariableDefs::RIGHT_SHIFTS);
    const size_t mapSize = solver->varToIdxMap.size();

    auto updateMinShift = [&](const unsigned idx) {
        // Use coefficient's trailing zeros, not the multiplied values'
        // This is correct because any input in the range could have 0 trailing zeros
        const unsigned coeffTZ = scm.coefficientTrailingZeros[idx];
        node.minShiftSavings[idx] = std::min(node.minShiftSavings[idx], coeffTZ);
    };

    auto computeMuxCost = [&](const unsigned bitsCount, const unsigned idx) {
        const unsigned int nbBits = node.variableBitWidths[idx] - node.minShiftSavings[idx];
        return 11u * (bitsCount - 1) * nbBits;
    };

    size_t bitPos = 0;
    unsigned int adderIdx = 0;
    unsigned int paramGlobalIdx = 0;
    unsigned int nbEncodingBits = 0;

    // Iterate layers → adders → variables
    for (const auto& layer : solver->layers) {
        for (const auto& adder : layer.adders) {
            // Track plus-minus flag
            node.isPlusMinus[adderIdx] =
                node.isPlusMinus[adderIdx] || node.rscm.isMinus[adderIdx] != scm.isMinus[adderIdx];

            unsigned paramInAdderIdx = 0;
            for (const auto& param : adder.variables) {
                // 2) Update shift savings early
                updateMinShift(paramGlobalIdx);

                // 3) Count the number of bits at one for this variable
                unsigned bitCount = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (node.rscm.set.test(bitPos + j)) ++bitCount;
                }
                bitPos += param.possibleValuesFusion.size();

                if (bitCount > 1)
                {
                    nbEncodingBits += std::ceil(std::log2(bitCount));
                    // 4) If more than one bit and not an adder -> we have a multiplexer
                    if (solver->idxToVarMap.at(paramInAdderIdx) != VariableDefs::RIGHT_MULTIPLIER) {
                        fineGrainCost += computeMuxCost(bitCount, paramGlobalIdx);
                    }
                }

                // 5) Handle RIGHT_MULTIPLIER adder cost (an adder)
                if (solver->idxToVarMap.at(paramInAdderIdx) == VariableDefs::RIGHT_MULTIPLIER) {
                    // determine coefficient depending on the adder type
                    unsigned coefFA = 20u;
                    unsigned coefHA = 13u;
                    if (node.isPlusMinus[adderIdx]) {
                        coefFA = 29u;
                        coefHA = 21u;
                    } else if (node.rscm.isMinus[adderIdx]) {
                        coefFA = 24u;
                        coefHA = 16u;
                    }

                    // compute the minimum shift savings (number of trailing zeros)
                    const auto leftIdx  = leftShiftBase  + adderIdx * mapSize;
                    const auto rightIdx = rightShiftBase + adderIdx * mapSize;
                    const unsigned minShift = std::min(
                        node.minShiftSavings[leftIdx],
                        node.minShiftSavings[rightIdx]
                    );

                    // compute FA/HA costs
                    if (node.rscm.minOutputValue[rightIdx] != 0 && node.rscm.maxOutputValue[rightIdx] != 0) {
                        unsigned wa = node.variableBitWidths[leftIdx] - minShift;
                        unsigned wb = node.variableBitWidths[rightIdx] - minShift;
                        unsigned int shiftA = node.minShiftSavings[leftIdx] - minShift;
                        unsigned int shiftB = node.minShiftSavings[rightIdx] - minShift;
                        auto [fst, snd] = std::minmax(shiftA, shiftB);
                        const unsigned int diff = snd - fst;

                        unsigned fa = 0, ha = 0;
                        if (!node.isPlusMinus[adderIdx] && !node.rscm.isMinus[adderIdx]) {
                            if (node.rscm.minOutputValue[leftIdx] != 0 && node.rscm.maxOutputValue[leftIdx] != 0) {
                                fa = std::max(wa, wb) - diff - 1;
                                ha = 1;
                            }
                        } else {
                            if (node.rscm.minOutputValue[leftIdx] == 0 && node.rscm.maxOutputValue[leftIdx] == 0) { // left (a) is zero
                                ha = wb - node.minShiftSavings[rightIdx];
                            } else if (shiftB >= shiftA) {
                                fa = std::max(wa, wb) - diff;
                            } else if (node.minShiftSavings[leftIdx] >= wb) {
                                fa = std::max(wa, wb) - diff - 1;
                                ha = 1 + wb;
                            } else {
                                fa = std::max(wa, wb) - diff;
                                ha = diff;
                            }
                        }
                        fineGrainCost += coefFA * fa + coefHA * ha;
                    }
                }

                ++paramGlobalIdx;
                ++paramInAdderIdx;
            }
            ++adderIdx;
        }
    }

    // add decoding overhead if nbEncodingBits > nbMinEncodingBits_
    if (nbEncodingBits > solver->nbMinEncodingBits_)
    {
        fineGrainCost += 442u + 1.35 * (nbEncodingBits * (1 << solver->nbMinEncodingBits_));
    }

    // Store result
    node.cost = fineGrainCost;
    return fineGrainCost;
}

unsigned int CostComputer::BitLength(const int maxConstant)
{
    return static_cast<unsigned int>(std::ceil(std::log2(std::max(std::abs(maxConstant + 1), std::abs(maxConstant))))) + 1;
}