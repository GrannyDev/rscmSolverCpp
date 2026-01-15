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

    auto muxTreeLuts = [&](const unsigned inputs, const unsigned bw) {
        if (inputs <= 1 || bw == 0) return 0u;
        unsigned tokens = inputs + static_cast<unsigned>(std::ceil(std::log2(static_cast<double>(inputs))));
        unsigned luts6 = 0;
        unsigned total = 0;
        while (tokens > solver->lutWidth_) {
            ++luts6;
            tokens = tokens - (solver->lutWidth_ - 1);
        }
        if (tokens == 6) ++luts6; // one last LUT
        else total += (bw + 2 - 1) / 2; // decompose in luts5
        total += luts6 * bw;
        return total; // final LUT handles the rest
    };

    size_t bitPos = 0;
    unsigned int adderIdx = 0;
    unsigned int paramGlobalIdx = 0;
    unsigned int nbEncodingBits = 0;
    std::array<std::pair<unsigned int, unsigned int>, 6> muxTracker; // indexed by static_cast<size_t>(VariableDefs)

    // Iterate layers → adders → variables
    for (const auto& layer : solver->layers) {
        for (const auto& adder : layer.adders) {
            // zero the tracker cheaply
            for (auto& entry : muxTracker) entry = {0u, 0u};
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

                if (bitCount > 1) nbEncodingBits += std::ceil(std::log2(bitCount));

                // 5) Handle RIGHT_MULTIPLIER adder cost (an adder)
                if (solver->idxToVarMap.at(paramInAdderIdx) == VariableDefs::RIGHT_MULTIPLIER) {
                    const size_t mapSize = solver->varToIdxMap.size();
                    const auto rightMultIdx =  solver->varToIdxMap.at(VariableDefs::RIGHT_MULTIPLIER) + adderIdx * mapSize;
                    lutCost += node.variableBitWidths[rightMultIdx];
                    if (node.minShiftSavings[rightMultIdx] == 0) lutCost--; // no need for a carry lut
                } else {
                    const auto varType = solver->idxToVarMap.at(paramInAdderIdx);
                    muxTracker[static_cast<size_t>(varType)] = {bitCount > 1 ? bitCount : 0, node.variableBitWidths[paramGlobalIdx]};
                }

                ++paramGlobalIdx;
                ++paramInAdderIdx;
            }

            // handle mux merge logic
            // 1) find the biggest one that can be merged into the add or sub
            if (!node.isPlusMinus[adderIdx]) // can't merge if it's both
            {
                unsigned int maxLUTs = 0;
                auto maxIdx = static_cast<size_t>(VariableDefs::OUTPUTS_SHIFTS); // default to output (won't be used)
                for (size_t i = 0; i < muxTracker.size(); ++i) {
                    if (i == static_cast<size_t>(VariableDefs::OUTPUTS_SHIFTS)) continue; // never merge output mux into adder
                    const auto [nbInputs, bw] = muxTracker[i];
                    if (nbInputs <= 1 || bw == 0) continue;
                    const auto selBits = static_cast<unsigned>(std::ceil(std::log2(static_cast<double>(nbInputs))));
                    if (nbInputs + selBits <= solver->nbInputsLeftByAdderLut_ && bw * (nbInputs + selBits) > maxLUTs) {
                        maxLUTs = muxTreeLuts(nbInputs + selBits, bw);
                        maxIdx = i;
                    }
                }
                // 2) remove the merged mux if there's one
                if (maxLUTs > 0)
                {
                    muxTracker[maxIdx].first = 0;
                }
            }

            // 3) compute the number of LUTs per right/left path
            auto fuseCost = [&](const unsigned aCnt, const unsigned aBw, const unsigned bCnt, const unsigned bBw) {
                const unsigned separate = muxTreeLuts(aCnt, aBw) + muxTreeLuts(bCnt, bBw);
                const unsigned combined = muxTreeLuts(aCnt + bCnt, bBw);
                return std::min(separate, combined);
            };

            const unsigned leftCost = fuseCost(
                muxTracker[static_cast<size_t>(VariableDefs::LEFT_INPUTS)].first,  muxTracker[static_cast<size_t>(VariableDefs::LEFT_INPUTS)].second,
                muxTracker[static_cast<size_t>(VariableDefs::LEFT_SHIFTS)].first,  muxTracker[static_cast<size_t>(VariableDefs::LEFT_SHIFTS)].second
            );
            const unsigned rightCost = fuseCost(
                muxTracker[static_cast<size_t>(VariableDefs::RIGHT_INPUTS)].first, muxTracker[static_cast<size_t>(VariableDefs::RIGHT_INPUTS)].second,
                muxTracker[static_cast<size_t>(VariableDefs::RIGHT_SHIFTS)].first, muxTracker[static_cast<size_t>(VariableDefs::RIGHT_SHIFTS)].second
            );

            const auto outputPathInputs = muxTracker[static_cast<size_t>(VariableDefs::OUTPUTS_SHIFTS)].first;
            const auto outputPathBw = muxTracker[static_cast<size_t>(VariableDefs::OUTPUTS_SHIFTS)].second;

            lutCost += leftCost;
            lutCost += rightCost;
            lutCost += muxTreeLuts(outputPathInputs, outputPathBw);

            ++adderIdx;
        }
    }

    // add decoding overhead if nbEncodingBits > nbMinEncodingBits_
    if (nbEncodingBits > solver->nbMinEncodingBits_)
    {
        lutCost +=  muxTreeLuts(solver->nbMinEncodingBits_, nbEncodingBits);
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