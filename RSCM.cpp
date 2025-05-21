//
// Created by smith on 23/02/25.
//

#include "RSCM.h"

RSCM::RSCM(const size_t nbBitsToAlloc, const size_t nbTargets, const size_t nbAdders, const size_t nbVariables, const size_t nbMuxes)
    : rscm(nbBitsToAlloc, nbAdders, nbMuxes)
{
    // Pre-allocate memory as their size is known in advance.
    scmIndexes.resize(nbTargets);
    isPlusMinus.resize(nbAdders, false);
    minShiftSavings.resize(nbVariables, std::numeric_limits<unsigned int>::max());
    cost = 0;
}

void RSCM::InitializeMinShiftSavings(const std::vector<Layer>& layers)
{
    unsigned int totalVariables = 0;
    for (const auto& l: layers) // Iterate through each layer.
    {
        for (const auto & adder : l.adders) // Iterate through each adder in the layer.
        {
            for ([[maybe_unused]] const auto& _: adder.variables) // Iterate through each variable of the adder.
            {
                // Compute the minimum number of trailing zeros for the variable output.
                // This gives the number of bits that can be saved by shifting the output.
                minShiftSavings[totalVariables] = static_cast<unsigned>(
                    rscm.maxOutputValue[totalVariables] ? __builtin_ctz(rscm.maxOutputValue[totalVariables]) : 0
                );
                totalVariables++;
            }
        }
    }
}