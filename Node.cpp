//
// Created by smith on 23/02/25.
//

#include "Node.h"

Node::Node(const size_t nbBitsToAlloc, const size_t nbTargets, const size_t nbAdders, const size_t nbParameters, const size_t nbMuxes) : rscm(nbBitsToAlloc, nbAdders, nbMuxes)
{
    scmIndexes.resize(nbTargets);
    isPlusMinus.resize(nbAdders, false);
    minShiftSavings.resize(nbParameters, std::numeric_limits<unsigned int>::max());
    cost = 0;
}

void Node::InitializeMinShiftSavings(const std::vector<Layer>& layers)
{
    unsigned int totalParameters = 0;
    for (const auto& l: layers)
    {
        for (const auto & adder : l.adders)
        {
            for (const auto& p: adder.parameters)
            {
                minShiftSavings[totalParameters] = static_cast<unsigned>(rscm.maxMuxOutputValue[totalParameters] ? __builtin_ctz(rscm.maxMuxOutputValue[totalParameters]) : 0);
                totalParameters++;
            }
        }
    }
}
