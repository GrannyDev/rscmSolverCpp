//
// Created by smith on 23/02/25.
//

#ifndef NODE_H
#define NODE_H
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "Layer.h"
#include "Rcm.h"

class Node {
public:
    Node(size_t nbBitsToAlloc, size_t nbTargets, size_t nbAdders, size_t nbParameters, size_t nbMuxes);
    ~Node() = default;

    void InitializeMinShiftSavings(const std::vector<Layer>& layers);

    Rcm rscm;
    std::vector<unsigned int> scmIndexes;
    std::vector<bool> isPlusMinus;
    std::vector<unsigned int> minShiftSavings;
    unsigned int cost;
};



#endif //NODE_H
