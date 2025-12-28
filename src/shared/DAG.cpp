//
// Created by smith on 30/03/25.
//

#include "DAG.h"

DAG::DAG(const unsigned int nbBitPerNode, const unsigned int nbAdders, const unsigned int nbVariables) {
    set = boost::dynamic_bitset<>(nbBitPerNode);
    maxOutputValue.reserve(nbVariables);
    isMinus.reserve(nbAdders);
}