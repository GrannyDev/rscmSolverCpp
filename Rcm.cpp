//
// Created by smith on 30/03/25.
//

#include "Rcm.h"

Rcm::Rcm(const unsigned int nbBitPerNode, const unsigned int nbAdders, const unsigned int nbParameters) {
    set = boost::dynamic_bitset<>(nbBitPerNode);
    maxMuxOutputValue.reserve(nbParameters);
    isMinus.reserve(nbAdders);
}
