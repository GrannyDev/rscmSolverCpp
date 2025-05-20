//
// Created by smith on 30/03/25.
//

#ifndef RCM_H
#define RCM_H
#include <boost/dynamic_bitset/dynamic_bitset.hpp>


class Rcm {
public:
    explicit Rcm(unsigned int nbBitPerNode, unsigned int nbAdders, unsigned int nbParameters);

    boost::dynamic_bitset<> set;
    std::vector<int> maxMuxOutputValue;
    std::vector<bool> isMinus;
};



#endif //RCM_H
