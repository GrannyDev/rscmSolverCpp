//
// Created by smith on 30/03/25.
//

#ifndef RCM_H
#define RCM_H

#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <vector>

/**
 * @class DAG
 * @brief Represents a class for managing DAGs and parameters necessary for computing the fine-grained cost function.
 */
class DAG {
public:
    /**
     * @brief Constructs a Rcm object with the specified parameters.
     *
     * @param nbBitPerNode The number of bits to represent the DAG as a bitset.
     * @param nbAdders The number of adders in the DAG.
     * @param nbVariables The number of variables of the DAG.
     */
    explicit DAG(unsigned int nbBitPerNode, unsigned int nbAdders, unsigned int nbVariables);

    /**
     * @brief A bitset representing the DAG.
     */
    boost::dynamic_bitset<> set;

    /**
     * @brief A vector storing parameters maximum output values.
     */
    std::vector<int> maxOutputValue;

    /**
     * @brief A vector of boolean flags indicating whether the corresponding adder is negative.
     */
    std::vector<bool> isMinus;
};

#endif // RCM_H