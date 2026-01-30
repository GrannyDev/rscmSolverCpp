//
// Created by smith on 23/02/25.
//

#ifndef NODE_H
#define NODE_H
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "../shared/Layer.h"
#include "../shared/DAG.h"

/**
 * @class RSCM
 * @brief Represents an RSCM.
 */
class RSCM {
public:
    /**
     * @brief Constructs a RSCM object with the specified parameters.
     *
     * @param nbBitsToAlloc The number of bits to represent the underlying DAG as a bitset.
     * @param nbTargets The number of target constants.
     * @param nbAdders The number of adders in the underlying DAG.
     * @param nbVariables The number of variables of the underlying DAG.
     * @param nbMuxes The number of multiplexers in the RSCM.
     */
    RSCM(size_t nbBitsToAlloc, size_t nbTargets, size_t nbAdders, size_t nbVariables, size_t nbMuxes);

    /**
     * @brief Default destructor for the Node class.
     */
    ~RSCM() = default;

    /**
     * @brief Initializes the minimum shift savings (how many LSB per variable output are at 0)
     * based on the provided layers.
     *
     * @param layers A vector of Layer describing the adders.
     */
    void InitializeMinShiftSavings(const std::vector<Layer>& layers);

    /**
     * @brief The underlying DAG.
     */
    DAG rscm;

    /**
     * @brief A vector storing the SCM DAGs indexes merged into this RSCM (sorted by first to last).
     */
    std::vector<unsigned int> scmIndexes;

    /**
     * @brief A vector of boolean indicating whether the corresponding adder is an adder/subtractor.
     */
    std::vector<bool> isPlusMinus;

    /**
     * @brief A vector storing the minimum number of trailing zeros for each variable output.
     */
    std::vector<unsigned int> minShiftSavings;

    /**
     * @brief A vector representig the bitwidth of each variable output.
     */
    std::vector<unsigned int> variableBitWidths;

    /**
     * @brief The cost associated with the RSCM.
     */
    unsigned int cost;
};

#endif //NODE_H