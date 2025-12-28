//
// Created by smith on 28/12/2025.
//

#ifndef RSCMSOLVER_ICOSTCOMPUTER_H
#define RSCMSOLVER_ICOSTCOMPUTER_H
#include "RSCM.h"

class Solver;

/**
 * @struct ICostComputer
 * @brief Abstract base class for cost computation strategies.
 */
struct ICostComputer {
    /**
     * @brief Constructor for ICostComputer.
     * @param solver Pointer to the Solver instance.
     */
    explicit ICostComputer(Solver* solver) : solver(solver) {}

    /**
     * @brief Virtual destructor.
     */
    virtual ~ICostComputer() = default;

    /**
     * @brief Computes the cost of adding a DAG to a thread node.
     * @param threadNode The RSCM to compute the cost for (once merged).
     * @param scmToMerge The DAG to merge.
     * @return The computed cost.
     */
    virtual unsigned int merge(RSCM& threadNode, DAG const& scmToMerge) const = 0;

    Solver* solver; ///< Pointer to the Solver instance.
};

#endif //RSCMSOLVER_ICOSTCOMPUTER_H
