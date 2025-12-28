//
// Created by smith on 28/12/2025.
//

#ifndef RSCMSOLVER_COSTCOMPUTER_H
#define RSCMSOLVER_COSTCOMPUTER_H
#include "ICostComputer.h"


class CostComputer
{
public:
    /**
     * @struct AreaCostComputer
     * @brief Computes costs with fine-grain granularity.
     */
    struct AreaCostComputer final : ICostComputer {
        using ICostComputer::ICostComputer;

        /**
         * @brief Computes the cost of adding a DAG to a node with fine-grain granularity.
         * @param node The RSCM to compute the cost for (once merged).
         * @param scm The DAG to merge.
         * @return The computed cost.
         */
        unsigned int merge(RSCM& node, DAG const& scm) const override;
    };

    /**
     * @struct MuxCountComputer
     * @brief Computes costs based on multiplexer count.
     */
    struct MuxCountComputer final : ICostComputer {
        using ICostComputer::ICostComputer;

        /**
         * @brief Computes the cost of adding a DAG to a node based on multiplexer count.
         * @param node The RSCM to compute the cost for (once merged).
         * @param scm The DAG to merge.
         * @return The computed cost.
         */
        unsigned int merge(RSCM& node, DAG const& scm) const override;
    };

    struct MuxBitsComputer final : ICostComputer {
        using ICostComputer::ICostComputer;
        unsigned int merge(RSCM& node, DAG const& scm) const override;
    };

    struct LutsCostComputer final : ICostComputer {
        using ICostComputer::ICostComputer;
        unsigned int merge(RSCM& node, DAG const& scm) const override;
    };

    struct FPGADelayComputer final : ICostComputer {
        using ICostComputer::ICostComputer;
        unsigned int merge(RSCM& node, DAG const& scm) const override;
    };

    struct ASICDelayComputer final : ICostComputer {
        using ICostComputer::ICostComputer;
        unsigned int merge(RSCM& node, DAG const& scm) const override;
    };

    /**
     * @brief Computes the bit length of a given constant.
     * @param maxConstant The constant to compute the bit length for.
     * @return The bit length.
     */
    static unsigned int BitLength(int maxConstant);
};


#endif //RSCMSOLVER_COSTCOMPUTER_H