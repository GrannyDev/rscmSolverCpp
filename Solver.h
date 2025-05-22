//
// Created by smith on 22/02/25.
//

#ifndef PROBLEM_H
#define PROBLEM_H

#include <unordered_set>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <fstream>

#include "Layer.h"
#include "RSCM.h"
#include "DAG.h"

/**
 * @class Solver
 * @brief The class finding the best RSCM for a given set of target constants.
 */
class Solver {
public:
    /**
     * @brief Constructor for the Solver class.
     * @param layout Vector specifying where the adders and layers are: {2,1} means 2 adders in layer 1 and 1 in layer 2.
     * @param maxCoef Maximum (intermediate) coefficient value.
     * @param minCoef Minimum (intermediate) coefficient value.
     * @param targets The set of target constants.
     * @param nbInputBits Number of bits encoding the input (X).
     * @param useFineGrainCost Flag to enable fine-grained cost function instead of coarse-grained.
     */
    Solver(
        std::vector<int> const& layout,
        int maxCoef, int minCoef,
        std::vector<int> const& targets,
        size_t nbInputBits,
        bool useFineGrainCost = false
    );

    /**
     * @brief Step 1: Find all possible SCMs for the given target constants.
     */
    void CPSolve();

    /**
     * @brief Step 2: DFS Branch & Prune to find the best RSCM.
     */
    void Solve();

    /**
     * @brief Prints the solution in a human-readable format.
     * @param solutionNode The solution RSCM to print.
     */
    void PrettyPrinter(const RSCM& solutionNode);

    /**
     * @brief Gets the current best cost of the solution.
     * @return The current best cost.
     */
    unsigned int GetCurrentCost() const;

    /**
     * @brief Retrieves the configuration bits assignment for each target constant.
     */
    void SolveConfigToMuxMapping() const;

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

    /**
     * @struct FineGrainCostComputer
     * @brief Computes costs with fine-grain granularity.
     */
    struct FineGrainCostComputer final : ICostComputer {
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

    std::unique_ptr<ICostComputer> fuseCostComputer; ///< Pointer to the cost computation strategy.
    std::vector<int> layout; ///< Layout for the RSCM.
    int maxCoef; ///< Maximum (intermediate) coefficient value.
    int minCoef; ///< Minimum (intermediate) coefficient value.
    unsigned int maxAbsCoef; ///< Maximum absolute (intermediate) coefficient value.
    size_t nbBitsPerSCM; ///< Number of bits to encode an SCM / RSCM.
    size_t nbAdders; ///< Number of adders in the layout.
    size_t nbPossibleVariables; ///< Number of variables in the layout.
    size_t nbInputBits; ///< Number of input bits of X.
    std::vector<VariableDefs> varDefs; ///< Variable definitions.
    std::vector<int> targets; ///< Target values.
    std::vector<Layer> layers; ///< Layers in the layout.
    std::vector<std::pair<int, std::vector<DAG>>> scmDesigns; ///< SCM designs for each target.
    RSCM solution; ///< The solution RSCM.
    std::unordered_map<VariableDefs, unsigned int> varToIdxMap; ///< Map from variables to indices.
    std::unordered_map<unsigned int, VariableDefs> idxToVarMap; ///< Map from indices to variables.

private:
    /**
     * @brief Runs the solver for a specific target coefficient.
     * @param coef The coefficient to find all SCM for.
     * @param completedJobs Atomic counter for completed jobs.
     * @param progressMutex Mutex for progress updates.
     * @param pushBackMutex Mutex for thread-safe push into the SCMi set.
     */
    void RunSolver(int coef, std::atomic<int> & completedJobs, std::mutex & progressMutex, std::mutex & pushBackMutex);

    /**
     * @brief Computes a branch of the merge tree.
     * @param depth The current depth in the tree.
     * @param threadNb The thread number.
     * @param startIndex The starting index in the SCM list for Ti.
     * @param currentCost The current cost of the RSCM.
     */
    void ComputeBranch(int depth, int threadNb, unsigned int startIndex, unsigned int currentCost);

    /**
     * @brief Prints the solution with the given cost.
     * @param cost The cost of the solution.
     */
    void PrintSolution(unsigned int cost);

    /**
     * @brief Compares two integers in two's complement representation.
     * @param a The first integer.
     * @param b The second integer.
     * @return Comparison result.
     */
    static int CompareTwoComplement(int a, int b);

    /**
     * @brief Checks if a number is a power of two.
     * @param x The number to check.
     * @return True if the number is a power of two, false otherwise.
     */
    static bool IsPowerOfTwo(int x);

    /**
     * @brief Computes the bit length of a given constant.
     * @param maxConstant The constant to compute the bit length for.
     * @return The bit length.
     */
    unsigned int BitLength(int maxConstant) const;

    std::mutex pushBackMutex_; ///< Mutex for thread-safe operations.
    std::mutex progressMutex_; ///< Mutex for progress updates.
    std::mutex solutionMutex_; ///< Mutex for solution updates.
    std::atomic<unsigned int> bestCost_; ///< Atomic variable when updating the best cost.
    std::vector<std::vector<RSCM>> threadedNodes_; ///< Threaded nodes for parallel merging.
    unsigned int nbAvailableThreads_; ///< Number of available threads.
    std::vector<std::vector<std::vector<unsigned int>>> threadedIndexes_; ///< Threaded indexes for parallel computation.
    /// contains the shuffled indexes of the SCMs for each target constant
};

#endif //PROBLEM_H