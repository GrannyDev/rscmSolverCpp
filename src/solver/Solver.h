//
// Created by smith on 22/02/25.
//

#ifndef PROBLEM_H
#define PROBLEM_H

#include <unordered_set>
#include <unordered_map>
#include <optional>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <fstream>

#include "../shared/Layer.h"
#include "../shared/DAG.h"
#include "RSCM.h"
#include "CostModel.h"
#include "ICostComputer.h"
#include "../writers/SnapshotIO.h"

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
     * @param costModel Cost model to use during solving.
     * @param lutWidth LUT width for LUT-based cost models
     */
    Solver(
        std::vector<int> const& layout,
        int maxCoef, int minCoef,
        std::vector<int> const& targets,
        size_t nbInputBits,
        CostModel costModel = CostModel::MuxCount,
        size_t lutWidth = 6
    );

    /**
     * @brief Step 1: Find all possible SCMs for the given target constants.
     */
    void CPSolve();
    void CPSolve(std::optional<unsigned int> heuristic);

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
     * @brief Writes the solution in Verilog format.
     * @param solutionNode The RSCM solution node to convert to Verilog.
     * @param outputUri The file to write to.
     * @param overwrite If the existing file is overwritten, will throw an error if false and the file exists
     */
    void Verilog(
        const RSCM& solutionNode,
        const std::string& outputUri,
        bool overwrite
    ) const;

    void DumpJSON(
        const RSCM& solutionNode,
        const std::string& outputUri,
        bool overwrite
    ) const;

    void DumpSnapshot(
        const RSCM& solutionNode,
        const std::string& outputUri,
        bool overwrite
    ) const;

    std::unordered_map<std::string, std::optional<unsigned int>> ComputeAllCosts(const RSCM& solutionNode) const;

    /**
     * @brief Gets the current best cost of the solution.
     * @return The current best cost.
     */
    unsigned int GetCurrentCost() const;

    /**
     * @brief Retrieves the configuration bits assignment for each target constant.
     */
    void SolveConfigToMuxMapping() const;

    std::unique_ptr<ICostComputer> fuseCostComputer; ///< Pointer to the cost computation strategy.
    std::vector<int> layout; ///< Layout for the RSCM.
    int maxCoef; ///< Maximum (intermediate) coefficient value.
    int minCoef; ///< Minimum (intermediate) coefficient value.
    unsigned int maxAbsCoef; ///< Maximum absolute (intermediate) coefficient value.
    size_t nbBitsPerSCM; ///< Number of bits to encode an SCM / RSCM.
    size_t nbAdders; ///< Number of adders in the layout.
    size_t nbPossibleVariables; ///< Number of variables in the layout.
    size_t nbInputBits; ///< Number of input bits of X.
    size_t lutWidth_; ///< LUT width for LUT-based cost models
    size_t nbMinEncodingBits_; ///< Number of bits required to encode each target constant
    size_t maxMuxInputPerLut_; ///< The maximum number of mux inputs that can be implemented in a single LUT.
    size_t nbInputsLeftByAdderLut_; ///< What's the biggest mux LUT we can merge into an adder LUT.
    CostModel costModel_; ///< Selected cost model.
    std::vector<VariableDefs> varDefs; ///< Variable definitions.
    std::vector<int> targets; ///< Target values.
    std::vector<Layer> layers; ///< Layers in the layout.
    std::vector<std::pair<int, std::vector<DAG>>> scmDesigns; ///< SCM designs for each target.
    RSCM solution; ///< The solution RSCM.
    std::unordered_map<std::string, std::optional<unsigned int>> solutionCosts_; ///< Cached costs for the best solution.
    std::unordered_map<VariableDefs, unsigned int> varToIdxMap; ///< Map from variables to indices.
    std::unordered_map<unsigned int, VariableDefs> idxToVarMap; ///< Map from indices to variables.

private:
    /**
     * @brief Runs the solver for a specific target coefficient.
     * @param coef The coefficient to find all SCM for.
     * @param completedJobs Atomic counter for completed jobs.
     * @param progressMutex Mutex for progress updates.
     * @param pushBackMutex Mutex for thread-safe push into the SCMi set.
     * @param heuristic Optional limit on number of SCMs to enumerate.
     */
    void RunSolver(int coef, std::atomic<int> & completedJobs, std::mutex & progressMutex, std::mutex & pushBackMutex, std::optional<unsigned int> heuristic = std::nullopt);

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
    static unsigned int BitLength(int maxConstant) ;

    /**
     * @brief Apply the normalization shift to the output shifts of the last adder in the given node.
     * @param node The RSCM node to update.
     */
    void ApplyNormalizationShift(RSCM& node) const;

    /**
     * @brief Compute cost for a given model using replay nodes when needed.
     * @param solutionNode Node to evaluate.
     * @param model Cost model to use.
     * @return Optional cost (empty if not implemented).
     */
    std::optional<unsigned int> EvaluateCost(const RSCM& solutionNode, CostModel model) const;

    /**
     * @brief Compute all costs for implemented models.
     * @param solutionNode Node to evaluate.
     * @return Map of model name to optional cost.
     */
    std::unordered_map<std::string, std::optional<unsigned int>> GetAllCosts(const RSCM& solutionNode) const;

    std::mutex pushBackMutex_; ///< Mutex for thread-safe operations.
    std::mutex progressMutex_; ///< Mutex for progress updates.
    std::mutex solutionMutex_; ///< Mutex for solution updates.
    std::atomic<unsigned int> bestCost_; ///< Atomic variable when updating the best cost.
    std::vector<std::vector<RSCM>> threadedNodes_; ///< Threaded nodes for parallel merging.
    unsigned int nbAvailableThreads_; ///< Number of available threads.
    unsigned int normShift_; ///< Normalization shift for the SCMs.
    std::vector<std::vector<std::vector<unsigned int>>> threadedIndexes_; ///< Threaded indexes for parallel computation.
    /// contains the shuffled indexes of the SCMs for each target constant
};

#endif //PROBLEM_H
