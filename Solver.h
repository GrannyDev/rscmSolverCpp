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
#include "Node.h"
#include "Rcm.h"

class Solver {
public:
    Solver(
        std::vector<int> const& addersPerLayer,
        int maxCoef, int minCoef,
        std::vector<int> const& targets,
        size_t nbInputBits,
        bool useFineGrainCost = false
        );
    void CPSolve();
    void Solve();
    void PrettyPrinter(Node & solutionNode, unsigned int cost);
    unsigned int GetCurrentCost() const;
    void RscmToCsv(std::string const& fileUri);
    void SolveConfigToMuxMapping() const;

    struct ICostComputer {
        explicit ICostComputer(Solver* solver) : solver(solver) {}
        virtual ~ICostComputer() = default;
        virtual unsigned int compute(Node& threadNode, Rcm const& scmToAdd) const = 0;
        Solver* solver;
    };

    struct FineGrainCostComputer final : ICostComputer {
        using ICostComputer::ICostComputer;
        unsigned int compute(Node& node, Rcm const& scm) const override;
    };

    struct MuxCountComputer final : ICostComputer {
        using ICostComputer::ICostComputer;
        unsigned int compute(Node& node, Rcm const& scm) const override;
    };

    std::unique_ptr<ICostComputer> fuseCostComputer;
    std::vector<int> addersPerLayer;
    int maxCoef;
    int minCoef;
    unsigned int maxAbsCoef;
    size_t nbBitsPerNode;
    size_t optimalNbMuxes;
    size_t nbAdders;
    size_t nbPossibleMuxes;
    size_t nbInputBits;
    bool forceSpaceExploration;
    std::vector<ParamDefs> paramDefs;
    std::vector<int> targets;
    std::vector<Layer> layers;
    std::vector<std::pair<int, std::vector<Rcm>>> scmDesigns;
    Node solution;
    std::unordered_map<ParamDefs, unsigned int> paramToIndexMap;
    std::unordered_map<unsigned int, ParamDefs> indexToParamMap;

private:
    void RunSolver(int coef, std::atomic<int> & completedJobs, std::mutex & progressMutex, std::mutex & pushBackMutex);
    void ComputeBranch(int depth, int threadNb, unsigned int startIndex, unsigned int currentCost);
    void PrintSolution(unsigned int cost);
    static int CompareTwoComplement(int a, int b);
    static bool IsPowerOfTwo(int x);
    unsigned int BitLength(int maxConstant) const;

    std::mutex pushBackMutex_;
    std::mutex progressMutex_;
    std::mutex solutionMutex_;
    std::atomic<unsigned int> bestCost_;
    std::vector<std::vector<Node>> threadedNodes_;
    unsigned int nbAvailableThreads_;
    std::vector<std::vector<std::vector<unsigned int>>> threadedIndexes_;
};



#endif //PROBLEM_H
