//
// Created by smith on 22/02/25.
//

#include "Solver.h"

#include <cmath>
#include <algorithm>
#include <random>
#include <ranges>
#include <stdexcept>
#include <thread>
#include <optional>
#include <unordered_map>
#include <string>
#include <chrono>

#include "../cp/CPModel.h"
#include "../writers/JSONDumper.h"
#include "../writers/Verilog.h"
#include "CostComputer.h"

// Constructor initializes solver parameters, layout configuration, and variable definitions
Solver::Solver(
    std::vector<int> const& layout,       // Defines number of adders per layer
    const int maxCoef,                    // Maximum coefficient allowed
    const int minCoef,                    // Minimum coefficient allowed
    std::vector<int> const& targets,      // Target constants to synthesize
    const size_t nbInputBits,             // Input bit-width
    const CostModel costModel,            // Cost model selection
    const size_t lutWidth,                // LUT width for LUT-based cost models
    const bool isSymmetric                // symmetric mux support (-T)
) : solution(0, 0, 0, 0, 0)               // Initialize empty solution node
{
    this->layout = layout;
    this->maxCoef = maxCoef;
    this->minCoef = minCoef;
    this->maxAbsCoef = std::pow(2, nbInputBits);
    this->targets = targets;
    this->nbInputBits = nbInputBits;
    this->nbAvailableThreads_ = std::thread::hardware_concurrency();
    this->costModel_ = costModel;
    this->lutWidth_ = lutWidth;
    this->isSymmetric_ = isSymmetric;
    this->nbMinEncodingBits_ = std::ceil(std::log2(targets.size()));

    // Find largest m such that m + ceil(log2(m)) <= lutWidth
    // At minimum, a LUT can implement a 2:1 mux.
    this->maxMuxInputPerLut_ = 2;
    for (unsigned m = 2; m <= lutWidth; ++m) { // safe upper cap
        const auto s = static_cast<unsigned>(std::ceil(std::log2(static_cast<double>(m))));
        if (m + s <= lutWidth) this->maxMuxInputPerLut_ = m;
        else break; // m grows => m+s grows, so we can stop
    }

    this->nbInputsLeftByAdderLut_ = lutWidth_ - 2; // add/sub LUT only has two inputs even if doing both

    // Select cost model implementation
    switch (costModel_) {
    case CostModel::AreaCost:
        fuseCostComputer = std::make_unique<CostComputer::AreaCostComputer>(this);
        break;
    case CostModel::MuxCount:
        fuseCostComputer = std::make_unique<CostComputer::MuxCountComputer>(this);
        break;
    case CostModel::MuxBits:
        fuseCostComputer = std::make_unique<CostComputer::MuxBitsComputer>(this);
        break;
    case CostModel::LutsCost:
        fuseCostComputer = std::make_unique<CostComputer::LutsCostComputer>(this);
        break;
    case CostModel::FPGADelay:
        fuseCostComputer = std::make_unique<CostComputer::FPGADelayComputer>(this);
        break;
    case CostModel::ASICDelay:
        fuseCostComputer = std::make_unique<CostComputer::ASICDelayComputer>(this);
        break;
    default:
        throw std::invalid_argument("Unknown cost model");
    }

    // Normalize targets constants
    normShift_ = std::numeric_limits<unsigned int>::max();
    for (const auto& t : targets) {
        if (t != 0) {
            const auto tz = static_cast<unsigned int>(__builtin_ctz(static_cast<unsigned int>(std::abs(t))));
            normShift_ = tz > normShift_ ? normShift_ : tz;
        }
    }
    if (normShift_ > 0) {
        std::cout << "Normalizing targets by right-shifting " << normShift_ << " bits." << std::endl;
        for (auto& t : this->targets) {
            t = t >> normShift_;
        }
    }

    // Initial cost bound (infinite)
    bestCost_ = std::numeric_limits<unsigned int>::max();

    // Map internal variable types to and from indexes
    varDefs.push_back(VariableDefs::LEFT_INPUTS);
    idxToVarMap[0] = VariableDefs::LEFT_INPUTS;
    varToIdxMap[VariableDefs::LEFT_INPUTS] = 0;

    varDefs.push_back(VariableDefs::LEFT_SHIFTS);
    idxToVarMap[1] = VariableDefs::LEFT_SHIFTS;
    varToIdxMap[VariableDefs::LEFT_SHIFTS] = 1;

    varDefs.push_back(VariableDefs::RIGHT_INPUTS);
    idxToVarMap[2] = VariableDefs::RIGHT_INPUTS;
    varToIdxMap[VariableDefs::RIGHT_INPUTS] = 2;

    varDefs.push_back(VariableDefs::RIGHT_SHIFTS);
    idxToVarMap[3] = VariableDefs::RIGHT_SHIFTS;
    varToIdxMap[VariableDefs::RIGHT_SHIFTS] = 3;

    varDefs.push_back(VariableDefs::OUTPUTS_SHIFTS);
    idxToVarMap[4] = VariableDefs::OUTPUTS_SHIFTS;
    varToIdxMap[VariableDefs::OUTPUTS_SHIFTS] = 4;

    varDefs.push_back(VariableDefs::LEFT_MULTIPLIER);
    idxToVarMap[5] = VariableDefs::LEFT_MULTIPLIER;
    varToIdxMap[VariableDefs::LEFT_MULTIPLIER] = 5;

    varDefs.push_back(VariableDefs::RIGHT_MULTIPLIER);
    idxToVarMap[6] = VariableDefs::RIGHT_MULTIPLIER;
    varToIdxMap[VariableDefs::RIGHT_MULTIPLIER] = 6;

    // Build each layer with specified number of adders
    int precedingLayerLastAdderNb = -1;
    int maxShift = std::ceil(std::log2(std::abs(minCoef)));
    for (int i = 0; i < layout.size(); i++) {
        layers.emplace_back(i, precedingLayerLastAdderNb, layout[i], varDefs, maxShift);
        precedingLayerLastAdderNb = layers[i].alpha;
    }

    // Calculate bit requirements for storing SCM configuration
    nbBitsPerSCM = 0;
    nbAdders = 0;
    nbPossibleVariables = 0;
    for (auto const& layer : layers) {
        for (auto const& adder : layer.adders) {
            nbAdders++;
            if (layer.layerIdx == 0) {
                nbPossibleVariables += 3; // initial layer uses only subset of parameters
            } else {
                nbPossibleVariables += 5;
            }
            for (auto const& parameter : adder.variables) {
                nbBitsPerSCM += parameter.possibleValuesFusion.size();
            }
        }
    }

    // Create an empty solution RSCM node with all required size fields
    solution = RSCM(nbBitsPerSCM, targets.size(), nbAdders,
        nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);

    if (costModel_ == CostModel::MuxCount || costModel_ == CostModel::MuxBits) {
        signAgnosticMask_.resize(nbBitsPerSCM);
        size_t bitPos = 0;
        for (const auto& layer : layers) {
            for (const auto& adder : layer.adders) {
                for (const auto& paramDef : varDefs) {
                    const auto& param = adder.variables[varToIdxMap.at(paramDef)];
                    const size_t sz = param.possibleValuesFusion.size();
                    if (paramDef != VariableDefs::LEFT_MULTIPLIER &&
                        paramDef != VariableDefs::RIGHT_MULTIPLIER) {
                        for (size_t i = 0; i < sz; ++i) {
                            signAgnosticMask_.set(bitPos + i);
                        }
                    }
                    bitPos += sz;
                }
            }
        }
    }
}

// Solve each target constant using CP (constraint programming)
void Solver::CPSolve()
{
    CPSolve(std::nullopt);
}

void Solver::CPSolve(std::optional<unsigned int> heuristic)
{
    // Start timeout countdown at beginning of solving if configured
    if (bbTimeoutSeconds_.has_value() && !bbDeadline_.has_value()) {
        bbDeadline_ = std::chrono::steady_clock::now() + std::chrono::seconds(*bbTimeoutSeconds_);
    }

    std::atomic completedJobs(0);
    boost::asio::thread_pool pool(std::thread::hardware_concurrency());

    for (auto const& t : targets) {
        auto coef = static_cast<short>(t);
        post(pool, [this, coef, &completedJobs, heuristic] {
            RunSolver(coef, completedJobs, progressMutex_, pushBackMutex_, heuristic);
        });
    }
    pool.join(); // Wait for all threads to complete
}

void Solver::SetBranchTimeoutSeconds(const std::optional<unsigned int> seconds)
{
    bbTimeoutSeconds_ = seconds.has_value() && seconds.value() > 0 ? seconds : std::optional<unsigned int>{};
}

// Internal call to CP solver for a single target coefficient
void Solver::RunSolver(const int coef, std::atomic<int>& completedJobs,
    std::mutex& progressMutex, std::mutex& pushBackMutex, const std::optional<unsigned int> heuristic)
{
    const CPModel cpModel(
        minCoef,
        maxCoef,
        static_cast<int>(-1 * std::pow(2, nbInputBits - 1)),
        static_cast<int>(std::pow(2, nbInputBits - 1) - 1)
    );
    cpModel.SolveFor(coef, scmDesigns, pushBackMutex, layers, nbBitsPerSCM,
                     varToIdxMap, varDefs, heuristic);

    if (costModel_ == CostModel::MuxCount || costModel_ == CostModel::MuxBits) {
        std::lock_guard lock(pushBackMutex);
        auto it = std::find_if(scmDesigns.begin(), scmDesigns.end(),
            [coef](const auto& entry) { return entry.first == coef; });
        if (it != scmDesigns.end()) {
            auto& scms = it->second;
            std::unordered_set<std::string> seen;
            seen.reserve(scms.size());
            std::vector<DAG> dedup;
            dedup.reserve(scms.size());
            for (const auto& scm : scms) {
                boost::dynamic_bitset<> signature = scm.set & signAgnosticMask_;
                std::string key;
                boost::to_string(signature, key);
                if (seen.insert(key).second) {
                    dedup.push_back(scm);
                }
            }
            scms.swap(dedup);
        }
    }
    ++completedJobs;

    // Progress bar output
    std::lock_guard lock(progressMutex);
    const float progress = static_cast<float>(completedJobs) / static_cast<float>(targets.size());
    std::cout << "\rProgress: " << std::setw(3) << static_cast<int>(progress * 100) << "% ["
              << std::string(static_cast<int>(progress * 50), '=') << std::string(50 - static_cast<int>(progress * 50), ' ')
              << "] " << completedJobs << '/' << targets.size() << std::flush;
}

// Final stage of solving: compute best combination of SCMs across all targets
void Solver::Solve()
{
    timeoutTriggered_.store(false, std::memory_order_relaxed);

    // 1. Select the target with the number of SCMs closest in size to number of available threads
    // (this is usefull to balance the load across threads)
    unsigned int closestIndex = 0;
    unsigned int bestDiff = UINT_MAX;
    for (int i = 0; i < scmDesigns.size(); i++) {
        const int diff = static_cast<int>(scmDesigns[i].second.size()) - static_cast<int>(nbAvailableThreads_);
        if (diff >= 0 && diff < bestDiff) {
            bestDiff = diff;
            closestIndex = i;
        }
    }

    // 2. Move best-fit target to index 0
    if (closestIndex != 0) {
        std::swap(scmDesigns[0], scmDesigns[closestIndex]);
    }

    // 3. Sort remaining SCMs by size (ascending) (a heuristic that appears to accelerate convergence)
    if (scmDesigns.size() > 1) {
        std::ranges::sort(scmDesigns.begin() + 1, scmDesigns.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.second.size() < rhs.second.size();
            });
    }

    // Shuffle SCMs for each target constants to accelerate convergence
    // As the SCMs yielded after each other by the CP solver are "close" to each other
    // Shuffling them will help the B&P to explore the search space more evenly
    threadedIndexes_.resize(scmDesigns[0].second.size());
    for (int i = 0; i < scmDesigns[0].second.size(); i++) {
        for (int j = 1; j < scmDesigns.size(); j++) {
            auto temp = std::vector<unsigned int>(scmDesigns[j].second.size());
            for (int k = 0; k < scmDesigns[j].second.size(); k++) {
                temp[k] = k;
            }
            std::ranges::shuffle(temp, std::mt19937(std::random_device()()));
            threadedIndexes_[i].push_back(temp);
        }
    }

    // Allocate memory for each thread’s solution path
    const unsigned int nbAvailableThreads = std::thread::hardware_concurrency();
    for (int i = 0; i < nbAvailableThreads; i++) {
        threadedNodes_.emplace_back();
        for (int j = 0; j < targets.size() + 1; j++) {
            threadedNodes_[i].emplace_back(nbBitsPerSCM, targets.size(), nbAdders,
                nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
        }
    }

    // Launch branch and prune threads
    // The first target SCMs are split across threads
    std::atomic<size_t> index{0};
    boost::asio::thread_pool pool(nbAvailableThreads);
    for (int thread = 0; thread < nbAvailableThreads; thread++) {
        post(pool, [this, thread, &index] {
            while (true) {
                if (timeoutTriggered_.load(std::memory_order_relaxed)) break;
                constexpr int depth = 1;
                const size_t i = index.fetch_add(1, std::memory_order_relaxed);
                if (i >= scmDesigns[0].second.size()) {
                    break;
                }
                threadedNodes_[thread][0].rscm = scmDesigns[0].second[i];
                threadedNodes_[thread][0].scmIndexes[0] = i;
                threadedNodes_[thread][0].InitializeMinShiftSavings(layers);
                ComputeBranch(depth, thread, i, 42);
            }
        });
    }
    pool.join(); // Wait for all branches to finish

    // Reshift solution SCMs to original scale
    if (normShift_ > 0) {
        for (auto& t : this->targets) {
            t = t << normShift_;
        }
        for (auto& fst : this->scmDesigns | std::views::keys) {
            fst = fst << normShift_;
        }
        ApplyNormalizationShift(solution);
    }

    if (timeoutTriggered_.load(std::memory_order_relaxed)) {
        std::cout << "Branch-and-bound terminated due to timeout; best-so-far solution is not guaranteed optimal." << std::endl;
    }

    // Compute all cost models once for the final solution
    solutionCosts_ = GetAllCosts(solution);

    if (costModel_ != CostModel::AreaCost)
    {
        // more complicated, we have to replay the whole merging process to compute the fine-grained cost
        RSCM replayNode(nbBitsPerSCM, targets.size(), nbAdders,
            nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
        // Start by copying the first SCM directly (not merging), like the native execution does
        replayNode.rscm = scmDesigns[0].second[solution.scmIndexes[0]];
        replayNode.scmIndexes = solution.scmIndexes;
        replayNode.InitializeMinShiftSavings(layers);
        const CostComputer::AreaCostComputer fineGrainCostComputer(this);
        // Now merge the remaining SCMs starting from depth 1
        for (int depth = 1; depth < targets.size(); depth++)
        {
            fineGrainCostComputer.merge(replayNode, scmDesigns[depth].second[solution.scmIndexes[depth]]);
        }
        ApplyNormalizationShift(replayNode);
        solution = replayNode;
    }
}

// Recursive branch and prune search to find minimum-cost solution
void Solver::ComputeBranch(const int depth, const int threadNb, const unsigned int startIndex, unsigned int currentCost)
{
    // If we have reached the last layer, check if the cost is lower than the best found
    if (depth == targets.size()) {
        if (currentCost < bestCost_.load(std::memory_order_relaxed)) {
            std::lock_guard lock(solutionMutex_);
            if (currentCost < bestCost_) { // update only if still the best
                bestCost_ = currentCost;
                solution = threadedNodes_[threadNb][depth - 1];
                PrintSolution(currentCost);
            }
        }
        return;
    }

    if (timeoutTriggered_.load(std::memory_order_relaxed) ||
        (bbDeadline_.has_value() && std::chrono::steady_clock::now() >= *bbDeadline_)) {
        timeoutTriggered_.store(true, std::memory_order_relaxed);
        return;
        }

    // Recurse to next layer if cost is below current best
    for (const auto& index : threadedIndexes_[startIndex][depth - 1]) {
        threadedNodes_[threadNb][depth] = threadedNodes_[threadNb][depth - 1];
        currentCost = fuseCostComputer->merge(threadedNodes_[threadNb][depth], scmDesigns[depth].second[index]);
        if (currentCost >= bestCost_.load(std::memory_order_relaxed))
        {
            continue;
        }
        threadedNodes_[threadNb][depth].scmIndexes[depth] = index;
        ComputeBranch(depth + 1, threadNb, startIndex, currentCost);
    }
}

bool Solver::IsPowerOfTwo(const int x) {
    return x > 0 && (x & x - 1) == 0;
}

unsigned int Solver::BitLength(const int maxConstant)
{
    return static_cast<unsigned int>(std::ceil(std::log2(std::max(std::abs(maxConstant + 1), std::abs(maxConstant))))) + 1;
}

void Solver::ApplyNormalizationShift(RSCM& node) const
{
    if (normShift_ == 0) return;

    size_t bitPos = 0;
    unsigned int adderIdx = 0;

    for (const auto& layer : layers) {
        for (const auto& adder : layer.adders) {
            unsigned int paramIdx = 0;
            for (const auto& param : adder.variables) {
                if (adderIdx == nbAdders - 1 && idxToVarMap.at(paramIdx) == VariableDefs::OUTPUTS_SHIFTS) {
                    const size_t baseBit = bitPos + param.zeroPoint;
                    std::vector<int> selectedValues;

                    // Collect selected shifts first to avoid clobbering when we rewrite bits
                    for (const int v : param.possibleValuesFusion) {
                        const size_t bitIndex = baseBit + static_cast<size_t>(v);
                        if (node.rscm.set.test(bitIndex)) {
                            selectedValues.push_back(v);
                        }
                        node.rscm.set.reset(bitIndex);
                    }

                    const int maxShift = *std::ranges::max_element(param.possibleValuesFusion
                    );

                    for (const int v : selectedValues) {
                        const int shiftedValue = v + static_cast<int>(normShift_);
                        if (shiftedValue > maxShift) continue; // out of domain, skip

                        const size_t newBitIndex = baseBit + static_cast<size_t>(shiftedValue);
                        node.rscm.set.set(newBitIndex);
                    }

                    const size_t globalParamIdx =
                        varToIdxMap.at(VariableDefs::OUTPUTS_SHIFTS) + adderIdx * adder.variables.size();
                    const int scale = 1 << normShift_;
                    node.rscm.maxOutputValue[globalParamIdx] *= scale;
                    node.rscm.minOutputValue[globalParamIdx] *= scale;
                    node.variableBitWidths[globalParamIdx] = std::max(
                        BitLength(node.rscm.maxOutputValue[globalParamIdx]),
                        BitLength(node.rscm.minOutputValue[globalParamIdx])
                    );
                    return;
                }
                bitPos += param.possibleValuesFusion.size();
                ++paramIdx;
            }
            ++adderIdx;
        }
    }
}

std::optional<unsigned int> Solver::EvaluateCost(const RSCM& solutionNode, const CostModel model) const
{
    auto* solverPtr = const_cast<Solver*>(this);
    try {
        switch (model) {
        case CostModel::AreaCost: {
            // Recompute fine-grain cost by replaying merges
            RSCM replayNode(nbBitsPerSCM, targets.size(), nbAdders,
                nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
            replayNode.rscm = scmDesigns[0].second[solutionNode.scmIndexes[0]];
            replayNode.InitializeMinShiftSavings(layers);
            const CostComputer::AreaCostComputer areaCostComputer(solverPtr);
            for (int depth = 1; depth < targets.size(); depth++) {
                areaCostComputer.merge(replayNode, scmDesigns[depth].second[solutionNode.scmIndexes[depth]]);
            }
            ApplyNormalizationShift(replayNode);
            return replayNode.cost;
        }
        case CostModel::MuxCount: {
            DAG emptySCM(nbBitsPerSCM, nbAdders, nbPossibleVariables);
            emptySCM.set.reset();
            RSCM replayNode = solutionNode;
            const CostComputer::MuxCountComputer muxCostComputer(solverPtr);
            ApplyNormalizationShift(replayNode);
            return muxCostComputer.merge(replayNode, emptySCM);
        }
        case CostModel::MuxBits: {
            DAG emptySCM(nbBitsPerSCM, nbAdders, nbPossibleVariables);
            emptySCM.set.reset();
            RSCM replayNode = solutionNode;
            const CostComputer::MuxBitsComputer muxBitsComputer(solverPtr);
            ApplyNormalizationShift(replayNode);
            return muxBitsComputer.merge(replayNode, emptySCM);
        }
        case CostModel::LutsCost: {
            // Recompute fine-grain cost by replaying merges
            RSCM replayNode(nbBitsPerSCM, targets.size(), nbAdders,
                nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
            replayNode.rscm = scmDesigns[0].second[solutionNode.scmIndexes[0]];
            replayNode.InitializeMinShiftSavings(layers);
            const CostComputer::LutsCostComputer lutsCostComputer(solverPtr);
            for (int depth = 1; depth < targets.size(); depth++) {
                lutsCostComputer.merge(replayNode, scmDesigns[depth].second[solutionNode.scmIndexes[depth]]);
            }
            ApplyNormalizationShift(replayNode);
            return replayNode.cost;
        }
        case CostModel::FPGADelay: {
            // Recompute fine-grain cost by replaying merges
            RSCM replayNode(nbBitsPerSCM, targets.size(), nbAdders,
                nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
            replayNode.rscm = scmDesigns[0].second[solutionNode.scmIndexes[0]];
            replayNode.InitializeMinShiftSavings(layers);
            const CostComputer::FPGADelayComputer fpgaDelayComputer(solverPtr);
            for (int depth = 1; depth < targets.size(); depth++) {
                fpgaDelayComputer.merge(replayNode, scmDesigns[depth].second[solutionNode.scmIndexes[depth]]);
            }
            ApplyNormalizationShift(replayNode);
            return replayNode.cost;
        }
        case CostModel::ASICDelay: {
            RSCM replayNode = solutionNode;
            const CostComputer::ASICDelayComputer asicDelayComputer(solverPtr);
            return asicDelayComputer.merge(replayNode, scmDesigns[0].second[solutionNode.scmIndexes[0]]);
        }
        default:
            return std::nullopt;
        }
    } catch (const std::logic_error&) {
        return std::nullopt;
    }
}

std::unordered_map<std::string, std::optional<unsigned int>> Solver::GetAllCosts(const RSCM& solutionNode) const
{
    const std::vector<std::pair<CostModel, std::string>> models = {
        {CostModel::AreaCost,   "area_cost"},
        {CostModel::MuxCount,   "mux_count"},
        {CostModel::MuxBits,    "mux_bits"},
        {CostModel::LutsCost,   "luts"},
        {CostModel::FPGADelay,  "fpga_delay"},
        {CostModel::ASICDelay,  "asic_delay"},
    };

    std::unordered_map<std::string, std::optional<unsigned int>> costs;
    for (const auto& [model, name] : models) {
        costs[name] = EvaluateCost(solutionNode, model);
    }
    return costs;
}

void Solver::PrintSolution(unsigned int const cost)
{
    const auto paramDefsToString = VarDefsToString();
    std::cout << "Solution cost: " << cost << " : ";
    unsigned int bitIndex = 0;
    for (const auto& layer : layers) {
        std::cout << "Layer " << layer.layerIdx << " : ";
        for (const auto& adder : layer.adders) {
            std::cout << "Adder " << adder.adderIdx << " : ";
            for (int p = 0; p < adder.variables.size(); p++) {
                std::cout << paramDefsToString(idxToVarMap[p]) << " : ";
                for (const auto& v : adder.variables[p].possibleValuesFusion) {
                    if (solution.rscm.set.test(bitIndex + v + adder.variables[p].zeroPoint)) {
                        std::cout << v << " ";
                    }
                }
                bitIndex += adder.variables[p].possibleValuesFusion.size();
            }
        }
    }
    std::cout << std::endl;
}

void Solver::PrettyPrinter(const RSCM& solutionNode)
{
    std::cout << "Costs:\n";
    for (const auto& [name, value] : solutionCosts_) {
        std::cout << " - " << name << ": ";
        if (value.has_value()) {
            std::cout << *value;
        } else {
            std::cout << "not implemented";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << std::endl << "Parameters: " << std::endl;
    const auto paramDefsToString = VarDefsToString();
    unsigned int bitIndex = 0;
    for (const auto& layer : layers) {
        std::cout << "Layer " << layer.layerIdx << ":" << std::endl;
        for (const auto& adder : layer.adders) {
            std::cout << "\tAdder " << adder.adderIdx << ":" << std::endl;
            for (int p = 0; p < adder.variables.size(); p++) {
                std::cout << "\t\t" << paramDefsToString(idxToVarMap[p]) << " : { ";
                for (const auto& v : adder.variables[p].possibleValuesFusion) {
                    if (solutionNode.rscm.set.test(bitIndex + v + adder.variables[p].zeroPoint)) {
                        if (idxToVarMap[p] == VariableDefs::LEFT_INPUTS || idxToVarMap[p] == VariableDefs::RIGHT_INPUTS)
                        {
                            if (v < 0)
                            {
                                if (v == -1)
                                {
                                    std::cout << "X" << " ";
                                } else
                                {
                                    std::cout << std::abs(v+1) << "X" << " ";
                                }
                            }
                            else
                            {
                                std::cout << "Adder" << v << " ";
                            }
                        }
                        else if (idxToVarMap[p] == VariableDefs::RIGHT_MULTIPLIER ||
                                 idxToVarMap[p] == VariableDefs::LEFT_MULTIPLIER)
                        {
                            if (v == -1)
                            {
                                std::cout << "-" << " ";
                            }
                            else
                            {
                                std::cout << "+" << " ";
                            }
                        }
                        else
                        {
                            std::cout << v << " ";
                        }
                    }
                }
                bitIndex += adder.variables[p].possibleValuesFusion.size();
                std::cout << "}" << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

void Solver::Verilog(const RSCM& solutionNode, const std::string& outputUri, const bool overwrite) const
{
    VerilogGenerator verilog(solutionNode, outputUri, layers, idxToVarMap, varToIdxMap, nbInputBits, targets, scmDesigns, overwrite);
}

void Solver::DumpJSON(const RSCM& solutionNode, const std::string& outputUri, const bool overwrite) const
{
    const auto costs = GetAllCosts(solutionNode);
    JSONDumper JSONDumper(solution, outputUri, layers, idxToVarMap, varToIdxMap, nbInputBits, targets, scmDesigns, costs, isSymmetric_, overwrite);
}

void Solver::DumpSnapshot(const RSCM& solutionNode, const std::string& outputUri, const bool overwrite) const
{
    WriteSnapshot(*this, solutionNode, outputUri, overwrite);
}

std::unordered_map<std::string, std::optional<unsigned int>> Solver::ComputeAllCosts(const RSCM& solutionNode) const
{
    return GetAllCosts(solutionNode);
}

void Solver::SolveConfigToMuxMapping() const
{
    std::cout << "------------------" << std::endl;

    const auto paramDefsToString = VarDefsToString();
    unsigned int bitIndex = 0;
    unsigned int parameterIndex = 0;
    std::vector<std::vector<unsigned int>> muxEnconding(targets.size());
    std::vector<unsigned int> muxNames;
    for (const auto& layer : layers)
    {
        for (const auto& adder : layer.adders) {
            for (const auto & parameter : adder.variables) {
                unsigned int nbPossibleValues = 0;
                // get number of possible values for this parameter
                for (size_t v = 0; v < parameter.possibleValuesFusion.size(); v++) {
                    if (solution.rscm.set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                        nbPossibleValues++;
                    }
                }
                // if more than one bit is set then we have a multiplexer
                if (nbPossibleValues > 1) {
                    unsigned int valueNumber = 0;
                    muxNames.push_back(parameterIndex);
                    for (size_t v = 0; v < parameter.possibleValuesFusion.size(); v++) {
                        if (solution.rscm.set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint)) {
                            for (int targetIndex = 0; targetIndex < targets.size(); targetIndex++)
                            {
                                if (scmDesigns[targetIndex].second[solution.scmIndexes[targetIndex]].set.test(bitIndex + parameter.possibleValuesFusion[v] + parameter.zeroPoint))
                                {
                                    muxEnconding[targetIndex].push_back(valueNumber);
                                }
                            }
                            valueNumber++;
                        }
                    }
                }
                bitIndex += parameter.possibleValuesFusion.size();
                parameterIndex++;
            }
        }
    }

    // print mux names
    std::cout << "Muxes : ";
    for (const auto& mux : muxNames)
    {
        std::cout << paramDefsToString(idxToVarMap.at(mux % layers[0].adders[0].variables.size())) << " of adder " << mux / layers[0].adders[0].variables.size() << " | ";
    }
    std::cout << std::endl;

    // print mux encoding
    for (int i = 0; i < targets.size(); i++)
    {
        std::cout << "Target " << scmDesigns[i].first << "\t : \t";
        for (const auto& mux : muxEnconding[i])
        {
            std::cout << mux << " ";
        }
        std::cout << std::endl;
    }
}

unsigned int Solver::GetCurrentCost() const
{
    return bestCost_.load();
}
