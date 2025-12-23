//
// Created by smith on 22/02/25.
//

#include "Solver.h"

#include <cmath>
#include <algorithm>
#include <random>
#include <ranges>
#include <thread>

#include "CPModel.h"
#include "JSONDumper.h"
#include "Verilog.h"

// Constructor initializes solver parameters, layout configuration, and variable definitions
Solver::Solver(
    std::vector<int> const& layout,       // Defines number of adders per layer
    const int maxCoef,                    // Maximum coefficient allowed
    const int minCoef,                    // Minimum coefficient allowed
    std::vector<int> const& targets,      // Target constants to synthesize
    const size_t nbInputBits,             // Input bit-width
    const bool useFineGrainCost           // Flag to choose between cost models
) : solution(0, 0, 0, 0, 0)               // Initialize empty solution node
{
    this->layout = layout;
    this->maxCoef = maxCoef;
    this->minCoef = minCoef;
    this->maxAbsCoef = std::pow(2, nbInputBits);
    this->targets = targets;
    this->nbInputBits = nbInputBits;
    this->nbAvailableThreads_ = std::thread::hardware_concurrency();
    this->useFineGrainCost = useFineGrainCost;

    // Select cost model implementation
    if (useFineGrainCost) {
        fuseCostComputer = std::make_unique<FineGrainCostComputer>(this);
    } else {
        fuseCostComputer = std::make_unique<MuxCountComputer>(this);
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

    varDefs.push_back(VariableDefs::RIGHT_MULTIPLIER);
    idxToVarMap[5] = VariableDefs::RIGHT_MULTIPLIER;
    varToIdxMap[VariableDefs::RIGHT_MULTIPLIER] = 5;

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
}

// Solve each target constant using CP (constraint programming)
void Solver::CPSolve()
{
    std::atomic completedJobs(0);
    boost::asio::thread_pool pool(std::thread::hardware_concurrency());

    for (auto const& t : targets) {
        auto coef = static_cast<short>(t);
        post(pool, [this, coef, &completedJobs] {
            RunSolver(coef, completedJobs, progressMutex_, pushBackMutex_);
        });
    }
    pool.join(); // Wait for all threads to complete
}

// Internal call to CP solver for a single target coefficient
void Solver::RunSolver(const int coef, std::atomic<int>& completedJobs,
    std::mutex& progressMutex, std::mutex& pushBackMutex)
{
    const CPModel cpModel(
        minCoef,
        maxCoef,
        static_cast<int>(-1 * std::pow(2, nbInputBits - 1)),
        static_cast<int>(std::pow(2, nbInputBits - 1) - 1)
    );
    cpModel.SolveFor(coef, scmDesigns, pushBackMutex, layers, nbBitsPerSCM,
                     varToIdxMap, varDefs);
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
        ApplyNormalizationShift(solution, true);
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

unsigned int Solver::BitLength(const int value)
{
    return static_cast<unsigned int>(std::ceil(std::log2(std::max(std::abs(value + 1), std::abs(value))))) + 1;
}

void Solver::ApplyNormalizationShift(RSCM& node, const bool logShifts) const
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

                    const int maxShift = *std::max_element(
                        param.possibleValuesFusion.begin(),
                        param.possibleValuesFusion.end()
                    );

                    for (const int v : selectedValues) {
                        const int shiftedValue = v + static_cast<int>(normShift_);
                        if (shiftedValue > maxShift) continue; // out of domain, skip

                        const size_t newBitIndex = baseBit + static_cast<size_t>(shiftedValue);
                        node.rscm.set.set(newBitIndex);
                        if (logShifts) {
                            std::cout << v << " shifted to " << shiftedValue << std::endl;
                        }
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

unsigned int Solver::MuxCountComputer::merge(RSCM& node, DAG const& scm) const
{
    // 1) Merge resource sets and update bounds
    node.rscm.set |= scm.set;

    // Prepare to count
    size_t bitPos = 0;
    unsigned int totalParams = 0;
    unsigned int muxCount = 0;

    // Iterate layers → adders → variables
    for (const auto& layer : solver->layers) {
        for (const auto& adder : layer.adders) {
            for (const auto& param : adder.variables) {
                // Count how many bits in this parameter are set
                unsigned int bitsSet = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (node.rscm.set.test(bitPos + j)) {
                        ++bitsSet;
                    }
                }
                bitPos += param.possibleValuesFusion.size();

                // If more than one bit, we need (bitsSet - 1) multiplexers
                if (bitsSet > 1)
                {
                    // compute the number of bits to select instead
                    muxCount += std::ceil(std::log2(bitsSet));
                }

                ++totalParams;
            }
        }
    }

    // Store result
    node.cost = muxCount;
    return muxCount;
}

unsigned int Solver::FineGrainCostComputer::merge(RSCM& node, DAG const& scm) const
{
    unsigned int fineGrainCost = 0;

    // 1) Merge resource sets and update bounds
    node.rscm.set |= scm.set;
    std::transform(
        node.rscm.maxOutputValue.begin(),
        node.rscm.maxOutputValue.end(),
        scm.maxOutputValue.begin(),
        node.rscm.maxOutputValue.begin(),
        [](const int a, const int b) { return std::max(a, b); }
    ); // max output values
    std::transform(
        node.rscm.minOutputValue.begin(),
        node.rscm.minOutputValue.end(),
        scm.minOutputValue.begin(),
        node.rscm.minOutputValue.begin(),
        [](const int a, const int b) { return std::min(a, b); }
    ); // min output values

    // update bitwidths
    for (size_t i = 0; i < node.variableBitWidths.size(); ++i) {
        node.variableBitWidths[i] = std::max(BitLength(node.rscm.maxOutputValue[i]), BitLength(node.rscm.minOutputValue[i]));
    }

    const auto leftShiftBase  = solver->varToIdxMap.at(VariableDefs::LEFT_SHIFTS);
    const auto rightShiftBase = solver->varToIdxMap.at(VariableDefs::RIGHT_SHIFTS);
    const size_t mapSize = solver->varToIdxMap.size();

    auto updateMinShift = [&](const unsigned idx) {
        // Use coefficient's trailing zeros, not the multiplied values'
        // This is correct because any input in the range could have 0 trailing zeros
        const unsigned coeffTZ = scm.coefficientTrailingZeros[idx];
        node.minShiftSavings[idx] = std::min(node.minShiftSavings[idx], coeffTZ);
    };

    auto computeMuxCost = [&](const unsigned bitsCount, const unsigned idx) {
        const unsigned int nbBits = node.variableBitWidths[idx] - node.minShiftSavings[idx];
        return 14u * bitsCount * nbBits;
    };

    size_t bitPos = 0;
    unsigned int adderIdx = 0;
    unsigned int paramGlobalIdx = 0;

    // Iterate layers → adders → variables
    for (const auto& layer : solver->layers) {
        for (const auto& adder : layer.adders) {
            // Track plus-minus flag
            node.isPlusMinus[adderIdx] =
                node.isPlusMinus[adderIdx] || node.rscm.isMinus[adderIdx] != scm.isMinus[adderIdx];

            unsigned paramInAdderIdx = 0;
            for (const auto& param : adder.variables) {
                // 2) Update shift savings early
                updateMinShift(paramGlobalIdx);

                // 3) Count the number of bits at one for this variable
                unsigned bitCount = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (node.rscm.set.test(bitPos + j)) ++bitCount;
                }
                bitPos += param.possibleValuesFusion.size();

                // 4) If more than one bit and not an adder -> we have a multiplexer
                if (bitCount > 1 && solver->idxToVarMap.at(paramInAdderIdx) != VariableDefs::RIGHT_MULTIPLIER)
                {
                    fineGrainCost += computeMuxCost(bitCount, paramGlobalIdx);
                }

                // 5) Handle RIGHT_MULTIPLIER adder cost (an adder)
                if (solver->idxToVarMap.at(paramInAdderIdx) == VariableDefs::RIGHT_MULTIPLIER) {
                    // determine coefficient depending on the adder type
                    unsigned coeff = 67u;
                    if      (node.isPlusMinus[adderIdx])    coeff = 93u;
                    else if (node.rscm.isMinus[adderIdx])   coeff = 75u;

                    // compute the minimum shift savings (number of trailing zeros)
                    const auto leftIdx  = leftShiftBase  + adderIdx * mapSize;
                    const auto rightIdx = rightShiftBase + adderIdx * mapSize;
                    const unsigned minShift   = std::min(
                        node.minShiftSavings[leftIdx],
                        node.minShiftSavings[rightIdx]
                    );

                    // compute FA/HA costs
                    if (node.rscm.minOutputValue[rightIdx] != 0 && node.rscm.maxOutputValue[rightIdx] != 0) {
                        unsigned wa = node.variableBitWidths[leftIdx] - minShift;
                        unsigned wb = node.variableBitWidths[rightIdx] - minShift;
                        unsigned int shiftA = node.minShiftSavings[leftIdx] - minShift;
                        unsigned int shiftB = node.minShiftSavings[rightIdx] - minShift;
                        auto [fst, snd] = std::minmax(shiftA, shiftB);
                        const unsigned int diff = snd - fst;

                        unsigned fa = 0, ha = 0;
                        if (!node.isPlusMinus[adderIdx] && !node.rscm.isMinus[adderIdx]) {
                            if (node.rscm.minOutputValue[leftIdx] != 0 && node.rscm.maxOutputValue[leftIdx] != 0) {
                                fa = std::max(wa, wb) - diff - 1;
                                ha = 1;
                            }
                        } else {
                            if (node.rscm.minOutputValue[leftIdx] == 0 && node.rscm.maxOutputValue[leftIdx] == 0) { // left (a) is zero
                                ha = wb - node.minShiftSavings[rightIdx];
                            } else if (shiftB >= shiftA) {
                                fa = std::max(wa, wb) - diff;
                            } else if (node.minShiftSavings[leftIdx] >= wb) {
                                fa = std::max(wa, wb) - diff - 1;
                                ha = 1 + wb;
                            } else {
                                fa = std::max(wa, wb) - diff;
                                ha = diff;
                            }
                        }
                        fineGrainCost += static_cast<unsigned>(coeff * fa + 5.0/9.0 * coeff * ha);
                    }
                }

                ++paramGlobalIdx;
                ++paramInAdderIdx;
            }
            ++adderIdx;
        }
    }

    // Store result
    node.cost = fineGrainCost;
    return fineGrainCost;
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
    unsigned int nbMuxes;
    unsigned int fineGrainedCost;
    // Need to compute the cost of the cost function that was not used during the solving
    if (useFineGrainCost)
    {
        // compute mux cost by creating an empty DAG scm to merge to the solutionNode
        DAG emptySCM(nbBitsPerSCM, nbAdders, nbPossibleVariables);
        emptySCM.set.reset();
        const MuxCountComputer muxCostComputer(this);
        RSCM replayNode = solutionNode; // copy the solution node to avoid modifying it
        nbMuxes = muxCostComputer.merge(replayNode, emptySCM);
        fineGrainedCost = solutionNode.cost;
    } else
    {
        // more complicated, we have to replay the whole merging process to compute the fine-grained cost
        RSCM replayNode(nbBitsPerSCM, targets.size(), nbAdders,
            nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
        // Start by copying the first SCM directly (not merging), like the native execution does
        replayNode.rscm = scmDesigns[0].second[solutionNode.scmIndexes[0]];
        replayNode.InitializeMinShiftSavings(layers);
        const FineGrainCostComputer fineGrainCostComputer(this);
        // Now merge the remaining SCMs starting from depth 1
        for (int depth = 1; depth < targets.size(); depth++)
        {
            fineGrainCostComputer.merge(replayNode, scmDesigns[depth].second[solutionNode.scmIndexes[depth]]);
        }
        fineGrainedCost = replayNode.cost;
        nbMuxes = solutionNode.cost;
    }

    std::cout << "Mux COST: " << nbMuxes << std::endl;
    std::cout << "FINE-GRAINED COST: " << fineGrainedCost << std::endl << std::endl;
    
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
                        else if (idxToVarMap[p] == VariableDefs::RIGHT_MULTIPLIER)
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

void Solver::Verilog(const RSCM& solutionNode, const std::string& outputUri, const bool overwrite)
{
    if (!useFineGrainCost)
    {
        // more complicated, we have to replay the whole merging process to compute the fine-grained cost
        RSCM replayNode(nbBitsPerSCM, targets.size(), nbAdders,
            nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
        // Start by copying the first SCM directly (not merging), like the native execution does
        replayNode.rscm = scmDesigns[0].second[solutionNode.scmIndexes[0]];
        replayNode.InitializeMinShiftSavings(layers);
        const FineGrainCostComputer fineGrainCostComputer(this);
        // Now merge the remaining SCMs starting from depth 1
        for (int depth = 1; depth < targets.size(); depth++)
        {
            fineGrainCostComputer.merge(replayNode, scmDesigns[depth].second[solutionNode.scmIndexes[depth]]);
        }
        ApplyNormalizationShift(replayNode);
        VerilogGenerator verilog(replayNode, outputUri, layers, idxToVarMap, varToIdxMap, nbInputBits, targets, scmDesigns, overwrite);
    } else {
        VerilogGenerator verilog(solutionNode, outputUri, layers, idxToVarMap, varToIdxMap, nbInputBits, targets, scmDesigns, overwrite);
    }
}

void Solver::DumpJSON(const RSCM& solutionNode, const std::string& outputUri, const bool overwrite)
{
    if (!useFineGrainCost)
    {
        // more complicated, we have to replay the whole merging process to compute the fine-grained cost
        RSCM replayNode(nbBitsPerSCM, targets.size(), nbAdders,
            nbAdders * layers[0].adders[0].variables.size(), nbPossibleVariables);
        // Start by copying the first SCM directly (not merging), like the native execution does
        replayNode.rscm = scmDesigns[0].second[solutionNode.scmIndexes[0]];
        replayNode.scmIndexes = solutionNode.scmIndexes;
        replayNode.InitializeMinShiftSavings(layers);
        const FineGrainCostComputer fineGrainCostComputer(this);
        // Now merge the remaining SCMs starting from depth 1
        for (int depth = 1; depth < targets.size(); depth++)
        {
            fineGrainCostComputer.merge(replayNode, scmDesigns[depth].second[solutionNode.scmIndexes[depth]]);
        }
        ApplyNormalizationShift(replayNode);
        JSONDumper JSONDumper(replayNode, outputUri, layers, idxToVarMap, varToIdxMap, nbInputBits, targets, scmDesigns, overwrite);
    } else {
        JSONDumper JSONDumper(solutionNode, outputUri, layers, idxToVarMap, varToIdxMap, nbInputBits, targets, scmDesigns, overwrite);
    }
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
