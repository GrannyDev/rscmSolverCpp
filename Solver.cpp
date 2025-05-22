//
// Created by smith on 22/02/25.
//

#include "Solver.h"

#include <cmath>
#include <random>
#include <thread>

#include "CPModel.h"

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

    // Select cost model implementation
    if (useFineGrainCost) {
        fuseCostComputer = std::make_unique<FineGrainCostComputer>(this);
    } else {
        fuseCostComputer = std::make_unique<MuxCountComputer>(this);
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
    const CPModel cpModel(minCoef, maxCoef);
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

int Solver::CompareTwoComplement(int a, int b) {
    int absA = std::abs(a);
    int absB = std::abs(b);

    if (absA != absB)
        return absA > absB ? a : b;
    // When abs values are equal, check if the number is a power of two.
    // If so, the positive number requires one more bit.
    // (Note: 0 is a special case, but typically both representations are equivalent.)
    if (IsPowerOfTwo(absA))
        return a >= 0 ? a : b;
    // Otherwise, they need the same number of bits.
    return a;
    // or b—they're equivalent in bit width.
}

bool Solver::IsPowerOfTwo(const int x) {
    return x > 0 && (x & x - 1) == 0;
}

unsigned int Solver::BitLength(const int maxConstant) const
{
    if (maxConstant == 0)
    {
        return 0;
    }

    const bool isNegative = maxConstant < 0;
    auto bitWidth = static_cast<unsigned int>(std::ceil(std::log2(std::abs(maxConstant) * maxAbsCoef)));

    if (IsPowerOfTwo(maxConstant) && isNegative)
    {
        bitWidth += 1;
    }

    return bitWidth;
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
                    muxCount += bitsSet - 1;
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
        CompareTwoComplement
    );

    const auto leftShiftBase  = solver->varToIdxMap.at(VariableDefs::LEFT_SHIFTS);
    const auto rightShiftBase = solver->varToIdxMap.at(VariableDefs::RIGHT_SHIFTS);
    const size_t mapSize = solver->varToIdxMap.size();

    auto updateMinShift = [&](const unsigned idx) {
        const auto maxVal = scm.maxOutputValue[idx];
        const unsigned tz = maxVal ? __builtin_ctz(maxVal) : std::numeric_limits<unsigned>::max();
        node.minShiftSavings[idx] = std::min(node.minShiftSavings[idx], tz);
    };

    auto computeMuxCost = [&](const unsigned bitsCount, const unsigned idx) {
        const int constVal = node.rscm.maxOutputValue[idx] >> node.minShiftSavings[idx];
        const unsigned nbBits = solver->BitLength(constVal);
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
                if (bitCount > 1 &&
                    solver->idxToVarMap.at(paramInAdderIdx) != VariableDefs::RIGHT_MULTIPLIER)
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
                    const int a = node.rscm.maxOutputValue[leftIdx];
                    const int b = node.rscm.maxOutputValue[rightIdx];
                    if (b != 0) {
                        unsigned wa = solver->BitLength(a) - minShift;
                        unsigned wb = solver->BitLength(b) - minShift;
                        unsigned int shiftA = node.minShiftSavings[leftIdx] - minShift;
                        unsigned int shiftB = node.minShiftSavings[rightIdx] - minShift;
                        auto [fst, snd] = std::minmax(shiftA, shiftB);
                        const unsigned int diff = snd - fst;

                        unsigned fa = 0, ha = 0;
                        if (!node.isPlusMinus[adderIdx] && !node.rscm.isMinus[adderIdx]) {
                            if (a != 0) {
                                fa = std::max(wa, wb) - diff - 1;
                                ha = 1;
                            }
                        } else {
                            if (a == 0) {
                                ha = wb - static_cast<unsigned>(__builtin_ctz(b));
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
    std::cout << "SOLUTION COST: " << solutionNode.cost << std::endl << std::endl;

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
