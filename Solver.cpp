//
// Created by smith on 22/02/25.
//

#include "Solver.h"

#include <cmath>
#include <random>
#include <thread>

#include "CPModel.h"

Solver::Solver(
    std::vector<int> const& addersPerLayer,
    const int maxCoef,
    const int minCoef,
    std::vector<int> const& targets,
    const size_t nbInputBits,
    const bool useFineGrainCost
    ):
    solution(0, 0, 0, 0, 0)
{
    this->addersPerLayer = addersPerLayer;
    this->maxCoef = maxCoef;
    this->minCoef = minCoef;
    this->maxAbsCoef = std::pow(2, nbInputBits);
    this->optimalNbMuxes = optimalNbMuxes;
    this->targets = targets;
    this->nbInputBits = nbInputBits;
    this->nbAvailableThreads_ = std::thread::hardware_concurrency();
    if (useFineGrainCost) {
        fuseCostComputer = std::make_unique<FineGrainCostComputer>(this);
    }
    else {
        fuseCostComputer = std::make_unique<MuxCountComputer>(this);
    }

    // override bounding value if optimalNbMuxes is 0
    bestCost_ = std::numeric_limits<unsigned int>::max();

    paramDefs.push_back(VariableDefs::LEFT_INPUTS);
    indexToParamMap[0] = VariableDefs::LEFT_INPUTS;
    paramToIndexMap[VariableDefs::LEFT_INPUTS] = 0;

    paramDefs.push_back(VariableDefs::LEFT_SHIFTS);
    indexToParamMap[1] = VariableDefs::LEFT_SHIFTS;
    paramToIndexMap[VariableDefs::LEFT_SHIFTS] = 1;

    paramDefs.push_back(VariableDefs::RIGHT_INPUTS);
    indexToParamMap[2] = VariableDefs::RIGHT_INPUTS;
    paramToIndexMap[VariableDefs::RIGHT_INPUTS] = 2;

    paramDefs.push_back(VariableDefs::RIGHT_SHIFTS);
    indexToParamMap[3] = VariableDefs::RIGHT_SHIFTS;
    paramToIndexMap[VariableDefs::RIGHT_SHIFTS] = 3;

    paramDefs.push_back(VariableDefs::OUTPUTS_SHIFTS);
    indexToParamMap[4] = VariableDefs::OUTPUTS_SHIFTS;
    paramToIndexMap[VariableDefs::OUTPUTS_SHIFTS] = 4;

    paramDefs.push_back(VariableDefs::RIGHT_MULTIPLIER);
    indexToParamMap[5] = VariableDefs::RIGHT_MULTIPLIER;
    paramToIndexMap[VariableDefs::RIGHT_MULTIPLIER] = 5;

    // Creating layers
    int precedingLayerLastAdderNb = -1;
    int maxShift = std::ceil(std::log2(std::abs(minCoef)));
    for (int i = 0; i < addersPerLayer.size(); i++)
    {
        layers.emplace_back(i, precedingLayerLastAdderNb, addersPerLayer[i], paramDefs, maxShift);
        precedingLayerLastAdderNb = layers[i].alpha;
    }

    // Computing required number of bits in a node bitset
    nbBitsPerNode = 0;
    nbAdders = 0;
    nbPossibleMuxes = 0;
    for (auto const& layer : layers)
    {
        for (auto const& adder : layer.adders)
        {
            nbAdders++;
            if (layer.layerIdx == 0)
            {
                nbPossibleMuxes += 3;
            }
            else
            {
                nbPossibleMuxes += 5;
            }
            for (auto const& parameter : adder.variables)
            {
                nbBitsPerNode += parameter.possibleValuesFusion.size();
            }
        }
    }

    solution = RSCM(nbBitsPerNode, targets.size(), nbAdders, nbAdders * layers[0].adders[0].variables.size(), nbPossibleMuxes);
}

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
    pool.join();
}

void Solver::RunSolver(const int coef, std::atomic<int>& completedJobs, std::mutex& progressMutex, std::mutex& pushBackMutex)
{
    const CPModel cpModel(minCoef, maxCoef);
    cpModel.SolveFor(coef, scmDesigns, pushBackMutex, layers, nbBitsPerNode, paramToIndexMap, paramDefs);
    ++completedJobs;

    std::lock_guard lock(progressMutex);
    const float progress = static_cast<float>(completedJobs) / static_cast<float>(targets.size());
    std::cout << "\rProgress: " << std::setw(3) << static_cast<int>(progress * 100) << "% ["
              << std::string(static_cast<int>(progress * 50), '=') << std::string(50 - static_cast<int>(progress * 50), ' ')
              << "] " << completedJobs << '/' << targets.size() << std::flush;
}

void Solver::Solve()
{
    // 1. Find the element with design size closest and superior to 28.
    unsigned int closestIndex = 0;
    unsigned int bestDiff = UINT_MAX;

    for (int i = 0; i < scmDesigns.size(); i++) {
        const int diff = static_cast<int>(scmDesigns[i].second.size()) - static_cast<int>(nbAvailableThreads_);
        if (diff >= 0 && diff < bestDiff) {
            bestDiff = diff;
            closestIndex = i;
        }
    }

    // 2. If found, swap it with the element at index 0.
    if (closestIndex != 0) {
        std::swap(scmDesigns[0], scmDesigns[closestIndex]);
    }

    // 3. Then sort the rest (from index 1 onward) by ascending size.
    if (scmDesigns.size() > 1) {
        std::ranges::sort(scmDesigns.begin() + 1, scmDesigns.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.second.size() < rhs.second.size();
            });
    }

    // shuffle the indexes
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

    // Allocate threaded nodes
    const unsigned int nbAvailableThreads = std::thread::hardware_concurrency();
    for (int i = 0; i < nbAvailableThreads; i++) {
        threadedNodes_.emplace_back();
        for (int j = 0; j < targets.size() + 1; j++) { // Needs one more than max depth to handle solutions node in recursion
            threadedNodes_[i].emplace_back(nbBitsPerNode, targets.size(), nbAdders, nbAdders * layers[0].adders[0].variables.size(), nbPossibleMuxes); // allow for now memory reallocation during dfs
        }
    }

    std::atomic<size_t> index{0}; // index of the rscm to solve for the thread
    boost::asio::thread_pool pool(nbAvailableThreads);
    for (int thread = 0; thread < nbAvailableThreads; thread++) {
        post(pool, [this, thread, &index] {
            while (true)
            {
                constexpr int depth = 1;
                const size_t i = index.fetch_add(1, std::memory_order_relaxed);
                if (i >= scmDesigns[0].second.size()) {
                    break; // no more rscm to solve
                }
                threadedNodes_[thread][0].rscm = scmDesigns[0].second[i];
                threadedNodes_[thread][0].scmIndexes[0] = i;
                threadedNodes_[thread][0].InitializeMinShiftSavings(layers);
                ComputeBranch(depth, thread, i, 42);
            }
        });
    }
    pool.join();
}

void Solver::ComputeBranch(const int depth, const int threadNb, const unsigned int startIndex, unsigned int currentCost) // NOLINT(*-no-recursion)
{
    if (depth == targets.size())
    {
        if (currentCost < bestCost_.load(std::memory_order_relaxed))
        {
            std::lock_guard lock(solutionMutex_); // high performance cost, done at last possible moment
            if (currentCost < bestCost_)
            {
                bestCost_ = currentCost;
                solution = threadedNodes_[threadNb][depth - 1];
                // PrintSolution(currentCost);
            }
        }
        return;
    }

    for (const auto& index : threadedIndexes_[startIndex][depth - 1])
    {
        threadedNodes_[threadNb][depth] = threadedNodes_[threadNb][depth - 1];
        currentCost = fuseCostComputer->compute(threadedNodes_[threadNb][depth], scmDesigns[depth].second[index]);
        // const unsigned int muxCost = fst;
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

unsigned int Solver::MuxCountComputer::compute(RSCM& node, DAG const& scm) const
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

unsigned int Solver::FineGrainCostComputer::compute(RSCM& node, DAG const& scm) const
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

    const auto leftShiftBase  = solver->paramToIndexMap.at(VariableDefs::LEFT_SHIFTS);
    const auto rightShiftBase = solver->paramToIndexMap.at(VariableDefs::RIGHT_SHIFTS);
    const size_t mapSize = solver->paramToIndexMap.size();

    size_t bitPos = 0;
    unsigned int adderIdx = 0;
    unsigned int paramGlobalIdx = 0;

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

                // 3) Count bits in one pass
                unsigned bitCount = 0;
                for (size_t j = 0; j < param.possibleValuesFusion.size(); ++j) {
                    if (node.rscm.set.test(bitPos + j)) ++bitCount;
                }
                bitPos += param.possibleValuesFusion.size();

                // 4) If more than one bit, add fine-grain mux cost
                if (bitCount > 1 &&
                    solver->indexToParamMap.at(paramInAdderIdx) != VariableDefs::RIGHT_MULTIPLIER)
                {
                    fineGrainCost += computeMuxCost(bitCount, paramGlobalIdx);
                }

                // 5) Handle RIGHT_MULTIPLIER adder cost
                if (solver->indexToParamMap.at(paramInAdderIdx) == VariableDefs::RIGHT_MULTIPLIER) {
                    // determine coefficient
                    unsigned coeff = 67u;
                    if      (node.isPlusMinus[adderIdx])    coeff = 93u;
                    else if (node.rscm.isMinus[adderIdx])   coeff = 75u;

                    // compute bit-width post-shift
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
                std::cout << paramDefsToString(indexToParamMap[p]) << " : ";
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

void Solver::PrettyPrinter(RSCM & solutionNode, const unsigned int cost)
{
    // copy solutionNode.rscm to a new node
    DAG rcmCopy = solutionNode.rscm;
    rcmCopy.set = boost::dynamic_bitset<>(rcmCopy.set.count());
    const unsigned int muxCost = MuxCountComputer(this).compute(solutionNode, rcmCopy);
    const unsigned int fineGrainCost = FineGrainCostComputer(this).compute(solutionNode, rcmCopy);

    std::cout << "SOLUTION COST: " << muxCost << " muxes, "  << fineGrainCost << " fine grain" << std::endl << std::endl;

    std::cout << std::endl << "Parameters: " << std::endl;
    const auto paramDefsToString = VarDefsToString();
    unsigned int bitIndex = 0;
    for (const auto& layer : layers) {
        std::cout << "Layer " << layer.layerIdx << ":" << std::endl;
        for (const auto& adder : layer.adders) {
            std::cout << "\tAdder " << adder.adderIdx << ":" << std::endl;
            for (int p = 0; p < adder.variables.size(); p++) {
                std::cout << "\t\t" << paramDefsToString(indexToParamMap[p]) << " : { ";
                for (const auto& v : adder.variables[p].possibleValuesFusion) {
                    if (solutionNode.rscm.set.test(bitIndex + v + adder.variables[p].zeroPoint)) {
                        if (indexToParamMap[p] == VariableDefs::LEFT_INPUTS || indexToParamMap[p] == VariableDefs::RIGHT_INPUTS)
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
                        else if (indexToParamMap[p] == VariableDefs::RIGHT_MULTIPLIER)
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
        std::cout << paramDefsToString(indexToParamMap.at(mux % layers[0].adders[0].variables.size())) << " of adder " << mux / layers[0].adders[0].variables.size() << " | ";
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

void Solver::RscmToCsv(std::string const& fileUri)
{
    std::ofstream csvFile(fileUri);
    if (!csvFile.is_open())
    {
        std::cerr << "Error opening file: " << fileUri << std::endl;
        return;
    }

    const auto paramDefsToString = VarDefsToString();
    csvFile << "Layer" << ",";
    csvFile << "Adder" << ",";
    for (const auto& parameter : paramDefs)
    {
        csvFile << paramDefsToString(parameter) << ",";
    }

    unsigned int bitIndex = 0;
    for (const auto& layer : layers) {
        for (const auto& adder : layer.adders) {
            csvFile << std::endl << layer.layerIdx << "," << adder.adderIdx << ",";
            for (int p = 0; p < adder.variables.size(); p++) {
                csvFile << "{ ";
                for (const auto& v : adder.variables[p].possibleValuesFusion) {
                    if (solution.rscm.set.test(bitIndex + v + adder.variables[p].zeroPoint)) {
                        if (indexToParamMap[p] == VariableDefs::LEFT_INPUTS || indexToParamMap[p] == VariableDefs::RIGHT_INPUTS)
                        {
                            if (v < 0)
                            {
                                if (v == -1)
                                {
                                    csvFile << "X" << " ";
                                } else
                                {
                                    csvFile << std::abs(v+1) << "X" << " ";
                                }
                            }
                            else
                            {
                                csvFile << "Adder" << v << " ";
                            }
                        }
                        else if (indexToParamMap[p] == VariableDefs::RIGHT_MULTIPLIER)
                        {
                            if (v == -1)
                            {
                                csvFile << "-" << " ";
                            }
                            else
                            {
                                csvFile << "+" << " ";
                            }
                        }
                        else
                        {
                            csvFile << v << " ";
                        }
                    }
                }
                bitIndex += adder.variables[p].possibleValuesFusion.size();
                csvFile << "}" << ",";
            }
        }
    }
    csvFile << std::endl;
}
