//
// Created by smith on 22/02/25.
//
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <optional>

#include "solver/Solver.h"
#include "writers/SnapshotIO.h"

int main(int argc, char* argv[]) {
    auto parseList = [](const std::string& s) {
        std::vector<int> vals;
        std::stringstream ss(s);
        std::string token;
        while (std::getline(ss, token, ',')) {
            if (token.empty()) continue;
            vals.push_back(std::stoi(token));
        }
        return vals;
    };

    auto parseCostModel = [](const std::string& s) -> std::optional<CostModel> {
        if (s == "area" || s == "area_cost") return CostModel::AreaCost;
        if (s == "mux_count") return CostModel::MuxCount;
        if (s == "mux_bits") return CostModel::MuxBits;
        if (s == "luts") return CostModel::LutsCost;
        if (s == "fpga_delay") return CostModel::FPGADelay;
        if (s == "asic_delay") return CostModel::ASICDelay;
        return std::nullopt;
    };

    size_t beta = 8; // maximum bit-width allowed for the constants (and intermediate values)
    size_t nbInputBits = 5; // number of bits of the input (to compute the fine-grained cost function)
    std::vector<int> targets = {-23, -22, -16, -15, -14, -13, -9, -8, -7, -6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 25, 26}; // target const set of the RSCM
    std::vector<int> layout = {1, 1}; // {1,1} describes the chosen layout, i.e. 1 adder on the first layer and 1 adder on the second layer
    std::optional<unsigned int> heuristic;
    std::optional<unsigned int> timeoutSeconds = 300;
    size_t lutWidth = 6;
    bool doJsonDump = true;
    std::string jsonPath = "dump.json";
    std::optional<std::string> snapshotOut;
    std::optional<std::string> snapshotIn;
    auto costModel = CostModel::LutsCost;
    bool isSymmetric = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout <<
                "Usage: " << argv[0] << " [options]\n"
                "Options:\n"
                "  --beta=<uint>            Set beta (max bit-width for constants), default 6\n"
                "  --nb-input-bits=<uint>   Set number of input bits, default 8\n"
                "  --targets=a,b,c          Comma-separated target constants\n"
                "  --layout=x,y             Comma-separated adders per layer\n"
                "  --heuristic=<uint>       Limit CP solutions per target\n"
                "  --cost=<model>           Cost model: area|mux_count|mux_bits|luts|fpga_delay|asic_delay\n"
                "  --lut-width=<uint>       Set LUT width for LUT-based cost models, default 6\n"
                "  --timeout=<uint>         Branch-and-bound timeout in seconds (CP phase not affected)\n"
                "  --json=<path>            Enable JSON dump to path (default dump.json)\n"
                "  --no-json                Disable JSON dump\n"
                "  --snapshot-out=<path>    Write a snapshot to recompute costs later\n"
                "  --recompute-snapshot=<path>  Recompute all costs from a snapshot file\n"
                "  -h, --help               Show this help\n";
            return 0;
        }
        try {
            if (arg.rfind("--beta=", 0) == 0) {
                beta = std::stoul(arg.substr(7));
            } else if (arg.rfind("--nb-input-bits=", 0) == 0) {
                nbInputBits = std::stoul(arg.substr(16));
            } else if (arg.rfind("--targets=", 0) == 0) {
                targets = parseList(arg.substr(10));
            } else if (arg.rfind("--layout=", 0) == 0) {
                layout = parseList(arg.substr(9));
            } else if (arg.rfind("--heuristic=", 0) == 0) {
                heuristic = static_cast<unsigned int>(std::stoul(arg.substr(12)));
            } else if (arg.rfind("--timeout=", 0) == 0) {
                // timeoutSeconds = static_cast<unsigned int>(std::stoul(arg.substr(10)));
            } else if (arg.rfind("--lut-width=", 0) == 0) {
                lutWidth = std::stoul(arg.substr(12));
            } else if (arg == "--no-json") {
                doJsonDump = false;
            } else if (arg.rfind("--json=", 0) == 0) {
                doJsonDump = true;
                jsonPath = arg.substr(7);
            } else if (arg.rfind("--cost=", 0) == 0) {
                auto parsed = parseCostModel(arg.substr(7));
                if (parsed.has_value()) {
                    costModel = *parsed;
                } else {
                    std::cerr << "Unknown cost model: " << arg.substr(7) << std::endl;
                }
            } else if (arg.rfind("--snapshot-out=", 0) == 0) {
                snapshotOut = arg.substr(15);
            } else if (arg.rfind("--recompute-snapshot=", 0) == 0) {
                snapshotIn = arg.substr(21);
            } else if (arg == "--is_symmetric") {
                isSymmetric = true;
            } else {
                std::cerr << "Ignoring unknown argument: " << arg << std::endl;
            }
        } catch (const std::exception&) {
            std::cerr << "Failed to parse argument: " << arg << std::endl;
        }
    }

    const auto maxCoef = static_cast<int>(std::pow(2, beta - 1) - 1);
    const auto minCoef = static_cast<int>(-1 * std::pow(2, beta - 1));

    if (snapshotIn.has_value()) {
        const auto snap = ReadSnapshot(*snapshotIn);
        if (!snap.has_value()) {
            std::cerr << "Failed to read snapshot: " << *snapshotIn << std::endl;
            return 1;
        }
        const auto& s = *snap;
        const bool sym = isSymmetric || s.isSymmetric;
        Solver problem(s.layout, s.maxCoef, s.minCoef, s.targets, s.nbInputBits, costModel, s.lutWidth, sym);
        problem.scmDesigns.clear();
        for (size_t i = 0; i < s.targets.size(); ++i) {
            problem.scmDesigns.emplace_back(s.targets[i], std::vector{ s.selectedScms[i] });
        }
        problem.solution = RSCM(s.nbBitsPerSCM, s.targets.size(), s.nbAdders, s.nbVariables, s.nbMuxes);
        problem.solution.scmIndexes = std::vector<unsigned int>(s.targets.size(), 0); // only one SCM per target in snapshot
        problem.solution.rscm.set = s.rscmSet;
        problem.solution.rscm.maxOutputValue = s.rscmMaxOutput;
        problem.solution.rscm.minOutputValue = s.rscmMinOutput;
        problem.solution.rscm.coefficientTrailingZeros = s.rscmCoeffTZ;
        problem.solution.rscm.isMinus = s.rscmIsMinus;
        problem.solution.minShiftSavings = s.minShiftSavings;
        problem.solution.variableBitWidths = s.variableBitWidths;
        problem.solution.isPlusMinus = s.isPlusMinus;
        problem.solution.InitializeMinShiftSavings(problem.layers);
        problem.solutionCosts_ = problem.ComputeAllCosts(problem.solution);
        std::cout << "_____________________RECOMPUTED COSTS_____________________" << std::endl;
        problem.PrettyPrinter(problem.solution);
        if (doJsonDump) {
            problem.DumpJSON(problem.solution, jsonPath, true);
        }
        return 0;
    }
    Solver problem(layout, maxCoef, minCoef, targets, nbInputBits, costModel, lutWidth, isSymmetric);

    problem.SetBranchTimeoutSeconds(timeoutSeconds);
    problem.CPSolve(heuristic); // step 1: solve the problem with the CPSolver

    // print the numbers of SCMs found for each target constant
    std::cout << std::endl;
    for (const auto& [fst, snd] : problem.scmDesigns) {
        std::cout << "Target: " << fst << " solutions: " << snd.size() << " | ";
    }
    std::cout << std::endl;

    // step 2: solve the problem with the DFS B&P
    problem.Solve();

    std::cout << "_____________________BEST SOLUTION_____________________" << std::endl;
    problem.PrettyPrinter(problem.solution);
    problem.SolveConfigToMuxMapping();
    if (doJsonDump) {
        problem.DumpJSON(problem.solution, jsonPath, true);
    }
    if (snapshotOut.has_value()) {
        problem.DumpSnapshot(problem.solution, *snapshotOut, true);
    }

    return 0;
}
