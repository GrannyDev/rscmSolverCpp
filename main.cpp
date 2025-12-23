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

#include "Solver.h"

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

    auto parseCostModel = [](const std::string& s) -> std::optional<Solver::CostModel> {
        if (s == "area" || s == "area_cost") return Solver::CostModel::AreaCost;
        if (s == "mux_count") return Solver::CostModel::MuxCount;
        if (s == "mux_bits") return Solver::CostModel::MuxBits;
        if (s == "luts") return Solver::CostModel::LutsCost;
        if (s == "fpga_delay") return Solver::CostModel::FPGADelay;
        if (s == "asic_delay") return Solver::CostModel::ASICDelay;
        return std::nullopt;
    };

    size_t beta = 6; // maximum bit-width allowed for the constants (and intermediate values)
    size_t nbInputBits = 8; // number of bits of the input (to compute the fine-grained cost function)
    std::vector<int> targets = {-20, -13, -8, -6, -5, -3, -2, -1, 0, 1, 2, 4, 5, 7, 12, 19}; // target const set of the RSCM
    std::vector<int> layout = {1, 1}; // {1,1} describes the chosen layout, i.e. 1 adder on the first layer and 1 adder on the second layer
    std::optional<unsigned int> heuristic;
    bool doJsonDump = true;
    std::string jsonPath = "dump.json";
    Solver::CostModel costModel = Solver::CostModel::MuxBits;

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
                "  --json=<path>            Enable JSON dump to path (default dump.json)\n"
                "  --no-json                Disable JSON dump\n"
                "  -h, --help               Show this help\n";
            return 0;
        }
        try {
            if (arg.rfind("--beta=", 0) == 0) {
                beta = static_cast<size_t>(std::stoul(arg.substr(7)));
            } else if (arg.rfind("--nb-input-bits=", 0) == 0) {
                nbInputBits = static_cast<size_t>(std::stoul(arg.substr(16)));
            } else if (arg.rfind("--targets=", 0) == 0) {
                targets = parseList(arg.substr(10));
            } else if (arg.rfind("--layout=", 0) == 0) {
                layout = parseList(arg.substr(9));
            } else if (arg.rfind("--heuristic=", 0) == 0) {
                heuristic = static_cast<unsigned int>(std::stoul(arg.substr(12)));
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
            } else {
                std::cerr << "Ignoring unknown argument: " << arg << std::endl;
            }
        } catch (const std::exception&) {
            std::cerr << "Failed to parse argument: " << arg << std::endl;
        }
    }

    const auto maxCoef = static_cast<int>(std::pow(2, beta - 1) - 1);
    const auto minCoef = static_cast<int>(-1 * std::pow(2, beta - 1));

    Solver problem(layout, maxCoef, minCoef, targets, nbInputBits, costModel);

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

    return 0;
}
