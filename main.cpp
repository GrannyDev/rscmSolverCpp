//
// Created by smith on 22/02/25.
//
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "Solver.h"

int main() {
    constexpr size_t beta = 6; // maximum bit-width allowed for the constants (and intermediate values)
    constexpr size_t nbInputBits = 8; // number of bits of the input (to compute the fine-grained cost function)
    constexpr auto maxCoef = static_cast<int>(std::pow(2, beta - 1) - 1);
    constexpr auto minCoef = static_cast<int>(-1 * std::pow(2, beta - 1));

    const auto targets = {-20, -13, -8, -6, -5, -3, -2, -1, 0, 1, 2, 4, 5, 7, 12, 19}; // target const set of the RSCM

    // {1,1} describes the chosen layout, i.e. 1 adder on the first layer and 1 adder on the second layer
    Solver problem({1,1}, maxCoef, minCoef, targets, nbInputBits, true);

    problem.CPSolve(); // step 1: solve the problem with the CPSolver

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

    return 0;
}
