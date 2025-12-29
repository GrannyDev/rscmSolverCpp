//
// Persist and reload a solved instance so costs can be recomputed later.
//

#ifndef RSCMSOLVER_SNAPSHOTIO_H
#define RSCMSOLVER_SNAPSHOTIO_H

#include <string>
#include <vector>
#include <optional>

#include "../solver/RSCM.h"
#include "../shared/DAG.h"

struct SnapshotData {
    std::vector<int> layout;
    std::vector<int> targets;
    size_t nbInputBits{};
    size_t lutWidth{};
    int maxCoef{};
    int minCoef{};
    size_t nbBitsPerSCM{};
    size_t nbAdders{};
    size_t nbVariables{};
    size_t nbMuxes{};
    bool isSymmetric{};

    std::vector<unsigned int> scmIndexes;

    boost::dynamic_bitset<> rscmSet;
    std::vector<int> rscmMaxOutput;
    std::vector<int> rscmMinOutput;
    std::vector<unsigned int> rscmCoeffTZ;
    std::vector<bool> rscmIsMinus;

    std::vector<unsigned int> minShiftSavings;
    std::vector<unsigned int> variableBitWidths;
    std::vector<bool> isPlusMinus;

    std::vector<DAG> selectedScms;
};

void WriteSnapshot(const class Solver& solver, const RSCM& solutionNode, const std::string& path, bool overwrite);
std::optional<SnapshotData> ReadSnapshot(const std::string& path);

#endif // RSCMSOLVER_SNAPSHOTIO_H
