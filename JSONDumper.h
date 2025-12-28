//
// Created by smith on 16/12/2025.
//

#ifndef RMCMSOLVER_JSONDUMPER_H
#define RMCMSOLVER_JSONDUMPER_H
#include <fstream>
#include <unordered_map>
#include <optional>

#include "RSCM.h"

class JSONDumper {
    public:
    /**
     * @brief Constructor for the Verilog class.
     * @param solutionNode The RSCM solution node to convert to Verilog.
     * @param outputUri The file to write to.
     * @param overwrite If the existing file is overwritten, will throw an error if false and the file exists
     * @param layers layers in the layout
     * @param idxToVarMap maps the param idx to its VariableDefs
     * @param varToIdxMap maps the VariableDefs to its param idx
     * @param inputBW x bitwidth
     * @param targets target constants for testbench generation
     * @param scmDesigns per-target SCM designs for testbench configuration
     * @param costs per-objective function computed cost
     */
    JSONDumper
    (
        const RSCM& solutionNode,
        const std::string& outputUri,
        const std::vector<Layer>& layers,
        const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap,
        const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap,
        unsigned int inputBW,
        const std::vector<int>& targets,
        const std::vector<std::pair<int, std::vector<DAG>>>& scmDesigns,
        const std::unordered_map<std::string, std::optional<unsigned int>>& costs,
        bool overwrite = false
    );

    void Dump();

private:
    std::ofstream outputFile_;
    const RSCM& solutionNode_;
    const std::vector<Layer>& layers_;
    const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap_;
    const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap_;
    const unsigned int inputBW_;
    const std::vector<int>& targets_;
    const std::vector<std::pair<int, std::vector<DAG>>>& scmDesigns_;
    const std::unordered_map<std::string, std::optional<unsigned int>>& costs_;
    const VarDefsToString varDefsTostring_;

    // Helper methods
    unsigned int GetParamBitWidth(const VariableDefs varType, const unsigned int adderIdx) const;
    unsigned int GetTotalAdderCount() const;
};


#endif //RMCMSOLVER_JSONDUMPER_H
