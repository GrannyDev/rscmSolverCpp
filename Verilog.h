//
// Created by smith on 09/12/2025.
//

#ifndef RMCMSOLVER_VERILOG_H
#define RMCMSOLVER_VERILOG_H
#include <fstream>

#include "RSCM.h"


class VerilogGenerator
{
public:
    /**
     * @brief Constructor for the Verilog class.
     * @param solutionNode The RSCM solution node to convert to Verilog.
     * @param outputUri The file to write to.
     * @param overwrite If the existing file is overwritten, will throw an error if false and the file exists
     * @param layers layers in the layout
     * @param idxToVarMap maps the param idx to its VariableDefs
     * @param varToIdxMap maps the VariableDefs to its param idx
     */
    VerilogGenerator
    (
        const RSCM& solutionNode,
        const std::string& outputUri,
        const std::vector<Layer>& layers,
        const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap,
        const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap,
        bool overwrite = false
    );

private:
    void PrintAdderModules();

    std::ofstream outputFile_;
    const RSCM& solutionNode_;
    const std::vector<Layer>& layers_;
    const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap_;
    const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap_;
};


#endif //RMCMSOLVER_VERILOG_H