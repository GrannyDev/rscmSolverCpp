//
// Created by smith on 09/12/2025.
//

#ifndef RMCMSOLVER_VERILOG_H
#define RMCMSOLVER_VERILOG_H
#include <fstream>

#include "../solver/RSCM.h"


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
     * @param inputBW x bitwidth
     * @param targets target constants for testbench generation
     * @param scmDesigns per-target SCM designs for testbench configuration
     */
    VerilogGenerator
    (
        const RSCM& solutionNode,
        const std::string& outputUri,
        const std::vector<Layer>& layers,
        const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap,
        const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap,
        unsigned int inputBW,
        const std::vector<int>& targets,
        const std::vector<std::pair<int, std::vector<DAG>>>& scmDesigns,
        bool overwrite = false
    );

private:
    void HandleModules();
    void HandleTopModule();
    std::stringstream HandleAdder(unsigned int layerIndex, unsigned int adderIdx, unsigned int nbAdder);
    std::vector<std::pair<std::string, unsigned int>> HandleParam(
        const Variables& param,
        unsigned int paramIdx,
        unsigned int paramGlobalIdx,
        unsigned int bitPos,
        unsigned int inputBw,
        unsigned int adderIdx,
        const VariableDefs& previousParamType,
        std::stringstream& inStream
    );
    std::stringstream& HandleMultiplexer(
        std::stringstream& muxStream,
        VariableDefs paramName,
        const Variables& param,
        const std::string& inputWire,
        const std::string& outputWire,
        unsigned int bitPos,
        unsigned int nbMuxInputs,
        unsigned int adderIdx,
        unsigned int outputBW
    );

    std::stringstream HandleAdderRouting(
        unsigned int layerIndex,
        unsigned int adderIdx,
        unsigned int nbAdder,
        unsigned int& bitPos
    ) const;

    static std::string PrintWireBitWidth(unsigned int bw);
    static std::string PrintBinaryMuxCode(unsigned int nbSelectBits, unsigned int selecValue);
    static unsigned int ComputeBitWidth(int value);
    unsigned int GetParamBitWidth(VariableDefs varType, unsigned int adderIdx) const;

    std::ofstream outputFile_;
    const RSCM& solutionNode_;
    const std::vector<Layer>& layers_;
    const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap_;
    const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap_;
    const unsigned int inputBW_;
    const std::vector<int>& targets_;
    const std::vector<std::pair<int, std::vector<DAG>>>& scmDesigns_;
    const VarDefsToString varDefsTostring_;
    unsigned int muxCounter_ = 0;
};


#endif //RMCMSOLVER_VERILOG_H