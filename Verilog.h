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
    void PrintModules();
    void PrintTestbench();
    std::string GenerateMuxCode(
        const Variables& param,
        size_t paramBitPos,
        unsigned int paramGlobalIdx,
        unsigned int nbPossibleValues,
        const std::string& muxName,
        unsigned int bitwidth,
        unsigned int adderIdx,
        unsigned int inputBitWidth,
        size_t nbVarsPerAdder
    ) const;
    
    // Helper methods
    static unsigned int ComputeBitWidth(int maxVal) ;
    unsigned int GetParamBitWidth(VariableDefs varType, unsigned int adderIdx, size_t nbVarsPerAdder) const;
    void ProcessParameter(
        const Variables& param,
        size_t& bitPos,
        unsigned int& paramGlobalIdx,
        unsigned int paramInAdderIdx,
        unsigned int adderIdx,
        size_t nbVarsPerAdder,
        std::unordered_map<unsigned int, unsigned int>& muxInputsMap,
        std::stringstream& left_shift_mux,
        std::stringstream& right_shift_mux,
        std::stringstream& outputShiftMux,
        std::stringstream& plusMinus,
        std::stringstream& left_input,
        std::stringstream& right_input,
        unsigned int leftInputPortBW,
        unsigned int rightInputPortBW,
        unsigned int adderOutputBW
    ) const;

    static std::string PrintWireBitWidth(unsigned int bw);
    static std::string PrintBinaryMuxCode(unsigned int nbSelectBits, unsigned int selecValue) ;

    std::ofstream outputFile_;
    const RSCM& solutionNode_;
    const std::vector<Layer>& layers_;
    const std::unordered_map<unsigned int, VariableDefs>& idxToVarMap_;
    const std::unordered_map<VariableDefs, unsigned int>& varToIdxMap_;
    const unsigned int inputBW_;
    const std::vector<int>& targets_;
    const std::vector<std::pair<int, std::vector<DAG>>>& scmDesigns_;
};


#endif //RMCMSOLVER_VERILOG_H