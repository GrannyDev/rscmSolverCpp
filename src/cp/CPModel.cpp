//
// Created by smith on 22/02/25.
//

#include "CPModel.h"

CPModel::CPModel(const int minCoef, const int maxCoef, const int minInputValue, const int maxInputValue)
{
    minCoef_ = minCoef;
    maxCoef_ = maxCoef;
    minInputValue_ = minInputValue;
    maxInputValue_ = maxInputValue;
}

void CPModel::SolveFor(
    int target,
    std::vector<std::pair<int, std::vector<DAG>>> & scmDesigns,
    std::mutex& pushBackMutex,
    std::vector<Layer> const& layers,
    size_t nbBitsPerDAG,
    std::unordered_map<VariableDefs, unsigned int> const& varToIndexMap,
    std::vector<VariableDefs> const& varDefs,
    std::optional<unsigned int> maxSolutions
) const
{
    // Initialize the constraint programming model and variable mappings
    operations_research::sat::CpModelBuilder cpModel;
    std::unordered_map<std::string, operations_research::sat::IntVar> variables;

    // Define domains for wires, shifts, and multipliers
    const operations_research::Domain wireDomain(minCoef_, maxCoef_);
    const operations_research::Domain shiftDomain = operations_research::Domain::FromValues(layers[0].adders[0].variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].possibleValuesCP);
    const operations_research::Domain multiplierDomain = operations_research::Domain::FromValues(layers[0].adders[0].variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].possibleValuesCP);

    // Initialize vectors for adders and shift constants
    std::vector<int> adders;
    std::vector<operations_research::sat::IntVar> shiftsConstants;
    for (int i = 0; i < layers[0].adders[0].variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].possibleValuesFusion.size(); i++)
    {
        adders.push_back(-1 * (i + 1));
        shiftsConstants.push_back(cpModel.NewConstant(layers[0].adders[0].variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].possibleValuesCP[i]));
    }
    std::vector<int> precedingAdders;

    // Iterate through layers and define variables and constraints
    int nbAdders = 0;
    for (size_t l = 0; l < layers.size(); l++)
    {
        precedingAdders.clear();
        int start = nbAdders;
        int end = nbAdders + static_cast<int>(layers[l].adders.size());
        for (int a = start; a < end; a++)
        {
            nbAdders++;
            precedingAdders.push_back(a);

            // Define variables for the current adder
            operations_research::sat::IntVar y = cpModel.NewIntVar(wireDomain).WithName("y" + std::to_string(a));
            variables["y" + std::to_string(a)] = y; // adder output
            operations_research::sat::IntVar mr = cpModel.NewIntVar(multiplierDomain).WithName("mr" + std::to_string(a));
            variables["mr" + std::to_string(a)] = mr; // adder/subtractor
            operations_research::sat::IntVar wl = cpModel.NewIntVar(shiftDomain).WithName("wl" + std::to_string(a));
            variables["wl" + std::to_string(a)] = wl; // left shift
            operations_research::sat::IntVar wr = cpModel.NewIntVar(shiftDomain).WithName("wr" + std::to_string(a));
            variables["wr" + std::to_string(a)] = wr; // right shift
            operations_research::sat::IntVar wo = cpModel.NewIntVar(shiftDomain).WithName("wo" + std::to_string(a));
            variables["wo" + std::to_string(a)] = wo; // output shift

            // Define constraints for the first layer
            if (l == 0)
            {
                // links the output of the adder to its inputs with respect to the shifts and the wire domains
                operations_research::sat::IntVar wrTIMESmr = cpModel.NewIntVar(wireDomain).WithName("wrTIMESmr" + std::to_string(a));
                cpModel.AddMultiplicationEquality(wrTIMESmr, {wr, mr});
                operations_research::sat::IntVar wrTIMESmrPLUSwl = cpModel.NewIntVar(wireDomain).WithName("wrTIMESmrPLUSwl" + std::to_string(a));
                cpModel.AddEquality(wrTIMESmrPLUSwl, operations_research::sat::LinearExpr::Sum({wl, wrTIMESmr}));
                operations_research::sat::IntVar wrTIMESmrPLUSwlTIMESwo = cpModel.NewIntVar(wireDomain).WithName("wrTIMESmrPLUSwlTIMESwo" + std::to_string(a));
                variables["ybeforewo" + std::to_string(a)] = wrTIMESmrPLUSwl; // needed to compute max output value in the cost function for adder output
                cpModel.AddMultiplicationEquality(wrTIMESmrPLUSwlTIMESwo, {wrTIMESmrPLUSwl, wo});
                cpModel.AddEquality(y, wrTIMESmrPLUSwlTIMESwo);
            }
            else
            {
                // Define constraints for subsequent layers (routing needed)
                operations_research::Domain adderDomain(-1 * static_cast<int64_t>(layers[0].adders[0].variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].possibleValuesCP.size()),
                    adders.back()); // Negative values for the possible shifts of X, positive values if it inputs the output of previous adders
                operations_research::sat::IntVar lia = cpModel.NewIntVar(adderDomain).WithName("lia" + std::to_string(a)); // the left input routing
                operations_research::sat::BoolVar isLiaNeg = cpModel.NewBoolVar().WithName("isLiaNeg" + std::to_string(a));
                cpModel.AddGreaterOrEqual(lia, 0).OnlyEnforceIf(isLiaNeg.Not());
                cpModel.AddLessThan(lia, 0).OnlyEnforceIf(isLiaNeg);
                operations_research::sat::IntVar li = cpModel.NewIntVar(wireDomain).WithName("li" + std::to_string(a)); // the left input value
                variables["lia" + std::to_string(a)] = lia;
                variables["li" + std::to_string(a)] = li;
                operations_research::sat::IntVar ria = cpModel.NewIntVar(adderDomain).WithName("ria" + std::to_string(a));
                operations_research::sat::BoolVar isRiaNeg = cpModel.NewBoolVar().WithName("isRiaNeg" + std::to_string(a));
                cpModel.AddGreaterOrEqual(ria, 0).OnlyEnforceIf(isRiaNeg.Not());
                cpModel.AddLessThan(ria, 0).OnlyEnforceIf(isRiaNeg);
                operations_research::sat::IntVar ri = cpModel.NewIntVar(wireDomain).WithName("ri" + std::to_string(a));
                variables["ria" + std::to_string(a)] = ria;
                variables["ri" + std::to_string(a)] = ri;
                cpModel.AddBoolOr({isLiaNeg.Not(), isRiaNeg.Not()}); // if both are negative, then this adder would be useless
                operations_research::sat::IntVar riTIMESwr = cpModel.NewIntVar(wireDomain).WithName("riTIMESwr" + std::to_string(a));
                cpModel.AddMultiplicationEquality(riTIMESwr, {ri, wr});
                variables["riTIMESwr" + std::to_string(a)] = riTIMESwr;
                operations_research::sat::IntVar liTIMESwl = cpModel.NewIntVar(wireDomain).WithName("liTIMESwl" + std::to_string(a));
                cpModel.AddMultiplicationEquality(liTIMESwl, {li, wl});
                variables["liTIMESwl" + std::to_string(a)] = liTIMESwl;
                operations_research::sat::IntVar riTIMESwrTIMESmr = cpModel.NewIntVar(wireDomain).WithName("liTIMESwlTIMESml" + std::to_string(a));
                cpModel.AddMultiplicationEquality(riTIMESwrTIMESmr, {mr, riTIMESwr});
                operations_research::sat::IntVar riTIMESwrTIMESmrPLUSliTIMESwl = cpModel.NewIntVar(wireDomain).WithName("riTIMESwrTIMESmrPLUSliTIMESwl" + std::to_string(a));
                cpModel.AddEquality(riTIMESwrTIMESmrPLUSliTIMESwl, operations_research::sat::LinearExpr::Sum({riTIMESwrTIMESmr, liTIMESwl}));
                operations_research::sat::IntVar riTIMESwrTIMESmrPLUSliTIMESwlTIMESwo = cpModel.NewIntVar(wireDomain).WithName("riTIMESwrTIMESmrPLUSliTIMESwlTIMESwo" + std::to_string(a));
                variables["ybeforewo" + std::to_string(a)] = riTIMESwrTIMESmrPLUSliTIMESwl;
                cpModel.AddMultiplicationEquality(riTIMESwrTIMESmrPLUSliTIMESwlTIMESwo, {riTIMESwrTIMESmrPLUSliTIMESwl, wo});
                cpModel.AddEquality(y, riTIMESwrTIMESmrPLUSliTIMESwlTIMESwo);
                for (int & adder : adders) // for each adder in the current layer
                {
                    operations_research::sat::BoolVar liaIsAdder = cpModel.NewBoolVar().WithName("liaIsAdder" + std::to_string(adder));
                    operations_research::sat::BoolVar riaIsAdder = cpModel.NewBoolVar().WithName("riaIsAdder" + std::to_string(adder));
                    cpModel.AddEquality(lia, adder).OnlyEnforceIf(liaIsAdder);
                    cpModel.AddNotEqual(lia, adder).OnlyEnforceIf(liaIsAdder.Not());
                    cpModel.AddEquality(ria, adder).OnlyEnforceIf(riaIsAdder);
                    cpModel.AddNotEqual(ria, adder).OnlyEnforceIf(riaIsAdder.Not());
                    if (adder <= -1) // if the input is X shifted
                    {
                        cpModel.AddEquality(li, shiftsConstants[abs(adder) - 1]).OnlyEnforceIf(liaIsAdder);
                        cpModel.AddEquality(ri, shiftsConstants[abs(adder) - 1]).OnlyEnforceIf(riaIsAdder);
                    } else // if the input is the output of a previous adder
                    {
                        cpModel.AddEquality(li, variables["y" + std::to_string(adder)]).OnlyEnforceIf(liaIsAdder);
                        cpModel.AddEquality(ri, variables["y" + std::to_string(adder)]).OnlyEnforceIf(riaIsAdder);
                    }
                }
            }
        }
        adders.insert(adders.end(), precedingAdders.begin(), precedingAdders.end());
    }
    cpModel.AddEquality(variables["y" + std::to_string(adders.back())], target);

    // Creating the DAG from each CP solution by gathering the assignment of each variable for each adder in each layer
    operations_research::sat::Model model;
    std::vector<DAG> DAGVector;
    size_t solutionCount = 0;
    model.Add(operations_research::sat::NewFeasibleSolutionObserver([&](const operations_research::sat::CpSolverResponse& r) {
        DAG scm(nbBitsPerDAG, nbAdders, nbAdders * layers[0].adders[0].variables.size());

        size_t currentBit = 0;
        for(auto const& l : layers)
        {
            for (auto const& a : l.adders)
            {
                scm.isMinus.push_back(static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])) == -1);
                for (auto const& p : varDefs)
                {
                    if (p == VariableDefs::RIGHT_MULTIPLIER)
                    {
                        // std::cout << "mr " << currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].zeroPoint + static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])) << " " << static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])) << std::endl;
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].zeroPoint + static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])));
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].possibleValuesFusion.size();
                        addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["ybeforewo" + std::to_string(a.adderIdx)])));
                    } else if (p == VariableDefs::LEFT_SHIFTS)
                    {
                        // std::cout << "wl " << currentBit + a.variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].zeroPoint + static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wl" + std::to_string(a.adderIdx)]))) << " " << static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wl" + std::to_string(a.adderIdx)]))) << std::endl;
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].zeroPoint + static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wl" + std::to_string(a.adderIdx)]))));
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].possibleValuesFusion.size();
                        if (l.layerIdx > 0)
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["liTIMESwl" + std::to_string(a.adderIdx)])));
                        } else
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["wl" + std::to_string(a.adderIdx)])));
                        }
                    } else if (p == VariableDefs::RIGHT_SHIFTS)
                    {
                        // std::cout << "wr " << currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_SHIFTS)].zeroPoint + static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wr" + std::to_string(a.adderIdx)]))) << " " << static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wr" + std::to_string(a.adderIdx)]))) << std::endl;
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_SHIFTS)].zeroPoint + static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wr" + std::to_string(a.adderIdx)]))));
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::RIGHT_SHIFTS)].possibleValuesFusion.size();
                        if (l.layerIdx > 0)
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["riTIMESwr" + std::to_string(a.adderIdx)])));
                        } else
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["wr" + std::to_string(a.adderIdx)])));
                        }
                    } else if (p == VariableDefs::RIGHT_INPUTS)
                    {
                        if (l.layerIdx > 0)
                        {
                            // std::cout << "ria " << a.variables[varToIndexMap.at(VariableDefs::RIGHT_INPUTS)].zeroPoint + SolutionIntegerValue(r, variables["ria" + std::to_string(a.adderIdx)]) << " " << SolutionIntegerValue(r, variables["ria" + std::to_string(a.adderIdx)]) << std::endl;
                            scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_INPUTS)].zeroPoint + SolutionIntegerValue(r, variables["ria" + std::to_string(a.adderIdx)]));
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["ri" + std::to_string(a.adderIdx)])));
                        } else
                        {
                            // std::cout << "ri " << currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_INPUTS)].zeroPoint - 1 << " " << -1 << std::endl;
                            scm.set.set(currentBit);
                            addOutputValues(scm, 1);
                        }
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::RIGHT_INPUTS)].possibleValuesFusion.size();
                    } else if (p == VariableDefs::LEFT_INPUTS)
                    {
                        if (l.layerIdx > 0)
                        {
                            // std::cout << "lia " << a.variables[varToIndexMap.at(VariableDefs::LEFT_INPUTS)].zeroPoint + SolutionIntegerValue(r, variables["lia" + std::to_string(a.adderIdx)]) << " " << SolutionIntegerValue(r, variables["lia" + std::to_string(a.adderIdx)]) << std::endl;
                            scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::LEFT_INPUTS)].zeroPoint + SolutionIntegerValue(r, variables["lia" + std::to_string(a.adderIdx)]));
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["li" + std::to_string(a.adderIdx)])));
                        } else
                        {
                            // std::cout << "li " << currentBit + a.variables[varToIndexMap.at(VariableDefs::LEFT_INPUTS)].zeroPoint - 1 << " " << -1 << std::endl;
                            scm.set.set(currentBit);
                            addOutputValues(scm, 1);
                        }
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::LEFT_INPUTS)].possibleValuesFusion.size();
                    } else if (p == VariableDefs::OUTPUTS_SHIFTS)
                    {
                        // std::cout << "wo " << currentBit + a.variables[varToIndexMap.at(VariableDefs::OUTPUTS_SHIFTS)].zeroPoint + static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wo" + std::to_string(a.adderIdx)]))) << " " << static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wo" + std::to_string(a.adderIdx)]))) << std::endl;
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::OUTPUTS_SHIFTS)].zeroPoint + static_cast<int>(std::log2(SolutionIntegerValue(r, variables["wo" + std::to_string(a.adderIdx)]))));
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::OUTPUTS_SHIFTS)].possibleValuesFusion.size();
                        addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["y" + std::to_string(a.adderIdx)])));
                    }
                }
            }
        }
        DAGVector.push_back(scm);
    }));

    operations_research::sat::SatParameters parameters;
    parameters.set_enumerate_all_solutions(true);
    model.Add(NewSatParameters(parameters));

    if (maxSolutions.has_value())
    {
        int num_solutions = 0;
        model.Add(operations_research::sat::NewFeasibleSolutionObserver([&](const operations_research::sat::CpSolverResponse&) {
          num_solutions++;
          if (num_solutions >= *maxSolutions) {
            StopSearch(&model);
          }
        }));
    }

    SolveCpModel(cpModel.Build(), &model);
    if (!DAGVector.empty())
    {
        std::lock_guard lock(pushBackMutex);
        scmDesigns.emplace_back(target, DAGVector);
    }
}
