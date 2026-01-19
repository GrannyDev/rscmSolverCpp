//
// Created by smith on 22/02/25.
//

#include "CPModel.h"

#include <algorithm>
#include <cstdlib>

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
    const auto& shiftParam = layers[0].adders[0].variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)];
    std::vector<int64_t> shiftValues(shiftParam.possibleValuesCP.begin(), shiftParam.possibleValuesCP.end());
    std::vector<int64_t> shiftExps;
    shiftExps.reserve(shiftParam.possibleValuesFusion.size());
    for (const int v : shiftParam.possibleValuesFusion) {
        shiftExps.push_back(v);
    }
    const operations_research::Domain shiftValueDomain = operations_research::Domain::FromValues(shiftValues);
    const operations_research::Domain shiftExpDomain = operations_research::Domain::FromValues(shiftExps);
    const operations_research::Domain multiplierDomain =
        operations_research::Domain::FromValues(layers[0].adders[0].variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].possibleValuesCP);

    // Pre-build shift constants for routing from X.
    std::vector<operations_research::sat::IntVar> shiftConstants;
    shiftConstants.reserve(shiftValues.size());
    for (const auto value : shiftValues) {
        shiftConstants.push_back(cpModel.NewConstant(value));
    }

    // Helpers for routing.
    auto buildRoutingValues = [](const Variables& param) {
        std::vector<int> xVals;
        std::vector<int> adderVals;
        for (const int v : param.possibleValuesFusion) {
            if (v < 0) xVals.push_back(v);
            else adderVals.push_back(v);
        }
        std::sort(xVals.begin(), xVals.end(), std::greater<int>());
        std::sort(adderVals.begin(), adderVals.end());
        std::vector<int> values;
        values.reserve(param.possibleValuesFusion.size());
        values.insert(values.end(), xVals.begin(), xVals.end());
        values.insert(values.end(), adderVals.begin(), adderVals.end());
        return values;
    };
    // Iterate through layers and define variables and constraints
    int nbAdders = 0;
    for (size_t l = 0; l < layers.size(); l++)
    {
        int start = nbAdders;
        int end = nbAdders + static_cast<int>(layers[l].adders.size());
        for (int a = start; a < end; a++)
        {
            nbAdders++;

            // Define variables for the current adder
            operations_research::sat::IntVar y = cpModel.NewIntVar(wireDomain).WithName("y" + std::to_string(a));
            variables["y" + std::to_string(a)] = y; // adder output
            operations_research::sat::IntVar ml = cpModel.NewIntVar(multiplierDomain).WithName("ml" + std::to_string(a));
            variables["ml" + std::to_string(a)] = ml; // left input sign
            operations_research::sat::IntVar mr = cpModel.NewIntVar(multiplierDomain).WithName("mr" + std::to_string(a));
            variables["mr" + std::to_string(a)] = mr; // right input sign
            operations_research::sat::IntVar wlExp = cpModel.NewIntVar(shiftExpDomain).WithName("wl" + std::to_string(a));
            variables["wl" + std::to_string(a)] = wlExp; // left shift exponent
            operations_research::sat::IntVar wlVal = cpModel.NewIntVar(shiftValueDomain).WithName("wlv" + std::to_string(a));
            cpModel.AddElement(wlExp, shiftValues, wlVal);
            variables["wlv" + std::to_string(a)] = wlVal;
            operations_research::sat::IntVar wrExp = cpModel.NewIntVar(shiftExpDomain).WithName("wr" + std::to_string(a));
            variables["wr" + std::to_string(a)] = wrExp; // right shift exponent
            operations_research::sat::IntVar wrVal = cpModel.NewIntVar(shiftValueDomain).WithName("wrv" + std::to_string(a));
            cpModel.AddElement(wrExp, shiftValues, wrVal);
            variables["wrv" + std::to_string(a)] = wrVal;
            operations_research::sat::IntVar woExp = cpModel.NewIntVar(shiftExpDomain).WithName("wo" + std::to_string(a));
            variables["wo" + std::to_string(a)] = woExp; // output shift exponent
            operations_research::sat::IntVar woVal = cpModel.NewIntVar(shiftValueDomain).WithName("wov" + std::to_string(a));
            cpModel.AddElement(woExp, shiftValues, woVal);
            variables["wov" + std::to_string(a)] = woVal;

            // Define constraints for the first layer
            if (l == 0)
            {
                // links the output of the adder to its inputs with respect to the shifts and the wire domains
                operations_research::sat::IntVar wrTIMESmr = cpModel.NewIntVar(wireDomain).WithName("wrTIMESmr" + std::to_string(a));
                cpModel.AddMultiplicationEquality(wrTIMESmr, {wrVal, mr});
                operations_research::sat::IntVar wlTIMESml = cpModel.NewIntVar(wireDomain).WithName("wlTIMESml" + std::to_string(a));
                cpModel.AddMultiplicationEquality(wlTIMESml, {wlVal, ml});
                variables["wlTIMESml" + std::to_string(a)] = wlTIMESml;
                operations_research::sat::IntVar wrTIMESmrPLUSwlTIMESml = cpModel.NewIntVar(wireDomain).WithName("wrTIMESmrPLUSwlTIMESml" + std::to_string(a));
                cpModel.AddEquality(wrTIMESmrPLUSwlTIMESml, operations_research::sat::LinearExpr::Sum({wlTIMESml, wrTIMESmr}));
                operations_research::sat::IntVar wrTIMESmrPLUSwlTIMESmlTIMESwo = cpModel.NewIntVar(wireDomain).WithName("wrTIMESmrPLUSwlTIMESmlTIMESwo" + std::to_string(a));
                variables["ybeforewo" + std::to_string(a)] = wrTIMESmrPLUSwlTIMESml; // needed to compute max output value in the cost function for adder output
                cpModel.AddMultiplicationEquality(wrTIMESmrPLUSwlTIMESmlTIMESwo, {wrTIMESmrPLUSwlTIMESml, woVal});
                cpModel.AddEquality(y, wrTIMESmrPLUSwlTIMESmlTIMESwo);
            }
            else
            {
                // Define constraints for subsequent layers (routing needed)
                const auto& adder = layers[l].adders[static_cast<size_t>(a - start)];
                const auto& leftInputParam = adder.variables[varToIndexMap.at(VariableDefs::LEFT_INPUTS)];
                const auto& rightInputParam = adder.variables[varToIndexMap.at(VariableDefs::RIGHT_INPUTS)];

                const auto leftRoutingValues = buildRoutingValues(leftInputParam);
                const auto rightRoutingValues = buildRoutingValues(rightInputParam);

                operations_research::sat::IntVar liaIdx =
                    cpModel.NewIntVar(operations_research::Domain(0, static_cast<int64_t>(leftRoutingValues.size() - 1)))
                        .WithName("lia" + std::to_string(a));
                operations_research::sat::IntVar riaIdx =
                    cpModel.NewIntVar(operations_research::Domain(0, static_cast<int64_t>(rightRoutingValues.size() - 1)))
                        .WithName("ria" + std::to_string(a));
                variables["lia" + std::to_string(a)] = liaIdx;
                variables["ria" + std::to_string(a)] = riaIdx;

                operations_research::sat::IntVar li = cpModel.NewIntVar(wireDomain).WithName("li" + std::to_string(a));
                operations_research::sat::IntVar ri = cpModel.NewIntVar(wireDomain).WithName("ri" + std::to_string(a));
                variables["li" + std::to_string(a)] = li;
                variables["ri" + std::to_string(a)] = ri;

                std::vector<operations_research::sat::LinearExpr> leftSources;
                leftSources.reserve(leftRoutingValues.size());
                for (const int v : leftRoutingValues) {
                    if (v < 0) {
                        leftSources.emplace_back(shiftConstants[std::abs(v) - 1]);
                    } else {
                        leftSources.emplace_back(variables["y" + std::to_string(v)]);
                    }
                }
                std::vector<operations_research::sat::LinearExpr> rightSources;
                rightSources.reserve(rightRoutingValues.size());
                for (const int v : rightRoutingValues) {
                    if (v < 0) {
                        rightSources.emplace_back(shiftConstants[std::abs(v) - 1]);
                    } else {
                        rightSources.emplace_back(variables["y" + std::to_string(v)]);
                    }
                }

                cpModel.AddElement(liaIdx, leftSources, li);
                cpModel.AddElement(riaIdx, rightSources, ri);

                const auto leftXCount = static_cast<int>(std::count_if(
                    leftRoutingValues.begin(),
                    leftRoutingValues.end(),
                    [](const int v) { return v < 0; }
                ));
                const auto rightXCount = static_cast<int>(std::count_if(
                    rightRoutingValues.begin(),
                    rightRoutingValues.end(),
                    [](const int v) { return v < 0; }
                ));
                const operations_research::sat::BoolVar liaIsX = cpModel.NewBoolVar();
                const operations_research::sat::BoolVar riaIsX = cpModel.NewBoolVar();
                if (leftXCount == 0) {
                    cpModel.AddEquality(liaIsX, 0);
                } else if (leftXCount == static_cast<int>(leftRoutingValues.size())) {
                    cpModel.AddEquality(liaIsX, 1);
                } else {
                    cpModel.AddLessOrEqual(liaIdx, leftXCount - 1).OnlyEnforceIf(liaIsX);
                    cpModel.AddGreaterOrEqual(liaIdx, leftXCount).OnlyEnforceIf(liaIsX.Not());
                }
                if (rightXCount == 0) {
                    cpModel.AddEquality(riaIsX, 0);
                } else if (rightXCount == static_cast<int>(rightRoutingValues.size())) {
                    cpModel.AddEquality(riaIsX, 1);
                } else {
                    cpModel.AddLessOrEqual(riaIdx, rightXCount - 1).OnlyEnforceIf(riaIsX);
                    cpModel.AddGreaterOrEqual(riaIdx, rightXCount).OnlyEnforceIf(riaIsX.Not());
                }
                cpModel.AddBoolOr({liaIsX.Not(), riaIsX.Not()}); // if both are X, then this adder would be useless

                operations_research::sat::IntVar riTIMESwr = cpModel.NewIntVar(wireDomain).WithName("riTIMESwr" + std::to_string(a));
                cpModel.AddMultiplicationEquality(riTIMESwr, {ri, wrVal});
                variables["riTIMESwr" + std::to_string(a)] = riTIMESwr;
                operations_research::sat::IntVar liTIMESwl = cpModel.NewIntVar(wireDomain).WithName("liTIMESwl" + std::to_string(a));
                cpModel.AddMultiplicationEquality(liTIMESwl, {li, wlVal});
                variables["liTIMESwl" + std::to_string(a)] = liTIMESwl;
                operations_research::sat::IntVar riTIMESwrTIMESmr = cpModel.NewIntVar(wireDomain).WithName("riTIMESwrTIMESmr" + std::to_string(a));
                cpModel.AddMultiplicationEquality(riTIMESwrTIMESmr, {mr, riTIMESwr});
                operations_research::sat::IntVar liTIMESwlTIMESml = cpModel.NewIntVar(wireDomain).WithName("liTIMESwlTIMESml" + std::to_string(a));
                cpModel.AddMultiplicationEquality(liTIMESwlTIMESml, {ml, liTIMESwl});
                variables["liTIMESwlTIMESml" + std::to_string(a)] = liTIMESwlTIMESml;
                operations_research::sat::IntVar riTIMESwrTIMESmrPLUSliTIMESwlTIMESml = cpModel.NewIntVar(wireDomain).WithName("riTIMESwrTIMESmrPLUSliTIMESwlTIMESml" + std::to_string(a));
                cpModel.AddEquality(riTIMESwrTIMESmrPLUSliTIMESwlTIMESml, operations_research::sat::LinearExpr::Sum({riTIMESwrTIMESmr, liTIMESwlTIMESml}));
                operations_research::sat::IntVar riTIMESwrTIMESmrPLUSliTIMESwlTIMESmlTIMESwo = cpModel.NewIntVar(wireDomain).WithName("riTIMESwrTIMESmrPLUSliTIMESwlTIMESmlTIMESwo" + std::to_string(a));
                variables["ybeforewo" + std::to_string(a)] = riTIMESwrTIMESmrPLUSliTIMESwlTIMESml;
                cpModel.AddMultiplicationEquality(riTIMESwrTIMESmrPLUSliTIMESwlTIMESmlTIMESwo, {riTIMESwrTIMESmrPLUSliTIMESwlTIMESml, woVal});
                cpModel.AddEquality(y, riTIMESwrTIMESmrPLUSliTIMESwlTIMESmlTIMESwo);
            }
        }
    }
    cpModel.AddEquality(variables["y" + std::to_string(nbAdders - 1)], target);

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
                scm.isLeftMinus.push_back(static_cast<int>(SolutionIntegerValue(r, variables["ml" + std::to_string(a.adderIdx)])) == -1);
                scm.isRightMinus.push_back(static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])) == -1);
                for (auto const& p : varDefs)
                {
                    if (p == VariableDefs::RIGHT_MULTIPLIER)
                    {
                        // std::cout << "mr " << currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].zeroPoint + static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])) << " " << static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])) << std::endl;
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].zeroPoint + static_cast<int>(SolutionIntegerValue(r, variables["mr" + std::to_string(a.adderIdx)])));
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::RIGHT_MULTIPLIER)].possibleValuesFusion.size();
                        addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["ybeforewo" + std::to_string(a.adderIdx)])));
                    } else if (p == VariableDefs::LEFT_MULTIPLIER)
                    {
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::LEFT_MULTIPLIER)].zeroPoint + static_cast<int>(SolutionIntegerValue(r, variables["ml" + std::to_string(a.adderIdx)])));
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::LEFT_MULTIPLIER)].possibleValuesFusion.size();
                        if (l.layerIdx > 0)
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["liTIMESwlTIMESml" + std::to_string(a.adderIdx)])));
                        } else
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["wlTIMESml" + std::to_string(a.adderIdx)])));
                        }
                    } else if (p == VariableDefs::LEFT_SHIFTS)
                    {
                        const int wlExp = static_cast<int>(SolutionIntegerValue(r, variables["wl" + std::to_string(a.adderIdx)]));
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].zeroPoint + wlExp);
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::LEFT_SHIFTS)].possibleValuesFusion.size();
                        if (l.layerIdx > 0)
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["liTIMESwl" + std::to_string(a.adderIdx)])));
                        } else
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["wlv" + std::to_string(a.adderIdx)])));
                        }
                    } else if (p == VariableDefs::RIGHT_SHIFTS)
                    {
                        const int wrExp = static_cast<int>(SolutionIntegerValue(r, variables["wr" + std::to_string(a.adderIdx)]));
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_SHIFTS)].zeroPoint + wrExp);
                        currentBit += a.variables[varToIndexMap.at(VariableDefs::RIGHT_SHIFTS)].possibleValuesFusion.size();
                        if (l.layerIdx > 0)
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["riTIMESwr" + std::to_string(a.adderIdx)])));
                        } else
                        {
                            addOutputValues(scm, static_cast<int>(SolutionIntegerValue(r, variables["wrv" + std::to_string(a.adderIdx)])));
                        }
                    } else if (p == VariableDefs::RIGHT_INPUTS)
                    {
                        if (l.layerIdx > 0)
                        {
                            const auto rightRoutingValues =
                                buildRoutingValues(a.variables[varToIndexMap.at(VariableDefs::RIGHT_INPUTS)]);
                            const int riaIdx = static_cast<int>(SolutionIntegerValue(r, variables["ria" + std::to_string(a.adderIdx)]));
                            scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::RIGHT_INPUTS)].zeroPoint + rightRoutingValues[riaIdx]);
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
                            const auto leftRoutingValues =
                                buildRoutingValues(a.variables[varToIndexMap.at(VariableDefs::LEFT_INPUTS)]);
                            const int liaIdx = static_cast<int>(SolutionIntegerValue(r, variables["lia" + std::to_string(a.adderIdx)]));
                            scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::LEFT_INPUTS)].zeroPoint + leftRoutingValues[liaIdx]);
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
                        const int woExp = static_cast<int>(SolutionIntegerValue(r, variables["wo" + std::to_string(a.adderIdx)]));
                        scm.set.set(currentBit + a.variables[varToIndexMap.at(VariableDefs::OUTPUTS_SHIFTS)].zeroPoint + woExp);
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
