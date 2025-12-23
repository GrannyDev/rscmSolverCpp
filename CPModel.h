//
// Created by smith on 22/02/25.
//

#ifndef CPMODEL_H
#define CPMODEL_H

#include <mutex>
#include <vector>
#include <optional>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"
#include "ortools/sat/sat_parameters.pb.h"

#include "Layer.h"
#include "DAG.h"

/**
 * @class CPModel
 * @brief Represents the constraint programming model for enumerating SCMs.
 */
class CPModel {
public:
    /**
     * @brief Constructs a CPModel object with the min and max intermediate coefficient values.
     * @param minCoef The minimum coefficient value.
     * @param maxCoef The maximum coefficient value.
     * @param minInputValue The minimum input value for the SCM.
     * @param maxInputValue The maximum input value for the SCM.
     */
    CPModel(int minCoef, int maxCoef, int minInputValue, int maxInputValue);

    /**
     * @brief Solves the model for a given target.
     *
     * @param target The target coefficient to solve for.
     * @param scmDesigns A reference to a vector of pairs, first is the target coefficient,
     * second is a vector of DAGs representing the possible designs.
     * @param pushBackMutex A mutex to ensure thread-safe updates to the scmDesigns vector.
     * @param layers A constant reference to a vector of Layer objects representing the layers.
     * @param nbBitsPerDAG The number of bits per DAG in the bitset.
     * @param varToIndexMap A constant reference to a map associating VariableDefs to their indices.
     * @param varDefs A constant reference to a vector of VariableDefs representing variable definitions.
     */
    void SolveFor(
        int target,
        std::vector<std::pair<int, std::vector<DAG>>> & scmDesigns,
        std::mutex& pushBackMutex,
        std::vector<Layer> const& layers,
        size_t nbBitsPerDAG,
        std::unordered_map<VariableDefs, unsigned int> const& varToIndexMap,
        std::vector<VariableDefs> const& varDefs,
        std::optional<unsigned int> maxSolutions = std::nullopt
    ) const;

private:
    int minCoef_; ///< The minimum coefficient value.
    int maxCoef_; ///< The maximum coefficient value.
    int minInputValue_; ///< The minimum input value for the SCM.
    int maxInputValue_; ///< The maximum input value for the SCM.

    /**
     * @brief Helper method to add output values to DAG.
     * @param scm The DAG to add values to.
     * @param multiplierAtStage The multiplier value at the current stage.
     */
    void addOutputValues(DAG& scm, const int multiplierAtStage) const {
        const int val1 = multiplierAtStage * maxInputValue_;
        const int val2 = multiplierAtStage * minInputValue_;
        scm.maxOutputValue.push_back(std::max(val1, val2));
        scm.minOutputValue.push_back(std::min(val1, val2));
        // Store coefficient's trailing zeros for accurate shift savings
        const unsigned int coeffTZ = multiplierAtStage ? __builtin_ctz(static_cast<unsigned>(std::abs(multiplierAtStage))) : std::numeric_limits<unsigned int>::max();
        scm.coefficientTrailingZeros.push_back(coeffTZ);
    }
};

#endif //CPMODEL_H
