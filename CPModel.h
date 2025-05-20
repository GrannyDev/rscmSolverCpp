//
// Created by smith on 22/02/25.
//

#ifndef CPMODEL_H
#define CPMODEL_H

#include <mutex>
#include <vector>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"
#include "ortools/sat/sat_parameters.pb.h"

#include "Layer.h"
#include "Rcm.h"

class CPModel {
public:
    CPModel(int minCoef, int maxCoef);
    void SolveFor(
        int target,
        std::vector<std::pair<int, std::vector<Rcm>>> & scmDesigns,
        std::mutex& pushBackMutex,
        std::vector<Layer> const& layers,
        size_t nbBitsPerNode,
        std::unordered_map<ParamDefs, unsigned int> const& paramToIndexMap,
        std::vector<ParamDefs> const& paramDefs
    ) const;

private:

    int minCoef_;
    int maxCoef_;
};

#endif //CPMODEL_H
