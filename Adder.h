//
// Created by smith on 22/02/25.
//

#ifndef ADDER_H
#define ADDER_H
#include <vector>

#include "ParamDefs.h"
#include "Parameter.h"


class Adder {
public:
    Adder(int layerNb, int adderNb, int precedingLayerLastAdderNb, std::vector<ParamDefs> const& paramDefs, int maxShift);
    ~Adder() = default;

    std::vector<Parameter> parameters;
    int layerNb;
    int adderNb;
};



#endif //ADDER_H
