//
// Created by smith on 22/02/25.
//

#ifndef LAYER_H
#define LAYER_H
#include <vector>

#include "Adder.h"


class Layer {
public:
    Layer(int nbLayer, int precedingLayerLastAdderNb, int nbAdders, std::vector<ParamDefs> const& paramDefs, int maxShift);
    ~Layer() = default;

    int lastAdderNb;
    int nbLayer;
    std::vector<Adder> adders;
};



#endif //LAYER_H
