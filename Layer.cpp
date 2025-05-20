//
// Created by smith on 22/02/25.
//

#include "Layer.h"

Layer::Layer(const int nbLayer, const int precedingLayerLastAdderNb, const int nbAdders, std::vector<ParamDefs> const& paramDefs, const int maxShift)
{
    this->nbLayer = nbLayer;

    for (int i = 0; i < nbAdders; i++)
    {
        adders.emplace_back(nbLayer, precedingLayerLastAdderNb + 1 + i, precedingLayerLastAdderNb, paramDefs, maxShift);
        lastAdderNb = precedingLayerLastAdderNb + 1 + i;
    }
}