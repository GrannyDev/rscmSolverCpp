//
// Created by smith on 22/02/25.
//

#include "Layer.h"

Layer::Layer(const int layerIdx, const int alpha, const int nbAdders, std::vector<VariableDefs> const& variableDefs, const int maxShift)
{
    this->layerIdx = layerIdx; // Set the index of the current layer.

    // Loop through the number of adders and initialize each one.
    for (int i = 0; i < nbAdders; i++)
    {
        // Create an Adder object and add it to the adders vector.
        adders.emplace_back(layerIdx, alpha + 1 + i, alpha, variableDefs, maxShift);

        // Update the alpha value for the current layer.
        this->alpha = alpha + 1 + i;
    }
}