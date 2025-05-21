//
// Created by smith on 22/02/25.
//

#ifndef LAYER_H
#define LAYER_H
#include <vector>

#include "Adder.h"

/**
 * @class Layer
 * @brief Represents a layer consisting of one or multiple adders.
 */
class Layer {
public:
    /**
     * @brief Constructs a Layer object.
     *
     * @param layerIdx The index of the current layer.
     * @param alpha The index of the last adder in the preceding layer.
     * @param nbAdders The number of adders in the current layer.
     * @param variableDefs A vector of variables definitions.
     * @param maxShift The maximum shift value allowed in the shift variables.
     */
    Layer(int layerIdx, int alpha, int nbAdders, std::vector<VariableDefs> const& variableDefs, int maxShift);

    /**
     * @brief Default destructor for the Layer class.
     */
    ~Layer() = default;

    int alpha; ///< The index of the last adder of the preceding layer
    int layerIdx; ///< The index of the current layer.
    std::vector<Adder> adders; ///< A collection of adders in the current layer.
};

#endif //LAYER_H