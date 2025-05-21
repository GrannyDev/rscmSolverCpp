//
// Created by smith on 22/02/25.
//

#ifndef ADDER_H
#define ADDER_H
#include <vector>

#include "VariableDefs.h"
#include "Variables.h"

/**
 * @class Adder
 * @brief Represents an adder.
 */
class Adder {
public:
    /**
     * @brief Constructs an Adder object with the specified parameters.
     *
     * @param layerIdx The index of the layer to which this adder belongs.
     * @param adderIdx The index of this adder within its layer.
     * @param alpha The index of the last adder in the preceding layer.
     * @param variableDefs A constant reference to a vector of variables definitions.
     * @param maxShift The maximum shift value allowed for this adder.
     */
    Adder(int layerIdx, int adderIdx, int alpha, std::vector<VariableDefs> const& variableDefs, int maxShift);

    /**
     * @brief Default destructor for the Adder class.
     */
    ~Adder() = default;

    /**
     * @brief A vector storing the variables associated with this adder.
     */
    std::vector<Variables> variables;

    /**
     * @brief The index of the layer to which this adder belongs.
     */
    int layerIdx;

    /**
     * @brief The index of this adder within its layer.
     */
    int adderIdx;
};

#endif //ADDER_H