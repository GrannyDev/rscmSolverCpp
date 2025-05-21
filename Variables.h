//
// Created by smith on 22/02/25.
//

#ifndef PARAMETER_H
#define PARAMETER_H
#include <vector>
#include <bits/stdint-intn.h>

/**
 * @class Variables
 * @brief Represents a model variable with its associated domain.
 */
class Variables {
public:
    /**
     * @brief Constructs a Variables object.
     * @param possibleValuesCP A vector of 64-bit integers representing possible values for CP.
     * @param possibleValuesFusion A vector of integers representing possible values for Merging (bitset).
     * @param zeroPoint An integer representing the zero point value.
     */
    Variables(std::vector<int64_t> const& possibleValuesCP, std::vector<int> const& possibleValuesFusion, int zeroPoint);

    /**
     * @brief Default destructor.
     */
    ~Variables() = default;

    /**
     * @brief A vector of 64-bit integers representing possible values for CP.
     */
    std::vector<int64_t> possibleValuesCP;

    /**
     * @brief A vector of integers representing possible values for Fusion.
     */
    std::vector<int> possibleValuesFusion;

    /**
     * @brief An integer representing the zero point value.
     */
    int zeroPoint;
};

#endif //PARAMETER_H