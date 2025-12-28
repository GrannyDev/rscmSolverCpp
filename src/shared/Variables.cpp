//
// Created by smith on 22/02/25.
//

#include "Variables.h"

Variables::Variables(std::vector<int64_t> const& possibleValuesCP, std::vector<int> const& possibleValuesFusion, const int zeroPoint) {
    this->possibleValuesCP = possibleValuesCP;
    this->possibleValuesFusion = possibleValuesFusion;
    this->zeroPoint = zeroPoint;
}