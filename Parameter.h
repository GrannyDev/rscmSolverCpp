//
// Created by smith on 22/02/25.
//

#ifndef PARAMETER_H
#define PARAMETER_H
#include <vector>
#include <bits/stdint-intn.h>


class Parameter {
public:
    Parameter(std::vector<int64_t> const& possibleValuesCP, std::vector<int> const& possibleValuesFusion, int zeroPoint);
    ~Parameter() = default;

    std::vector<int64_t> possibleValuesCP;
    std::vector<int> possibleValuesFusion;
    int zeroPoint;
};



#endif //PARAMETER_H
