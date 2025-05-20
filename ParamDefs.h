//
// Created by smith on 22/02/25.
//

#ifndef PARAMDEFINITIONS_H
#define PARAMDEFINITIONS_H
#include <string>

enum class ParamDefs
{
    RIGHT_SHIFTS,
    LEFT_INPUTS,
    LEFT_SHIFTS,
    RIGHT_INPUTS,
    OUTPUTS_SHIFTS,
    RIGHT_MULTIPLIER,
};

class ParamDefsToString
{
public:
    ParamDefsToString() = default;
    ~ParamDefsToString() = default;

    static std::string operator()(const ParamDefs parameter)
    {
        switch (parameter)
        {
        case ParamDefs::RIGHT_MULTIPLIER:
            return "MULTIPLIER";
        case ParamDefs::LEFT_SHIFTS:
            return "LEFT SHIFTS";
        case ParamDefs::RIGHT_SHIFTS:
            return "RIGHT SHIFTS";
        case ParamDefs::OUTPUTS_SHIFTS:
            return "OUTPUTS SHIFTS";
        case ParamDefs::LEFT_INPUTS:
            return "LEFT INPUTS";
        case ParamDefs::RIGHT_INPUTS:
            return "RIGHT INPUTS";
        default:
            return "";
        }
    }
};



#endif //PARAMDEFINITIONS_H
