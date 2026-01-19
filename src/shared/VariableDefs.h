//
// Created by smith on 22/02/25.
//

#ifndef PARAMDEFINITIONS_H
#define PARAMDEFINITIONS_H
#include <string>

/**
 * @enum VariableDefs
 * @brief Enum representing various variables definitions.
 */
enum class VariableDefs
{
    RIGHT_SHIFTS,      ///< Represents right shift operations.
    LEFT_INPUTS,       ///< Represents left input route.
    LEFT_SHIFTS,       ///< Represents left shift operations.
    RIGHT_INPUTS,      ///< Represents right input route.
    OUTPUTS_SHIFTS,    ///< Represents output shift operations.
    RIGHT_MULTIPLIER,  ///< Represents a multiplier parameter (adder/subtractor).
    LEFT_MULTIPLIER,   ///< Represents a multiplier parameter for the left input.
};

/**
 * @class VarDefsToString
 * @brief Utility class to convert VariableDefs enum values to their string representations.
 */
class VarDefsToString
{
public:
    /**
     * @brief Default constructor.
     */
    VarDefsToString() = default;

    /**
     * @brief Default destructor.
     */
    ~VarDefsToString() = default;

    /**
     * @brief Converts a VariableDefs enum value to its corresponding string representation.
     * @param variable The VariableDefs enum value to convert.
     * @return A string representation of the given VariableDefs value.
     */
    static std::string operator()(const VariableDefs variable)
    {
        switch (variable)
        {
        case VariableDefs::RIGHT_MULTIPLIER:
            return "ADD_SUB";
        case VariableDefs::LEFT_MULTIPLIER:
            return "LEFT_SIGN";
        case VariableDefs::LEFT_SHIFTS:
            return "LEFT_SHIFTS";
        case VariableDefs::RIGHT_SHIFTS:
            return "RIGHT_SHIFTS";
        case VariableDefs::OUTPUTS_SHIFTS:
            return "OUTPUTS_SHIFTS";
        case VariableDefs::LEFT_INPUTS:
            return "LEFT_INPUTS";
        case VariableDefs::RIGHT_INPUTS:
            return "RIGHT_INPUTS";
        default:
            return "";
        }
    }
};

#endif //PARAMDEFINITIONS_H
