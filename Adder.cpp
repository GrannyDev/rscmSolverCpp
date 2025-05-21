//
// Created by smith on 22/02/25.
//

#include "Adder.h"

#include <cmath>

Adder::Adder(const int layerIdx, const int adderIdx, const int alpha, std::vector<VariableDefs> const& variableDefs, const int maxShift) {
    this->layerIdx = layerIdx;
    this->adderIdx = adderIdx;

    // Iterate through each parameter definition and set up possible values for the variables.
    for (const auto& paramDef : variableDefs) {
        std::vector<int64_t> possibleValuesCP; // Possible values for the CP (uses a different representation for the shift).
        std::vector<int> possibleValuesFusion; // Possible values for the Fusion Path.

        // Handle the RIGHT_MULTIPLIER parameter definition.
        if (paramDef == VariableDefs::RIGHT_MULTIPLIER) {
            possibleValuesCP.push_back(-1); // subtraction
            possibleValuesCP.push_back(1); // addition
            possibleValuesFusion.push_back(-1); // subtraction
            possibleValuesFusion.push_back(0); // addition
            possibleValuesFusion.push_back(1); // both (TODO: might not be needed ...)
            variables.emplace_back(possibleValuesCP, possibleValuesFusion, 1);
        }
        // Handle the RIGHT_SHIFTS, LEFT_SHIFTS, and OUTPUTS_SHIFTS parameter definitions.
        else if (paramDef == VariableDefs::RIGHT_SHIFTS || paramDef == VariableDefs::LEFT_SHIFTS || paramDef == VariableDefs::OUTPUTS_SHIFTS) {
            for (int i = 0; i <= maxShift; i++) {
                possibleValuesCP.push_back(std::pow(2, i)); // Powers of 2 for CP.
                possibleValuesFusion.push_back(i); // Corresponding shift values for Fusion.
            }
            variables.emplace_back(possibleValuesCP, possibleValuesFusion, 0);
        }
        // Handle all other parameter definitions (routing).
        else {
            possibleValuesCP.push_back(-1);
            possibleValuesFusion.push_back(-1);
            int zeroPoint = 1;

            // If the layer index is greater than 0, routing variables are used.
            if (layerIdx > 0) {
                for (int i = 0; i <= alpha; i++) {
                    possibleValuesCP.push_back(i);
                    possibleValuesFusion.push_back(i);
                }
                for (int i = 2; i <= maxShift + 1; i++) {
                    possibleValuesCP.push_back(-i);
                    possibleValuesFusion.push_back(-i);
                }
                zeroPoint = maxShift + 1; // Adjust the zero point based on maxShift.
            }
            variables.emplace_back(possibleValuesCP, possibleValuesFusion, zeroPoint);
        }
    }
}