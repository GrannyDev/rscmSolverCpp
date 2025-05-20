//
// Created by smith on 22/02/25.
//

#include "Adder.h"

#include <cmath>

Adder::Adder(const int layerNb, const int adderNb, const int precedingLayerLastAdderNb, std::vector<ParamDefs> const& paramDefs, const int maxShift) {
    this->layerNb = layerNb;
    this->adderNb = adderNb;

    for (const auto& paramDef : paramDefs) {
        std::vector<int64_t> possibleValuesCP;
        std::vector<int> possibleValuesFusion;
        if (paramDef == ParamDefs::RIGHT_MULTIPLIER)
        {
            possibleValuesCP.push_back(-1);
            possibleValuesCP.push_back(1);
            possibleValuesFusion.push_back(-1);
            possibleValuesFusion.push_back(0);
            possibleValuesFusion.push_back(1);
            parameters.emplace_back(possibleValuesCP, possibleValuesFusion, 1);
        } else if (paramDef == ParamDefs::RIGHT_SHIFTS || paramDef == ParamDefs::LEFT_SHIFTS || paramDef == ParamDefs::OUTPUTS_SHIFTS)
        {
            for (int i = 0; i <= maxShift; i++)
            {
                possibleValuesCP.push_back(std::pow(2, i));
                possibleValuesFusion.push_back(i);
            }
            parameters.emplace_back(possibleValuesCP, possibleValuesFusion, 0);
        } else
        {
            possibleValuesCP.push_back(-1);
            possibleValuesFusion.push_back(-1);
            int zeroPoint = 1;
            if (layerNb > 0)
            {
                for (int i = 0; i <= precedingLayerLastAdderNb; i++)
                {
                    possibleValuesCP.push_back(i);
                    possibleValuesFusion.push_back(i);
                }
                for (int i = 2; i <= maxShift + 1; i++)
                {
                    possibleValuesCP.push_back(-i);
                    possibleValuesFusion.push_back(-i);
                }
                zeroPoint = maxShift + 1;
            }
            parameters.emplace_back(possibleValuesCP, possibleValuesFusion, zeroPoint);
        }
    }
}
