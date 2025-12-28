//
// Created by smith on 28/12/2025.
//

#ifndef RSCMSOLVER_COSTMODEL_H
#define RSCMSOLVER_COSTMODEL_H

/**
 * @enum CostModel
 * @brief Enumerates available cost computation strategies.
 */
enum class CostModel {
    AreaCost,
    MuxCount,
    MuxBits,
    LutsCost,
    FPGADelay,
    ASICDelay
};

#endif //RSCMSOLVER_COSTMODEL_H