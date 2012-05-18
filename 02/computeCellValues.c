#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "math.h"

void computeDensity(const double *const currentCell, double *density){
    int i;
    *density = 0;
    for(i = 0; i < Q; i++)
    {
        *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    int i;
    velocity[0] = 0;
    velocity[1] = 0;
    velocity[2] = 0;

    for (i = 0; i < Q; i++) {
        velocity[0] += currentCell[i]*LATTICEVELOCITIES[i][0]/(*density);
        velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1]/(*density);
        velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2]/(*density);
    }
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
    int i;
    double cispu; /* scalar product of ci and u */
    double uspu = velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]; /* scalar product of u and u */
    for(i = 0; i < Q; i++)
    {
        cispu = LATTICEVELOCITIES[i][0]*velocity[0] + LATTICEVELOCITIES[i][1]*velocity[1] + LATTICEVELOCITIES[i][2]*velocity[2];
        feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + cispu/(C_S*C_S)) + (cispu*cispu)/(2*C_S*C_S*C_S*C_S) - uspu/(2*C_S*C_S);
    }
}

