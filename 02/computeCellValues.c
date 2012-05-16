#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "math.h"

void computeDensity(const double *const currentCell, double *density){
    int i;
    for(i = 0; i < Q; i++)
    {
        density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    int i;
    for(i = 0; i<Q; i++)
    {
        velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1]/desity;
        velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2]/desity;
        velocity[3] += currentCell[i]*LATTICEVELOCITIES[i][3]/desity;
    }
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
    int i;
    double cs = sqrt(1/3);
    double cispu; //scalar product of ci and u
    double uspu = velocity[1]*velocity[1] + velocity[2]*velocity[2] + velocity[3]*velocity[3]; //scalar product of u and u
    for(i = 0; i < Q; i++)
    {
        cispu = LATTICEVELOCITIES[i][1]*velocity[1] + LATTICEVELOCITIES[i][2]*velocity[2] + LATTICEVELOCITIES[i][3]*velocity[3];
        feq[i] = LATTICEWEIGHTS[i] * density[i] * (1 + cispu/(cs*cs)) + (cispu*cispu)/(2*cs*cs*cs*cs) - uspu/(2*cs*cs);
    }
}

