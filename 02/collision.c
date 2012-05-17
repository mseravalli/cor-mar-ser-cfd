#include "collision.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    int i;
    for(i = 0; i < Q; i++)
    {
        currentCell[i] = currentCell[i] - 1/ (*tau) * (currentCell[i] - feq[i]);
    }
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
  int x, y, z;
  double density, velocity[3], feq[Q];
  
  for(z = 0; z < xlength+2; z++)
  {
    for(y = 0; y < xlength+2; y++)
    {
        for(x = 0; x < xlength+2; x++)
        {
            if(flagField[z*xlength*xlength + y*xlength + x] == FLUID)
            {
                computeDensity(&collideField[Q*(z*xlength*xlength + y*xlength + x)], &density);
                computeVelocity(&collideField[Q*(z*xlength*xlength + y*xlength + x)], &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(&collideField[Q*(z*xlength*xlength + y*xlength + x)], tau, feq);
            }
        }
    }
  }
}

