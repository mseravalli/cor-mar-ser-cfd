#include "collision.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    int i;
    for(i = 0; i < Q; i++)
    {
        currentCell[i] = currentCell[i] - ((currentCell[i] - feq[i]) / (*tau)) ;
    }
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
  int x, y, z;
  int pos;
  double density, velocity[3], feq[Q];
  
  for(z = 0; z < xlength+2; z++)
  {
    for(y = 0; y < xlength+2; y++)
    {
        for(x = 0; x < xlength+2; x++)
        {
            pos = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x ); 
            if(flagField[pos]== FLUID)
            {
                computeDensity(&collideField[Q*pos], &density);
                computeVelocity(&collideField[Q*pos], &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(&collideField[Q*pos], tau, feq);
            }
        }
    }
  }
}

