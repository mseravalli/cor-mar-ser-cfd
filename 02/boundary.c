#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){

    int x = 0;
    int y = 0;
    int z = 0;
    int i = 0;
    int pos = 0;
    int nb = 0;
    double density = 0;
    double cispu = 0;
    double wv = 0;

    for (z = 0; z < xlength + 2; ++z) {
        for (y = 0; y < xlength + 2; ++y) {
            for (x = 0; x < xlength + 2; ++x) {
            
                pos = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x ); 

                /* proceed only if we are at a boundary */
                if (flagField[pos] == NO_SLIP || flagField[pos] == MOVING_WALL) {

                    for (i = 0; i < Q; ++i) {

                        /* check if the neighbour is within the domain */
                        if (   (z+LATTICEVELOCITIES[i][2] > 0 && z+LATTICEVELOCITIES[i][2] < xlength + 1) 
                            && (y+LATTICEVELOCITIES[i][1] > 0 && y+LATTICEVELOCITIES[i][1] < xlength + 1)
                            && (x+LATTICEVELOCITIES[i][0] > 0 && x+LATTICEVELOCITIES[i][0] < xlength + 1) ) {

                            nb = ( (z+LATTICEVELOCITIES[i][2])*(xlength+2)*(xlength+2) 
                                 + (y+LATTICEVELOCITIES[i][1])*(xlength+2) 
                                 + (x+LATTICEVELOCITIES[i][0]) );

                            wv = 0;
                            if( flagField[pos] == MOVING_WALL) {
                               computeDensity(&collideField[Q*nb],&density);
                               cispu = LATTICEVELOCITIES[i][0]*wallVelocity[0] 
                                     + LATTICEVELOCITIES[i][1]*wallVelocity[1] 
                                     + LATTICEVELOCITIES[i][2]*wallVelocity[2];
                               wv = 2*LATTICEWEIGHTS[i]*density*cispu/(C_S*C_S);
                            }

                            collideField[Q * pos + i] = collideField[Q * nb + (Q-1-i)] + wv;
                        }
                    }
                }
            }
        }
    }
}


