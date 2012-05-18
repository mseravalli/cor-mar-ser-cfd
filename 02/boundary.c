#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void findFluid(int* flagField, int x, int y, int z, int xlength, int* direction){

    int pos = 0;

    direction[0] = 0;
    direction[1] = 0;
    direction[2] = 0;
    
    pos = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x ); 

    if( x == 0 ){
        if(flagField[pos + 1] == FLUID){
            direction[0] = 1;
            return;
        }
    }
    if( x == xlength + 1 ){
        if(flagField[pos - 1] == FLUID){
            direction[0] = -1;
            return;
        }
    }
    if( y == 0 ){
        if(flagField[pos + (xlength + 2)] == FLUID){
            direction[1] = 1;
            return;
        }
    }
    if( y == xlength + 1 ){
        if(flagField[pos - (xlength + 2)] == FLUID){
            direction[1] = -1;
            return;
        }
    }
    if( z == 0 ){
        if(flagField[pos + ((xlength+2)*(xlength+2))] == FLUID){
            direction[2] = 1;
            return;
        }
    }
    if( z == xlength + 1 ){
        if(flagField[pos - ((xlength+2)*(xlength+2))] == FLUID){
            direction[2] = -1;
            return;
        }
    }
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){

    int x = 0;
    int y = 0;
    int z = 0;
    int i = 0;
    int fluidDirection [3];
    int dirIndex = -1;
    int pos = 0;
    int nb = 0;
    double density = 0;
    double cispu = 0;
    double wv = 0;

    for (z = 0; z < xlength + 2; ++z) {
        for (y = 0; y < xlength + 2; ++y) {
            for (x = 0; x < xlength + 2; ++x) {
            
                pos = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x ); 

                if (flagField[pos] == NO_SLIP || flagField[pos] == MOVING_WALL) {

                        dirIndex = -1;
                    
                    findFluid(flagField, x, y, z, xlength, fluidDirection);

                    for (i = 0; i < 3; ++i) {
                        if (fluidDirection[i] != 0) {
                            dirIndex = i;
                        }
                    }

                    if(dirIndex == -1)
                        continue;

                    for (i = 0; i < Q; ++i) {
                        if(LATTICEVELOCITIES[i][dirIndex] == fluidDirection[dirIndex]){
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


