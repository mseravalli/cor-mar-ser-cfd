#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){
    
    double* collideField = NULL;
    double* streamField  = NULL;
    int*    flagField    = NULL;
    int     xlength;
    double  tau;
    double  velocityWall[3];
    int     timesteps;
    int     timestepsPerPlotting;
    int     readResult = 1;

    int t = 0;
    double* swap = NULL;
    
    /* read the parameters */
    readResult = readParameters(&xlength,             
                                &tau,              
                                velocityWall,     
                                &timesteps,           
                                &timestepsPerPlotting,
                                argc,                 
                                argv);

    if(readResult != 0){
        printf("could nor read some parameters, program will exit\n");
        return 1;
    }

    /* initialise the pointers */
    collideField = calloc(Q*(xlength+2)*(xlength+2)*(xlength+2), sizeof(double));
    streamField  = calloc(Q*(xlength+2)*(xlength+2)*(xlength+2), sizeof(double));
    streamField  = calloc((xlength+2)*(xlength+2)*(xlength+2), sizeof(int));
    
    for(t = 0; t < timesteps; ++t){
        doStreaming(collideField, streamField, flagField, xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;

        doCollision(collideField, flagField, &tau, xlength);
        treatBoundary(collideField, flagField, velocityWall, xlength);

        if(t % timestepsPerPlotting == 0){
            writeVtkOutput(collideField, flagField, argv[2], t, xlength);
        }
    }

    free(collideField);
    free(streamField); 
    free(flagField); 

    return 0;
}

#endif

