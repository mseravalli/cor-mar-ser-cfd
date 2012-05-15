#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength,                 /* reads domain size. Parameter name: "xlength" */
                   double *tau,                  /* relaxation parameter tau. Parameter name: "tau" */                   
                   double *velocityWall,         /* velocity of the lid. Parameter name: "characteristicvelocity" */     
                   int *timesteps,               /* number of timesteps. Parameter name: "timesteps" */                  
                   int *timestepsPerPlotting,    /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
                   int argc,                     /* number of arguments. Should equal 2 (program + name of config file */ 
                   char *argv[])                 /* argv[1] shall contain the path to the config file */
{
    if(argc < 2){
        printf("filename not passed \n");    
        return 1;
    }

    READ_INT(    argv[1], *xlength);
    READ_DOUBLE( argv[1], *tau );
    READ_DOUBLE( argv[1], *velocityWall );
    READ_INT(    argv[1], *timesteps );
    READ_INT(    argv[1], *timestepsPerPlotting );

    return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){

    int x = 0;
    int y = 0;
    int z = 0;  
    int i = 0;
    int pos = 0;

    /* from 0 to 20 */
    for(z = 0; z < xlength + 2; ++z){
        /* from 0 to 20 */
        for(y = 0; y < xlength + 2; ++y){
            /* from 0 to 20 */
            for(x = 0; x < xlength + 2; ++x){
                
                pos = ( z*xlength*xlength + y*xlength + x ); 
                /* 
                 * set NO_SLIP at the boundaries, MOVING_WALL on the top
                 * and FLUID everywhere else
                 */
                if(x == 0 || x == Q - 1 || y == 0 || y == Q - 1 || z == 0){
                    flagField[pos] = NO_SLIP; 
                } else if(z == Q - 1){
                    flagField[pos] = MOVING_WALL; 
                } else {
                    flagField[pos] = FLUID; 
                }

                for(i = 0; i < Q; ++i){
                    /* process only the fluid cells */
                    /* TODO: the boundary values ?? */
                    if(flagField[pos] == FLUID){
                        collideField[Q * pos + i] = LATTICEWEIGHTS[i];
                        streamField [Q * pos + i] = LATTICEWEIGHTS[i];
                    } else {                        
                        collideField[Q * pos + i] = LATTICEWEIGHTS[i];
                        streamField [Q * pos + i] = LATTICEWEIGHTS[i];
                    }
                }
            }
        }
    }
    
}

