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

    double W1;
    double W2;
    double W3;

    READ_INT(    argv[1], *xlength);
    READ_DOUBLE( argv[1], *tau );
    
    read_double( argv[1], "W1", &W1 );
    read_double( argv[1], "W2", &W2 );
    read_double( argv[1], "W3", &W3 );
    velocityWall[0] = W1;
    velocityWall[1] = W2;
    velocityWall[2] = W3;

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

    for(z = 0; z < xlength + 2; ++z){
        for(y = 0; y < xlength + 2; ++y){
            for(x = 0; x < xlength + 2; ++x){
                
                pos = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x ); 
                /* 
                 * set NO_SLIP at the boundaries, MOVING_WALL on the top
                 * and FLUID everywhere else
                 */
                if (z == xlength + 1) {
                    flagField[pos] = MOVING_WALL; 
                } else if (x == 0 || x == xlength + 1 || y == 0 || y == xlength + 1 || z == 0) {
                    flagField[pos] = NO_SLIP; 
                } else {
                    flagField[pos] = FLUID; 
                }

                for(i = 0; i < Q; ++i){
                    collideField[Q * pos + i] = LATTICEWEIGHTS[i];
                    streamField [Q * pos + i] = LATTICEWEIGHTS[i];
                }
            }
        }
    }
    
}

