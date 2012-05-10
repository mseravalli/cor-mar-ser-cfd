#include "initLB.h"

int readParameters(const char *szFileName,        /* name of file*/
                   int *xlength,                 /* reads domain size. Parameter name: "xlength" */
                   double *tau,                  /* relaxation parameter tau. Parameter name: "tau" */                   
                   double *velocityWall,         /* velocity of the lid. Parameter name: "characteristicvelocity" */     
                   int *timesteps,               /* number of timesteps. Parameter name: "timesteps" */                  
                   int *timestepsPerPlotting,    /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
                   int argc,                     /* number of arguments. Should equal 2 (program + name of config file */ 
                   char *argv[])                 /* argv[1] shall contain the path to the config file */
{
 READ_INT( szFileName, *xlength);
 READ_DOUBLE( szFileName, *tau );
 READ_DOUBLE( szFileName, *velocityWall );
 READ_INT( szFileName, *timesteps );
 READ_INT( szFileName, *timestepsPerPlotting );
 READ_INT( szFileName, argc );
 READ_STRING( szFileName, *argv);



  return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
  /* TODO */
}

