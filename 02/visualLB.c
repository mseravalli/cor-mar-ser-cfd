#include "visualLB.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void write_vtkHeader( FILE *fp, int xlength) {
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Null pointer in write_vtkHeader" );
    ERROR( szBuff );
    return;
  }

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"generated by CFD-lab course output (written by Marco Seravalli) \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"\n");	
  fprintf(fp,"DATASET STRUCTURED_GRID\n");
  fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength, xlength, xlength);
  fprintf(fp,"POINTS %i float\n", (xlength)*(xlength)*(xlength) );
  fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int xlength ) {
    double originX = 0.0;  
    double originY = 0.0;
    double originZ = 0.0;
    
    int x = 0;
    int y = 0;
    int z = 0;

    double delta = 1/(double)(xlength - 1);

    for(z = 0; z < xlength; ++z) {
        for(y = 0; y < xlength; ++y) {
            for(x = 0; x < xlength; ++x) {
                fprintf(fp, "%f %f %f\n", originX+(x*delta), originY+(y*delta), originZ+(z*delta) );
            }
        }
    }
}


void writeVtkOutput(const double * const collideField, 
                    const int * const flagField, 
                    const char *filename, 
                    int t,
                    int xlength)
{
/*int i,j;*/
    char szFileName[80];
    FILE *fp=NULL;

    int x = 0;
    int y = 0;
    int z = 0;
    int pos = 0;
    double density = 0;
    double vel[3];

    /* try to open the file, if it does not exists abort */
    sprintf( szFileName, "%s.%i.vtk", filename, t );
    fp = fopen( szFileName, "w");
    if( fp == NULL ){
      char szBuff[80];
      sprintf( szBuff, "Failed to open %s", szFileName );
      ERROR( szBuff );
      return;
    }

    write_vtkHeader(fp, xlength);
    write_vtkPointCoordinates(fp, xlength);

    fprintf(fp,"\n");
    fprintf(fp,"POINT_DATA %i \n", ((xlength) * (xlength) * (xlength)) );
    fprintf(fp, "VECTORS velocity float \n"); 
    for(z = 1; z <= xlength; ++z) {
        for(y = 1; y <= xlength; ++y) {
            for(x = 1; x <= xlength; ++x) {
                pos = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x ); 
                    computeDensity(&collideField[Q*pos], &density);
                    computeVelocity(&collideField[Q*pos], &density, vel);
                    fprintf(fp, "%f %f %f\n", vel[0], vel[1], vel[2] );
            }
        }
    }
  
}


