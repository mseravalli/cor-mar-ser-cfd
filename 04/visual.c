#include "helper.h"
#include "visual.h"
#include <stdio.h>
#include "string.h"

void write_vtkFile(const char *szProblem,
         		   int    timeStepNumber,
		           double xlength,
                   double ylength,
                   int    imax,
                   int    jmax,
		           double dx,
        		   double dy,
                   double **U,
                   double **V,
                   double **P) {
  
  int i,j;
  char szFileName[80];
  FILE *fp=NULL;
  sprintf( szFileName, "%s.%i.vtk", szProblem, timeStepNumber );
  fp = fopen( szFileName, "w");
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    ERROR( szBuff );
    return;
  }

  write_vtkHeader( fp, imax, jmax, dx, dy);
  write_vtkPointCoordinates(fp, imax, jmax, dx, dy);

  fprintf(fp,"POINT_DATA %i \n", (imax+1)*(jmax+1) );
	
  fprintf(fp,"\n");
  fprintf(fp, "VECTORS velocity float\n");
  for(j = 0; j < jmax+1; j++) {
    for(i = 0; i < imax+1; i++) {
      fprintf(fp, "%f %f 0\n", (U[i][j] + U[i][j+1]) * 0.5, (V[i][j] + V[i+1][j]) * 0.5 );
    }
  }

  fprintf(fp,"\n");
  fprintf(fp,"CELL_DATA %i \n", ((imax)*(jmax)) );
  fprintf(fp, "SCALARS pressure float 1 \n"); 
  fprintf(fp, "LOOKUP_TABLE default \n");
  for(j = 1; j < jmax+1; j++) {
    for(i = 1; i < imax+1; i++) {
      fprintf(fp, "%f\n", P[i][j] );
    }
  }

  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    ERROR( szBuff );
  }
}


void write_vtkHeader( FILE *fp, int imax, int jmax, 
                      double dx, double dy) {
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Null pointer in write_vtkHeader" );
    ERROR( szBuff );
    return;
  }

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"\n");	
  fprintf(fp,"DATASET STRUCTURED_GRID\n");
  fprintf(fp,"DIMENSIONS  %i %i 1 \n", imax+1, jmax+1);
  fprintf(fp,"POINTS %i float\n", (imax+1)*(jmax+1) );
  fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int imax, int jmax, 
                      double dx, double dy) {
  double originX = 0.0;  
  double originY = 0.0;

  int i = 0;
  int j = 0;

  for(j = 0; j < jmax+1; j++) {
    for(i = 0; i < imax+1; i++) {
      fprintf(fp, "%f %f 0\n", originX+(i*dx), originY+(j*dy) );
    }
  }
}

void parallelContainer(double** U, 
                       double** V, 
                       double** P, 
                       int omg_i,
                       int omg_j,
                       int imax, 
                       int jmax,
                       int iproc,
                       int jproc,
                       int nTS,             /* number of timestep */
                       char* outputFile) 
{

    char pieceExtend[256];
    char parallelCont[256];
    int i;
    int j;

    int il;
    int ir;
    int jt;
    int jb;
    
    FILE* parallelContFile = NULL;
    sprintf(parallelCont, "files/P%s.%i.vtk", outputFile, nTS);
    parallelContFile = fopen(parallelCont, "w");

    /* write header */
    fprintf(parallelContFile, "<VTKFile type=\"PStructuredGrid\"");
    fprintf(parallelContFile, 
            "<PStructuredGrid WholeExtent=\"1 imax 1 jmax 0 0\" GhostLevel=\"0\">");

    fprintf(parallelContFile, "<PPointData>");
    fprintf(parallelContFile, "</PPointData>");

    fprintf(parallelContFile, "<PCellData>");
    fprintf(parallelContFile, "</PCellData>");

    fprintf(parallelContFile, "<PPoints>");
    fprintf(parallelContFile, "</PPoints>");
    
    for (i = 0; i < iproc; ++i) {
        for (j = 0; j < jproc; ++j) {
            sprintf(pieceExtend, "%s_%i%i.%i", outputFile, i, j, nTS);

            il = (imax / iproc) * i + 1;
            jb = (jmax / jproc) * j + 1;

            if (i == iproc - 1) {
                ir = imax; 
            } else {
                ir = (i + 1) * (imax / iproc); 
            }

            if (j == jproc - 1) {
                jt = jmax; 
            } else {
                jt = (j + 1) * (jmax / jproc); 
            }

            fprintf(parallelContFile, 
                    "<Piece Extent=\"%i %i %i %i 0 0\" Source=\"%s.vts\"/>", 
                    il,
                    ir,
                    jb,
                    jt,
                    pieceExtend);

        }
    }

    fprintf(parallelContFile, "</PStructuredGrid>");
    fprintf(parallelContFile, "</VTKFile>");

}

void parallelHeader(FILE *fp,  
                    int il, 
                    int ir, 
                    int jb,
                    int jt,
                    int imax,
                    int jmax,
                    double dx, 
                    double dy)
{
    if( fp == NULL )		       
    {
      char szBuff[80];
      sprintf( szBuff, "Null pointer in write_vtkHeader" );
      ERROR( szBuff );
      return;
    }
    
    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");	
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i 1 \n", ir+1, jt+1);
    fprintf(fp,"POINTS %i float\n", (ir+1)*(jt+1) );
    fprintf(fp,"\n");

}

void parallelCoords(FILE *fp,
                    int il, 
                    int ir, 
                    int jb,
                    int jt,
                    double dx, 
                    double dy){

    double originX = il*dx;  
    double originY = jb*dx;
    
    int i = 0;
    int j = 0;
    
    for(j = jb; j < jt; j++) {
      for(i = il; i < ir; i++) {
        fprintf(fp, "%f %f 0\n", originX+(i*dx), originY+(j*dy) );
      }
    }
}


void output_vtk(double** U,
                double** V,
                double** P,
                int il,
                int ir,
                int jb,
                int jt,
                int imax, 
                int jmax,
                int omg_i,
                int omg_j,
                double dx,
                double dy,
                int timeStepNumber,
                char* outputFile
                )
{
    
  int i,j;
  char szFileName[80];
  FILE *fp=NULL;
  sprintf( szFileName, "files/%s_%i%i.%i.vtk", outputFile, omg_i, omg_j, timeStepNumber );
  fp = fopen( szFileName, "w");
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    ERROR( szBuff );
    return;
  }

  write_vtkHeader( fp, ir - il + 1, jt - jb + 1, dx, dy);
  write_vtkPointCoordinates(fp, ir - il + 1, jt - jb + 1, dx, dy);

  fprintf(fp,"POINT_DATA %i \n", (ir - il + 2)*(jt - jb + 2) );
	
  fprintf(fp,"\n");
  fprintf(fp, "VECTORS velocity float\n");
  for(j = 0; j < (jt - jb + 2); j++) {
    for(i = 0; i < (ir - il + 2); i++) {
      fprintf(fp, "%f %f 0\n", (U[i][j] + U[i][j+1]) * 0.5, (V[i][j] + V[i+1][j]) * 0.5 );
    }
  }

  fprintf(fp,"\n");
  fprintf(fp,"CELL_DATA %i \n", ((ir - il + 1)*(jt - jb + 1)) );
  fprintf(fp, "SCALARS pressure float 1 \n"); 
  fprintf(fp, "LOOKUP_TABLE default \n");
  for(j = 1; j < (jt - jb + 2); j++) {
    for(i = 1; i < (ir - il + 2); i++) {
      fprintf(fp, "%f\n", P[i][j] );
    }
  }

  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    ERROR( szBuff );
  }

}

