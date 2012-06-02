#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

/* definition of the MPI messagess */
#define OMGI   10
#define OMGJ   20

#define NBL    110
#define NBR    120
#define NBT    130
#define NBB    140

#define BOUNDL 210
#define BOUNDR 220
#define BOUNDT 230
#define BOUNDB 240

void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

void init_parallel(int iproc,
                    int jproc,
                    int imax,
                    int jmax,
                    int *myrank,
                    int *il,
                    int *ir,
                    int *jb,
                    int *jt,
                    int *rank_l,
                    int *rank_r,
                    int *rank_b,
                    int *rank_t,
                    int *omg_i,
                    int *omg_j,
                    int num_proc);


/* Exchanges pressure between nighbouring processes */
void pressure_comm(double **P,
                   int il,
                   int ir,
                   int jb,
                   int jt,
                   int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
                   double *bufRecv,
                   MPI_Status *status,
                   int chunk);

/* Exchanges velocities between nighbouring processes */
void uv_comm(double **U,
              double **V,
              int il,
              int ir,
              int jb,
              int jt,
              int rank_l,
              int rank_r,
              int rank_b,
              int rank_t,
              double *bufSend,
              double *bufRecv,
              MPI_Status *status,
              int chunk);
