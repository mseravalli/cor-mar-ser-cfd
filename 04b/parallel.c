#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}


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
                   int num_proc)
{
    int i = 0;
    int j = 0;
    int tmpNb = MPI_PROC_NULL;
    int tmpBound = 0;

    int setProc = 0;

    MPI_Status status;

    if ((*myrank) == 0) {
        for (j = 0; j < jproc; ++j) {
            for (i = 0; i < iproc; ++i) {
                MPI_Send(&i, 1, MPI_INT, setProc, OMGI, MPI_COMM_WORLD);
                MPI_Send(&j, 1, MPI_INT, setProc, OMGJ, MPI_COMM_WORLD);

                /* if in the bottom row */
                if (j == 0) {
                    tmpNb = MPI_PROC_NULL;
                } else {
                    tmpNb = setProc - iproc;
                }
                MPI_Send(&tmpNb, 1, MPI_INT, setProc, NBB, MPI_COMM_WORLD);
                tmpBound = (jmax / jproc) * j + 1;
                MPI_Send(&tmpBound, 1, MPI_INT, setProc, BOUNDB, MPI_COMM_WORLD);

                /* if in the top row */
                if (j == jproc - 1) {
                    tmpNb = MPI_PROC_NULL;
                    tmpBound = jmax; 
                } else {
                    tmpNb = setProc + iproc;
                    tmpBound = (j + 1) * (jmax / jproc); 
                }
                MPI_Send(&tmpNb, 1, MPI_INT, setProc, NBT, MPI_COMM_WORLD);
                MPI_Send(&tmpBound, 1, MPI_INT, setProc, BOUNDT, MPI_COMM_WORLD);

                /* if in the left column */
                if (i == 0) {
                    tmpNb = MPI_PROC_NULL;
                } else {
                    tmpNb = setProc - 1;
                }
                MPI_Send(&tmpNb, 1, MPI_INT, setProc, NBL, MPI_COMM_WORLD);
                tmpBound = (imax / iproc) * i + 1;
                MPI_Send(&tmpBound, 1, MPI_INT, setProc, BOUNDL, MPI_COMM_WORLD);

                /* if in the right column */
                if (i == iproc - 1) {
                    tmpNb = MPI_PROC_NULL;
                    tmpBound = imax; 
                } else {
                    tmpNb = setProc + 1;
                    tmpBound = (i + 1) * (imax / iproc); 
                }
                MPI_Send(&tmpNb, 1, MPI_INT, setProc, NBR, MPI_COMM_WORLD);
                MPI_Send(&tmpBound, 1, MPI_INT, setProc, BOUNDR, MPI_COMM_WORLD);

                ++setProc;
            }
        }
    }

    MPI_Recv(omg_i, 1, MPI_INT, 0, OMGI, MPI_COMM_WORLD, &status);
    MPI_Recv(omg_j, 1, MPI_INT, 0, OMGJ, MPI_COMM_WORLD, &status);

    MPI_Recv(rank_b, 1, MPI_INT, 0, NBB, MPI_COMM_WORLD, &status);
    MPI_Recv(jb, 1, MPI_INT, 0, BOUNDB, MPI_COMM_WORLD, &status);
    MPI_Recv(rank_t, 1, MPI_INT, 0, NBT, MPI_COMM_WORLD, &status);
    MPI_Recv(jt, 1, MPI_INT, 0, BOUNDT, MPI_COMM_WORLD, &status);
    MPI_Recv(rank_l, 1, MPI_INT, 0, NBL, MPI_COMM_WORLD, &status);
    MPI_Recv(il, 1, MPI_INT, 0, BOUNDL, MPI_COMM_WORLD, &status);
    MPI_Recv(rank_r, 1, MPI_INT, 0, NBR, MPI_COMM_WORLD, &status);
    MPI_Recv(ir, 1, MPI_INT, 0, BOUNDR, MPI_COMM_WORLD, &status);

}

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
                   int chunk)
{
}

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
             int chunk)
{
}
