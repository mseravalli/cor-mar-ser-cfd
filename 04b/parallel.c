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
                   int imax,
                   int jmax,
                   int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
                   double *bufRecv,
                   MPI_Status *status,
                   int chunk)
{
    int i;
    int j;

    /*** SEND LEFT, receive from right ***/
    /*copying data to be send prom pressure matrix to send buffer*/
    /*if there is no left neighbour there is not need to copy data, since it will not be sent anyway*/
    if(rank_l != MPI_PROC_NULL) 
    {
        for(j = 0; j <= jmax+1; j++)
        {
            bufSend[j] = P[1][j];
        }
    }
    /*trasfering data*/
    MPI_Sendrecv(bufSend, jmax+2, MPI_DOUBLE, rank_l, chunk, bufRecv, jmax+2, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);
    /*copying received data from receive buffer to pressure matric*/
    if(rank_r != MPI_PROC_NULL) /*if there is no right neighbour there is no need to copy data. Since this is a system boundary it will not be changed anyway*/
    {
        for(j = 0; j <= jmax+1; j++)
        {
            P[imax+1][j] = bufRecv[j];
        }
    }
    
    /*** SEND RIGHT, receive from left ***/
    if(rank_r != MPI_PROC_NULL)
    {
        for(j = 0; j <= jmax+1; j++)
        {
            bufSend[j] = P[imax][j];
        }
    }
    MPI_Sendrecv(bufSend, jmax+2, MPI_DOUBLE, rank_r, chunk, bufRecv, jmax+2, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);
    if(rank_l != MPI_PROC_NULL){
        for(j = 0; j <= jmax+1; j++)
        {
            P[0][j] = bufRecv[j];
        }
    }
    
    /*** SEND UP, receive from bottom ***/
    if(rank_t != MPI_PROC_NULL)
    {
        for(i = 0; i <= imax+1; i++)
        {
            bufSend[i] = P[i][jmax];
        }
    }
    MPI_Sendrecv(bufSend, imax+2, MPI_DOUBLE, rank_t, chunk, bufRecv, imax+2, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    if(rank_b != MPI_PROC_NULL)
    {
        for(i = 0; i <= imax+1; i++)
        {
            P[i][0] = bufRecv[i];
        }
    }
    
    /*** SEND DOWN, receive from up ***/
    if(rank_b != MPI_PROC_NULL)
    {
        for(i = 0; i <= imax+1; i++)
        {
            bufSend[i] = P[i][1];
        }
    }
    MPI_Sendrecv(bufSend, imax+2, MPI_DOUBLE, rank_b, chunk, bufRecv, imax+2, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status);
    if(rank_t != MPI_PROC_NULL)
    {
        for(i = 0; i <= imax+1; i++)
        {
            P[i][jmax+1] = bufRecv[i];
        }
    }

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
