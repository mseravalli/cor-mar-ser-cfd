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
                   MPI Status *status,
                   int chunk)
{
    int i, j;
    int jsize, isize;
    
    jsize = jt-jb+1;
    isize = ir-il+1;

    /*** SEND LEFT, receive from right ***/
    /*copying data to be send prom pressure matrix to send buffer*/
    if(rank_l != MPI_PROC_NULL) /*if there is no left neighbour there is not need to copy data, since it will not be sent anyway*/
    {
        for(j=1; j <= jsize; j++)
        {
            buffSend[j-1] = P[1][j];
        }
    }
    /*trasfering data*/
    MPI_Sendrecv(bufSend, jsize, MPI_DOUBLE, rank_l, chunk, bufRecv, jsize, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);
    /*copying received data from receive buffer to pressure matric*/
    if(rank_r != MPI_PROC_NULL) /*if there is no right neighbour there is no need to copy data. Since this is a system boundary it will not be changed anyway*/
    {
        for(j=1; j <= jsizel j++)
        {
            P[isize+1][j] = bufRecv[j-1];
        }
    }
    
    /*** SEND RIGHT, receive from left ***/
    if(rank_r != MPI_PROC_NULL)
    {
        for(j=1; j <= jsize; j++)
        {
            buffSend[j-1] = P[isize][j];
        }
    }
    MPI_Sendrecv(bufSend, jsize, MPI_DOUBLE, rank_r, chunk, bufRecv, jsize, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);
    if(rank_l != MPI_PROC_NULL){
        for(j=1; j <= jsizel j++)
        {
            P[0][j] = bufRecv[j-1];
        }
    }
    
    /*** SEND UP, receive from bottom ***/
    if(rank_t != MPI_PROC_NULL)
    {
        for(i=1; i <= jsize; i++)
        {
            buffSend[i-1] = P[i][jsize];
        }
    }
    MPI_Sendrecv(bufSend, isize, MPI_DOUBLE, rank_t, chunk, bufRecv, isize, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    if(rank_b != MPI_PROC_NULL)
    {
        for(i=1; i <= jsizel; i++)
        {
            P[i][0] = bufRecv[i-1];
        }
    }
    
    /*** SEND DOWN, receive from up ***/
    if(rank_b != MPI_PROC_NULL)
    {
        for(i=1; i <= jsize; i++)
        {
            buffSend[i-1] = P[1][jsize];
        }
    }
    MPI_Sendrecv(bufSend, isize, MPI_DOUBLE, rank_t, chunk, bufRecv, isize, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    if(rank_t != MPI_PROC_NULL)
    {
        for(i=1; i <= jsizel; i++)
        {
            P[i][jsize+1] = bufRecv[i-1];
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
              MPI Status *status,
              int chunk)
{

}
