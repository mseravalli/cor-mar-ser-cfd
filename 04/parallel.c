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

    /* TODO: check if the number of regions is equal to the number of procs */
    /* TODO: differentiate for master?? */

    if ((*myrank) == 0) {
        for (i = 0; i < iproc; ++i) {
            for (j = 0; j < jproc; ++j) {
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
                MPI_Send(&tmpBound, 2, MPI_INT, setProc, BOUNDB, MPI_COMM_WORLD);

                /* if in the top row */
                if (j == jproc - 1) {
                    tmpNb = MPI_PROC_NULL;
                    tmpBound = jmax; 
                } else {
                    tmpNb = setProc + iproc;
                    tmpBound = (j + 1) * (jmax / jproc); 
                }
                MPI_Send(&tmpNb, 1, MPI_INT, setProc, NBT, MPI_COMM_WORLD);
                MPI_Send(&tmpBound, 2, MPI_INT, setProc, BOUNDT, MPI_COMM_WORLD);

                /* if in the left column */
                if (i == 0) {
                    tmpNb = MPI_PROC_NULL;
                } else {
                    tmpNb = setProc - jproc;
                }
                MPI_Send(&tmpNb, 1, MPI_INT, setProc, NBL, MPI_COMM_WORLD);
                tmpBound = (imax / iproc) * i + 1;
                MPI_Send(&tmpBound, 2, MPI_INT, setProc, BOUNDL, MPI_COMM_WORLD);

                /* if in the right column */
                if (i == iproc - 1) {
                    tmpNb = MPI_PROC_NULL;
                    tmpBound = imax; 
                } else {
                    tmpNb = setProc + jproc;
                    tmpBound = (i + 1) * (imax / iproc); 
                }
                MPI_Send(&tmpNb, 1, MPI_INT, setProc, NBR, MPI_COMM_WORLD);
                MPI_Send(&tmpBound, 2, MPI_INT, setProc, BOUNDR, MPI_COMM_WORLD);

                ++setProc;
            }
        }
    }

    MPI_Recv(omg_i, 1, MPI_INT, 0, OMGI, MPI_COMM_WORLD, &status);
    MPI_Recv(omg_j, 1, MPI_INT, 0, OMGJ, MPI_COMM_WORLD, &status);

    MPI_Recv(rank_b, 1, MPI_INT, 0, NBB, MPI_COMM_WORLD, &status);
    MPI_Recv(rank_t, 1, MPI_INT, 0, NBT, MPI_COMM_WORLD, &status);
    MPI_Recv(rank_l, 1, MPI_INT, 0, NBL, MPI_COMM_WORLD, &status);
    MPI_Recv(rank_r, 1, MPI_INT, 0, NBR, MPI_COMM_WORLD, &status);

    MPI_Recv(jb, 1, MPI_INT, 0, BOUNDB, MPI_COMM_WORLD, &status);
    MPI_Recv(jt, 1, MPI_INT, 0, BOUNDT, MPI_COMM_WORLD, &status);
    MPI_Recv(il, 1, MPI_INT, 0, BOUNDL, MPI_COMM_WORLD, &status);
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
            bufSend[j-1] = P[1][j];
        }
    }
    /*trasfering data*/
    MPI_Sendrecv(bufSend, jsize, MPI_DOUBLE, rank_l, chunk, bufRecv, jsize, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);
    /*copying received data from receive buffer to pressure matric*/
    if(rank_r != MPI_PROC_NULL) /*if there is no right neighbour there is no need to copy data. Since this is a system boundary it will not be changed anyway*/
    {
        for(j=1; j <= jsize; j++)
        {
            P[isize+1][j] = bufRecv[j-1];
        }
    }
    
    /*** SEND RIGHT, receive from left ***/
    if(rank_r != MPI_PROC_NULL)
    {
        for(j=1; j <= jsize; j++)
        {
            bufSend[j-1] = P[isize][j];
        }
    }
    MPI_Sendrecv(bufSend, jsize, MPI_DOUBLE, rank_r, chunk, bufRecv, jsize, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);
    if(rank_l != MPI_PROC_NULL){
        for(j=1; j <= jsize; j++)
        {
            P[0][j] = bufRecv[j-1];
        }
    }
    
    /*** SEND UP, receive from bottom ***/
    if(rank_t != MPI_PROC_NULL)
    {
        for(i=1; i <= jsize; i++)
        {
            bufSend[i-1] = P[i][jsize];
        }
    }
    MPI_Sendrecv(bufSend, isize, MPI_DOUBLE, rank_t, chunk, bufRecv, isize, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    if(rank_b != MPI_PROC_NULL)
    {
        for(i=1; i <= jsize; i++)
        {
            P[i][0] = bufRecv[i-1];
        }
    }
    
    /*** SEND DOWN, receive from up ***/
    if(rank_b != MPI_PROC_NULL)
    {
        for(i=1; i <= jsize; i++)
        {
            bufSend[i-1] = P[1][jsize];
        }
    }
    MPI_Sendrecv(bufSend, isize, MPI_DOUBLE, rank_t, chunk, bufRecv, isize, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    if(rank_t != MPI_PROC_NULL)
    {
        for(i=1; i <= jsize; i++)
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
              MPI_Status *status,
              int chunk)
{
    int imax = ir - il + 1;
    int jmax = jt - jb + 1;
    int i = 0;
    int j = 0;

    /* Velocity U*/

    /* Send-Recv Left-Right */

    if (rank_l != MPI_PROC_NULL){
        for (j= 1;j <= jmax;j++)
        {
            bufSend[j-1] = U[2][j];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, chunk, bufRecv, jmax, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);

    if (rank_r !=MPI_PROC_NULL){
        for (j= 1;j <= jmax;j++){
             U[imax+1][j]=bufRecv[j-1];
        }
    }

     /* Send-Recv Right-Left */

    if (rank_r != MPI_PROC_NULL){
        for (j= 1;j <= jmax;j++)
        {
            bufSend[j-1] = U[imax-1][j];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, chunk, bufRecv, jmax, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);

    if (rank_r !=MPI_PROC_NULL){
        for (j= 1;j <= jmax;j++){
             U[0][j]=bufRecv[j-1];
        }
    }

     /* Send-Recv Top-Bot */

    if (rank_t != MPI_PROC_NULL){
        for (i= 1; i <= imax; i++)
        {
            bufSend[i-1] = U[i][jmax];
        }
    }
       MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_t, chunk, bufRecv, jmax, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);

    if (rank_b !=MPI_PROC_NULL){
        for (i= 1;i <= imax;i++){
             U[i][0]=bufRecv[i-1];
        }
    }

    /* Send-Recv Bot-Top */

    if (rank_b != MPI_PROC_NULL){
        for (i= 1; i <= imax; i++)
        {
            bufSend[i-1] = U[i][1];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_b, chunk, bufRecv, jmax, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status);

    if (rank_t !=MPI_PROC_NULL){
        for (i= 1;i <= imax;i++){
             U[i][jmax+1]=bufRecv[i-1];
        }
    }


     /* Velocity V */

    /* Send-Recv Left-Right */

    if (rank_l != MPI_PROC_NULL){
        for (j= 0; j <= jmax + 1 ; j++)
        {
            bufSend[j] = V[1][j];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, chunk, bufRecv, jmax, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);

    if (rank_r !=MPI_PROC_NULL){
        for (j= 0; j <= jmax + 1; j++){
             V[imax+1][j]=bufRecv[j];
        }
    }
    
     /* Send-Recv Right-Left */

    if (rank_r != MPI_PROC_NULL){
        for (j= 0; j <= jmax + 1; j++)
        {
            bufSend[j] = V[imax][j];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, chunk, bufRecv, jmax, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);

    if (rank_r !=MPI_PROC_NULL){
        for (j = 0 ; j <= jmax +1; j++){
             V[0][j] = bufRecv[j];
        }
    }

     /* Send-Recv Top-Bot */

    if (rank_t != MPI_PROC_NULL){
        for (i = 1; i <= imax; i++)
        {
            bufSend[i-1] = V[i][jmax];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_t, chunk, bufRecv, jmax, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);

    if (rank_b != MPI_PROC_NULL){
        for (i= 1;i <= imax;i++){
             V[i][0] = bufRecv[i-1];
        }
    }

    /* Send-Recv Bot-Top */

    if (rank_b != MPI_PROC_NULL){
        for (i= 1; i <= imax; i++)
        {
            bufSend[i-1] = V[i][1];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_b, chunk, bufRecv, jmax, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status);

    if (rank_t !=MPI_PROC_NULL){
        for (i= 1;i <= imax;i++){
             V[i][jmax+1] = bufRecv[i-1];
        }
    }

}
