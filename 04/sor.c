#include "sor.h"
#include <math.h>
#include <mpi.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int il,
  int ir,
  int jb,
  int jt,
  int my_rank,
  int numProc,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double *bufSend,
  double *bufRecv,
  MPI Status *status,
  int chunk
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

	/* exchanges pressure between processes */
  pressure_comm(P,
                il,
                ir,
                jb,
                jt,
                rank_l,
                rank_r,
                rank_b,
                rank_t,
                bufSend,
                bufRecv,
                status,
                chunk);

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  rloc = rloc/(imax*jmax);
  rloc = sqrt(rloc);
  /* set residual */

  /* sums all local residuals and sends result to master */
  MPI_Reduce(&rloc, res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* set boundary values */
	/* bottom row */
	if(rank_b == MPI_PROC_NULL)
	{
		for(i = 1; i <= ir-il+1; i++)
		{
			P[i][0] = P[i][1];
		}
	}
	/* top row */
	if(rank_t == MPI_PROC_NULL)
	{
		for(i = 1; i <= imax; i++)
		{
			P[i][jmax+1] = P[i][jmax];
		}
	}
	/* left column */
	if(rank_l == MPI_PROC_NULL)
	{
		for(j = 0; j <= jmax; j++)
		{
			P[0][j] = P[1][j];
		}
	}
	/* right column */
	if(rank_r == MPI_PROC_NULL)
	{
		for(j = 0; j <= jmax; j++)
		{
			P[imax+1][j] = P[imax][j];
		}
	}
}

