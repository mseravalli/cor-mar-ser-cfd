#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include "mpi.h"
#include "parallel.h"


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){

    const char *szFileName = "cavity100.dat";       /* name of the file */
    double  Re;                /* reynolds number   */          
    double  UI;                /* velocity x-direction */
    double  VI;                /* velocity y-direction */
    double  PI;                /* pressure */
    double  GX;                /* gravitation x-direction */
    double  GY;                /* gravitation y-direction */
    double  t_end;             /* end time */
    double  xlength;           /* length of the domain x-dir.*/
    double  ylength;           /* length of the domain y-dir.*/
    double  dt;                /* time step */
    double  dx;                /* length of a cell x-dir. */
    double  dy;                /* length of a cell y-dir. */
    int     imax;                /* number of cells x-direction*/
    int     jmax;                /* number of cells y-direction*/
    double  alpha;             /* uppwind differencing factor*/
    double  omg;               /* relaxation factor */
    double  tau;               /* safety factor for time step*/
    int     itermax;             /* max. number of iterations  */
                                /* for pressure per time step */
    double  eps;               /* accuracy bound for pressure*/
    double  dt_value;           /* time for output */              

    double **U = NULL;
    double **V = NULL;
    double **P = NULL;
    double **F = NULL;
    double **G = NULL;
    double **RS = NULL;
    
    double t; 
    int n;
    int it;
    double res;
    char nextSORiter;
    
    int iproc;
    int jproc;
    int omg_i;
    int omg_j;

    int il;
    int ir;
    int jb;
    int jt;
    int max_length;

    int rank_l;
    int rank_r;
    int rank_b;
    int rank_t;
    
    int num_proc;
    int myrank;

    double *bufSend;
    double *bufRecv;

    /**** initialise all the variables ****/
    MPI_Init(&argn, &args);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    read_parameters(szFileName,
                    &Re,     
                    &UI,     
                    &VI,     
                    &PI,     
                    &GX,     
                    &GY,     
                    &t_end,  
                    &xlength,
                    &ylength,
                    &dt,     
                    &dx,     
                    &dy,     
                    &imax,   
                    &jmax,   
                    &alpha,  
                    &omg,    
                    &tau,    
                    &itermax,
                    &eps,    
                    &dt_value,
                    &iproc,
                    &jproc);
                    
    /* if the number of processes is not coherent with the subdomains, exit */
    if (num_proc % (iproc*jproc) != 0) {
        if (myrank == 0) {
            printf("cannot execute : number of processes and subdomains should be coherent\n");
        }
        MPI_Finalize();
        return 1;
    }

    init_parallel(iproc,
                  jproc,
                  imax,
                  jmax,
                  &myrank,
                  &il,
                  &ir,
                  &jb,
                  &jt,
                  &rank_l,
                  &rank_r,
                  &rank_b,
                  &rank_t,
                  &omg_i,
                  &omg_j,
                  num_proc);

    imax = ir - il + 1;
    jmax = jt - jb + 1;

    t = 0;
    n = 0;

    /**** allocate all the matrices and the vectores ****/
    U = matrix(0, imax + 2, 0, jmax + 1); 
    F = matrix(0, imax + 2, 0, jmax + 1);
    V = matrix(0, imax + 1, 0, jmax + 2); 
    G = matrix(0, imax + 1, 0, jmax + 2);
    P = matrix(0, imax + 1, 0, jmax + 1);
    RS = matrix(0, imax + 1, 0, jmax + 1); 

    max_length = imax+3;
    if(max_length < jmax+3)
    {
      max_length = jmax+3;    
    }
    
    bufSend = (double *) malloc((size_t)(max_length * sizeof(double)));
    bufRecv = (double *) malloc((size_t)(max_length * sizeof(double)));

    /**** initialise the matrices ****/
    init_matrix(U, 0, imax + 2, 0, jmax + 1, UI); 
    init_matrix(V, 0, imax + 1, 0, jmax + 2, VI); 
    init_matrix(P, 0, imax + 1, 0, jmax + 1, PI);

    while (t < t_end)
    {
        calculate_dt(Re,
                 tau,
                 &dt,
                 dx,
                 dy,
                 imax,
                 jmax,
                 U,
                 V);

        boundaryvalues(imax,
                       jmax,
                       rank_l,
                       rank_r,
                       rank_b,
                       rank_t,
                       U,
                       V);

        calculate_fg(Re,
                 GX,
                 GY,
                 alpha,
                 dt,
                 dx,
                 dy,
                 imax,
                 jmax,
                 rank_l,
                 rank_r,
                 rank_b,
                 rank_t,
                 U,
                 V,
                 F,
                 G);

        calculate_rs(dt,
                 dx,
                 dy,
                 imax,
                 jmax,
                 F,
                 G,
                 RS);
                 
        it = 0;
        res = eps + 1;
        nextSORiter = 1;

        while (nextSORiter)
        {
             sor(
                 omg,
                 dx,
                 dy,
                 imax,
                 jmax,
                 rank_l,
                 rank_r,
                 rank_b,
                 rank_t,
                 it,
                 bufSend,
                 bufRecv,
                 P,
                 RS,
                 &res);

            /*master decides if another iteration in required*/
            if(myrank == 0){
              nextSORiter = ((it < itermax && res > eps)? 1 : 0);
            }
            
            /*it broadcasts the decision to all other processes*/
            MPI_Bcast(&nextSORiter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                 
            it++;
        }

        break;
        
        calculate_uv(dt,
                 dx,
                 dy,
                 imax,
                 jmax,
                 U,
                 V,
                 F,
                 G,
                 P);

            write_vtkFile("files/file",
		                  n,
		                  xlength,
                          ylength,
                          imax,
                          jmax,
                		  dx,
		                  dy,
                          U,
                          V,
                          P);

        t += dt;
        n++;
    }

    write_vtkFile("files/file",
		                  n,
		                  xlength,
                          ylength,
                          imax,
                          jmax,
                		  dx,
		                  dy,
                          U,
                          V,
                          P);


    free(bufSend);
    free(bufRecv);

    free_matrix(U, 0, imax + 2, 0, jmax + 1);
    free_matrix(F, 0, imax + 2, 0, jmax + 1);
    free_matrix(V, 0, imax + 1, 0, jmax + 2);
    free_matrix(G, 0, imax + 1, 0, jmax + 2);
    free_matrix(P, 0, imax + 1, 0, jmax + 1);
    free_matrix(RS, 0, imax + 1, 0, jmax + 1);

    Programm_Stop("end correctly reached");

    MPI_Finalize();

    return 0;
}




