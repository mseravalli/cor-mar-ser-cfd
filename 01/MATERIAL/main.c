#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>


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
                    &dt_value);

    init_uvp(UI, VI, PI, imax, jmax, U, V, P);
    
    free_matrix(U, 0, imax + 1, 0, jmax + 1);

/*    free_matrix(V, 0, imax + 1, 0, jmax + 1);
    free_matrix(P, 0, imax, 0, jmax);
*/
    return 1;
}




