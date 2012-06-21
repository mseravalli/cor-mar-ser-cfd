#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include <string.h>


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

    double  Re;                /* reynolds number   */          
    double  UI;                /* velocity x-direction */
    double  VI;                /* velocity y-direction */
    double  PI;                /* pressure */
    double  C0;              /* initial values for substance concentration */
    double  C1;
    double  C2;
    double  C3;
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
    int     kmax = 4;
    double  alpha;             /* uppwind differencing factor*/
    double  omg;               /* relaxation factor */
    double  tau;               /* safety factor for time step*/
    int     itermax;             /* max. number of iterations  */
                                /* for pressure per time step */
    double  eps;               /* accuracy bound for pressure*/
    double  dt_value;           /* time for output */              

    int wl;
    int wr;
    int wt;
    int wb;

    double deltaP;

    char imageName[64];
    char cavityFile[64];

    double**  U = NULL;
    double**  V = NULL;
    double**  P = NULL;
    double*** C = NULL;
    double**  F = NULL;
    double**  G = NULL;
    double**  RS = NULL;
    char problem[64];
    int **Problem = NULL;
    int **Flag = NULL;
    
    double t; 
    int n;
    int it;
    double res;
    double D = 5.3;

    if(argn <= 1)
    {
        printf("ERROR: you need to specify a problem (karman, plane, step)\n");
        return 1;
    } else {
        if( !(   strcmp(args[1], "karman") == 0 
              || strcmp(args[1], "plane")  == 0
              || strcmp(args[1], "step")   == 0)){
            printf("ERROR: the passed argument was different from karman, plane or step\n");
            return 1;
        }
    }

    strcpy(problem, args[1]);

    strcpy(cavityFile, problem);
    strcat(cavityFile, "_cavity.dat");

    read_parameters(cavityFile,
                    &Re,     
                    &UI,     
                    &VI,     
                    &PI,     
                    &C0,
                    &C1,
                    &C2,
                    &C3,
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
                    &wl,
                    &wr,
                    &wt,
                    &wb,
                    &dt_value,
                    &deltaP);
                    
    t = 0;
    n = 0;

    strcpy(imageName, problem);
    strcat(imageName, ".pgm");
    Problem = read_pgm(imageName, &imax, &jmax);
    dx = xlength / (double)(imax);
    dy = ylength / (double)(jmax);
    
    Flag = imatrix(0, imax + 1, 0, jmax + 1);

    if(init_flag(Problem, imax, jmax, Flag) == 1)
    {
        /* 
         * if there was a forbidden obstacle it returns an error, 
         * frees everything and finishes the program
         */
        printf("ERROR: Invalid obstacle in .pgm file\n");

        free_imatrix(Problem, 0, imax + 1, 0, jmax + 1);
        free_imatrix(Flag, 0, imax + 1, 0, jmax + 1);

        return 1;
    }

    U   = matrix(0, imax + 1, 0, jmax + 1); 
    V   = matrix(0, imax + 1, 0, jmax + 1); 
    P   = matrix(0, imax + 1, 0, jmax + 1);
    F   = matrix(0, imax + 1, 0, jmax + 1);
    G   = matrix(0, imax + 1, 0, jmax + 1);
    RS  = matrix(0, imax + 1, 0, jmax + 1);
    C   = matrix3(0, imax + 1, 0, jmax + 1, kmax);
    init_uvp(UI, 
             VI, 
             PI, 
             C0,
             C1,
             C2,
             C3,
             imax, 
             jmax, 
             problem,
             U, 
             V, 
             P, 
             C);

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
                 V,
                 Flag,
                 D);
                 
        boundaryvalues(imax,
                       jmax,
                       wl,
                       wr,
                       wt,
                       wb,
                       U,
                       V,
                       F,
                       G,
                       P,
                       Flag);

        spec_boundary_val(problem,
                          imax,
                          jmax,
                          dx,
                          dy,
                          Re,
                          deltaP,
                          U,
                          V,
                          P);

        calculate_fg(Re,
                 GX,
                 GY,
                 alpha,
                 dt,
                 dx,
                 dy,
                 imax,
                 jmax,
                 U,
                 V,
                 F,
                 G,
                 Flag,
                 wl,
                 wr,
                 wt,
                 wb);
        
        calculate_rs(dt,
                 dx,
                 dy,
                 imax,
                 jmax,
                 F,
                 G,
                 RS,
                 Flag);
                 
        it = 0;
        res = eps + 1;
        
        while (it < itermax && res > eps)
        {

             sor(
                 omg,
                 dx,
                 dy,
                 imax,
                 jmax,
                 deltaP,
                 P,
                 RS,
                 Flag,
                 &res,
                 problem);
                 
            it++;
        }
        
        calculate_uv(dt,
                 dx,
                 dy,
                 imax,
                 jmax,
                 U,
                 V,
                 F,
                 G,
                 P,
                 Flag);

        if( ((int)t) % ((int)dt_value) == 0 
            && t > n*dt_value){
            write_vtkFile("vtk/cavity",
		                  n,
		                  xlength,
                          ylength,
                          imax,
                          jmax,
                          kmax,
                		  dx,
		                  dy,
                          U,
                          V,
                          P,
                          C);
            write_vtkConcentrations("vtk/concentration",
                                    n,
                                    imax,
                                    jmax,
                                    kmax,
		                            dx,
		                            dy,
                                    C);
            n++;
        }

        t += dt;
    }

    write_vtkFile("vtk/cavity",
		          n,
		          xlength,
                  ylength,
                  imax,
                  jmax,
                  kmax,
                  dx,
		          dy,
                  U,
                  V,
                  P,
                  C);
    write_vtkConcentrations("vtk/concentration",
                            n,
                            imax,
                            jmax,
                            kmax,
                            dx,
                            dy,
                            C);


    free_matrix(U,  0, imax + 1, 0, jmax + 1);
    free_matrix(V,  0, imax + 1, 0, jmax + 1);
    free_matrix(P,  0, imax + 1, 0, jmax + 1);
    free_matrix(F,  0, imax + 1, 0, jmax + 1);
    free_matrix(G,  0, imax + 1, 0, jmax + 1);
    free_matrix(RS, 0, imax + 1, 0, jmax + 1);
    free_imatrix(Flag, 0, imax + 1, 0, jmax + 1);
    free_imatrix(Problem, 0, imax + 1, 0, jmax + 1);
    free_matrix3(C , 0, imax + 1, 0, jmax + 1, kmax);

    return 0;
}




