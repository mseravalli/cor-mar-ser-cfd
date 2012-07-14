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
    double  Pr;                /* prandtl number    */
    double  UI;                /* velocity x-direction */
    double  VI;                /* velocity y-direction */
    double  PI;                /* pressure */
    double  TI;
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
    int     kmax;
    int     productsNum;
    int     reactantsNum;
    double  alpha;             /* uppwind differencing factor*/
    double  beta;              /* thermal expansion coefficient */
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

    int cl;
    int cr;
    int cb;
    int ct;

    double deltaP;

    char imageName[64];
    char cavityFile[64];

    double**  U = NULL;
    double**  V = NULL;
    double**  P = NULL;
    double**  F = NULL;
    double**  G = NULL;
    double**  RS = NULL;
    double**  T = NULL;      /*  Temperature  */
    char problem[64];
    int **Problem = NULL;
    int **Flag = NULL;
    int **Sources = NULL;
    double*** C  = NULL;
    double*   C0 = NULL;
    double*** Q  = NULL;
    double**   K  = NULL;
    
    double t; 
    int n;
    int it;
    double res;
    double D;      /*Diffusion Cons.*/
    double ki;      /* Kinetic Cons.irreversible reaction */ 
    double kr;      /* Kinetic Cons. reversible reaction */
    double Ei;      /* Reaction Energy irreversible */
    double Er;      /* Reaction Energy reversible */
    double catRate;
    
    if(argn <= 1)
    {
        printf("ERROR: you need to specify a problem (karman, plane, step)\n");
        return 1;
    } else {
        if( !(   strcmp(args[1], "karman") == 0
              || strcmp(args[1], "baffle") == 0 
              || strcmp(args[1], "semibaffle") == 0 
              || strcmp(args[1], "diffusion") == 0 
              || strcmp(args[1], "plane")  == 0
              || strcmp(args[1], "step")   == 0)){
            printf("ERROR: the passed argument was different from karman, plane or step\n");
            return 1;
        }
    }

    strcpy(problem, args[1]);

    sprintf(cavityFile, "scenarios/%s_parameters.dat", problem);

    read_parameters(cavityFile,
                    &Re,
                    &Pr,     
                    &UI,     
                    &VI,     
                    &PI,
                    &TI,     
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
                    &beta,  
                    &omg,    
                    &tau,    
                    &itermax,
                    &eps,    
                    &wl,
                    &wr,
                    &wb,
                    &wt,
                    &cl,
                    &cr,
                    &cb,
                    &ct,
                    &dt_value,
                    &deltaP,
                    &D,
                    &kmax,
                    &ki,
                    &kr,
                    &Ei,
                    &Er,
                    &catRate);
                    
    t = 0;
    n = 0;

    sprintf(imageName, "scenarios/%s.pgm", problem);
    Problem = read_pgm(imageName, &imax, &jmax);
    dx = xlength / (double)(imax);
    dy = ylength / (double)(jmax);
    
    Flag = imatrix(0, imax + 1, 0, jmax + 1);
    Sources = imatrix(0, imax + 1, 0, jmax+1);

    if(init_flag(Problem, imax, jmax, Flag, Sources) == 1)
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
    T   = matrix(0, imax + 1, 0, jmax + 1);
    C   = matrix3(0, imax + 1, 0, jmax + 1, kmax);
    C0  = (double*) malloc((size_t) (kmax * sizeof(double)));
    Q   = matrix3(0, imax + 1, 0, jmax + 1, kmax);
    K = matrix(0, 3, 0, kmax);
    init_uvp(UI, 
             VI, 
             PI,
             TI, 
             imax, 
             jmax, 
             problem,
             U, 
             V, 
             P, 
             C,
             T,
             kmax);
    
    init_C0K(cavityFile, kmax, C0, K, ki, kr, &reactantsNum, &productsNum);

    while (t < t_end)
    {
        calculate_dt(Re,
                 Pr,
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
                       kmax,
                       dx,
                       dy,
                       wl,
                       wr,
                       wb,
                       wt,
                       cl,
                       cr,
                       cb,
                       ct,
                       U,
                       V,
                       F,
                       G,
                       P,
                       C,
		       T,
                       Flag,
                       Sources,
                       C0);

        spec_boundary_val(problem,
                          imax,
                          jmax,
			  kmax,
                          dx,
                          dy,
                          Re,
                          deltaP,
                          U,
                          V,
                          P,
                          C);
          calculate_t(dt,
                    dx,
                    dy,
                    alpha,
                    imax,
                    jmax,
                    Re,
                    Pr,
                    U,
                    V,
                    T,
                    Flag);



         calculate_q(K,
                    imax,
                    jmax,
                    kmax,
                    Q,
                    C,
                    Flag,
                    Ei,
                    Er,
                    T,
                    reactantsNum,
                    productsNum,
                    catRate);

          calculate_c(dt,
                    dx,
                    dy,
                    alpha,
                    D,
                    imax,
                    jmax,
                    kmax,
                    U,
                    V,
                    Q,
                    C,
                    Flag);


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
                 T,
                 beta);
        
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
            write_vtkTemperature("vtk/temperature",
                                    n,
                                    imax,
                                    jmax,
		                    dx,
		                    dy,
                                    T);
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
    write_vtkTemperature("vtk/temperature",
                            n,
                            imax,
                            jmax,
                            dx,
                            dy,
                            T);


    free_matrix(U,  0, imax + 1, 0, jmax + 1);
    free_matrix(V,  0, imax + 1, 0, jmax + 1);
    free_matrix(P,  0, imax + 1, 0, jmax + 1);
    free_matrix(F,  0, imax + 1, 0, jmax + 1);
    free_matrix(G,  0, imax + 1, 0, jmax + 1);
    free_matrix(RS, 0, imax + 1, 0, jmax + 1);
    free_matrix(T, 0, imax + 1, 0, jmax + 1);
    free_imatrix(Flag, 0, imax + 1, 0, jmax + 1);
    free_imatrix(Sources, 0, imax+1, 0, jmax+1);
    free_imatrix(Problem, 0, imax + 1, 0, jmax + 1);
    free_matrix3(C , 0, imax + 1, 0, jmax + 1, kmax);
    free(C0);
    free_matrix3(Q,  0, imax + 1, 0, jmax + 1, kmax);
    free(K);

    return 0;
}



