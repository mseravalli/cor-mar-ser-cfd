#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q           19
#define C_S         0.5773502691896258
#define FLUID       0
#define NO_SLIP     1
#define MOVING_WALL 2

/* TODO check LATTICEVELOCITIES */

  static const int LATTICEVELOCITIES[19][3] ={ { 0, -1, -1},    /* 0 */ 
                                               {-1,  0, -1},    /* 1 */                   
                                               { 0,  0, -1},    /* 2 */                   
                                               { 1,  0, -1},    /* 3 */                   
                                               { 0,  1, -1},    /* 4 */
                                               {-1, -1,  0},    /* 5 */
                                               { 0, -1,  0},    /* 6 */
                                               { 1, -1,  0},    /* 7 */
                                               {-1,  0,  0},    /* 8 */
                                               { 0,  0,  0},    /* 9 */
                                               { 1,  0,  0},    /* 10 */
                                               {-1,  1,  0},    /* 11 */
                                               { 0,  1,  0},    /* 12 */
                                               { 1,  1,  0},    /* 13 */
                                               { 0, -1,  1},    /* 14 */
                                               {-1,  0,  1},    /* 15 */
                                               { 0,  0,  1},    /* 16 */
                                               { 1,  0,  1},    /* 17 */
                                               { 0,  1,  1} };  /* 18 */
  static const double LATTICEWEIGHTS[19] = {0.027777777777777778,    /* 0 */
                                            0.027777777777777778,    /* 1 */                   
                                            0.055555555555555556,    /* 2 */                   
                                            0.027777777777777778,    /* 3 */                   
                                            0.027777777777777778,    /* 4 */
                                            0.027777777777777778,    /* 5 */
                                            0.055555555555555556,    /* 6 */
                                            0.027777777777777778,    /* 7 */
                                            0.055555555555555556,    /* 8 */
                                            0.333333333333333333,   /* 9 */
                                            0.055555555555555556,    /* 10 */
                                            0.027777777777777778,    /* 11 */
                                            0.055555555555555556,    /* 12 */
                                            0.027777777777777778,    /* 13 */
                                            0.027777777777777778,    /* 14 */
                                            0.027777777777777778,    /* 15 */
                                            0.055555555555555556,    /* 16 */
                                            0.027777777777777778,    /* 17 */
                                            0.027777777777777778};   /* 18 */

#endif

