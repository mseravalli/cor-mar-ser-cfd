#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

/* TODO check LATTICEVELOCITIES */

  static const int Q = 19;
  static const int FLUID = 0;
  static const int NO_SLIP = 1;
  static const int MOVING_WALL = 2;
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
  static const double LATTICEWEIGHTS[19] = {1/36,    /* 0 */
                                            1/36,    /* 1 */                   
                                            2/36,    /* 2 */                   
                                            1/36,    /* 3 */                   
                                            1/36,    /* 4 */
                                            1/36,    /* 5 */
                                            2/36,    /* 6 */
                                            1/36,    /* 7 */
                                            2/36,    /* 8 */
                                            0,       /* 9 */
                                            2/36,    /* 10 */
                                            1/36,    /* 11 */
                                            2/36,    /* 12 */
                                            1/36,    /* 13 */
                                            1/36,    /* 14 */
                                            1/36,    /* 15 */
                                            2/36,    /* 16 */
                                            1/36,    /* 17 */
                                            1/36};   /* 18 */
  static const double C_S;

#endif

