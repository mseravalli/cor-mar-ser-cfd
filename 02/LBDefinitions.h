#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

/* TODO check LATTICEVELOCITIES */

  static const int LATTICEVELOCITIES[19][3] ={ { 0, -1, -1},    /* 0 */ 
                                               {-1,  0, -1},    /* 1 */                   
                                               { 0,  0, -1},    /* 2 */                   
                                               { 0,  0, -1},    /* 3 */                   
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
  static const double LATTICEWEIGHTS[19];
  static const double C_S;

#endif

