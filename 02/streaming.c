#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){

    int x = 0;
    int y = 0;
    int z = 0;
    int i = 0;
    int pos = 0;

    /* from 0 to 20 */
    for(z = 0; z < Q + 2; ++z){
        /* from 0 to 20 */
        for(y = 0; y < Q + 2; ++y){
            /* from 0 to 20 */
            for(x = 0; x < Q + 2; ++x){
                pos = ( z*xlength*xlength + y*xlength + x ); 
                /* from 0 to 20 */
                for (i = 0; i < Q + 2; ++i) {
                    /* probably not correct, need to check */
                    streamField[Q * pos + i] = collideField[Q * pos + i];
                }
            }
        }
    }

}
 
