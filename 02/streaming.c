#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){

    int x = 0;
    int y = 0;
    int z = 0;
    int i = 0;
    int pos = 0;

    for(z = 0; z < xlength + 2; ++z){
        for(y = 0; y < xlength + 2; ++y){
            for(x = 0; x < xlength + 2; ++x){
                pos = ( z*xlength*xlength + y*xlength + x ); 
                /* from 0 to 19 */
                for (i = 0; i < Q; ++i) {
                    /* probably not correct, need to check */
                    streamField[Q * pos + i] = collideField[Q * pos + i];
                }
            }
        }
    }

}
 
