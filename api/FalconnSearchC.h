/*
 * FalconnSearchC.h
 *
 *  Created on: 20 Apr 2017
 *      Author: joni
 */

#ifndef SOURCE_FALCONNSEARCHC_H_
#define SOURCE_FALCONNSEARCHC_H_

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    struct Result {
        int* Result;
        int ResultCount;
    };
    typedef struct Result Result;

    int buildIndexFromData(void* data, int featureCount, int featureBytes);

    Result search(void* searchData, int dataLength);

#ifdef __cplusplus
}
#endif

#endif
