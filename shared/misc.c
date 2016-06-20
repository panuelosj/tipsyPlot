#include <stdlib.h>
#include <stdio.h>
#include "../tipsyPlot.h"

// Brute force is probably faster, even if it has more comparisons
float findMaxVal(float* arrayIn, int len){
    int i;
    float max = arrayIn[0];
    for (i=1; i<len; i++)
        if (arrayIn[i] > max)
            max = arrayIn[i];
    return max;
}
float findMinVal(float* arrayIn, int len){
    int i;
    float min = arrayIn[0];
    for (i=1; i<len; i++)
        if (arrayIn[i] < min)
            min = arrayIn[i];
    return min;
}
