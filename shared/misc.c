#include <stdlib.h>
#include <stdio.h>
#include "../tipsyPlot.h"

// Brute force is probably faster, even if it has more comparisons
double findMaxVal(double* arrayIn, int len){
    int i;
    double max = arrayIn[0];
    for (i=1; i<len; i++)
        if (arrayIn[i] > max)
            max = arrayIn[i];
    return max;
}
double findMinVal(double* arrayIn, int len){
    int i;
    double min = arrayIn[0];
    for (i=1; i<len; i++)
        if (arrayIn[i] < min)
            min = arrayIn[i];
    return min;
}
