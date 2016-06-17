#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

profile* profileCreate(tipsy* tipsyIn, int nbins, float min, float max){
    /* Creates a new profile based on an input tipsy. Allocating all relevant
        memory based on the number of bins stated, as well as calculates the
        binned averages for all default variables. Structure for the profile is
        essentially exactly the same as a tipsy, using the same gas, dark, and
        star particle structs used, but each "particle" would represent the
        average for that bin.

        Parameters:
    */

    // Create
}
