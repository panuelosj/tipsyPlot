#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

profile* profileCreate(tipsy* tipsyIn, const int nbins, const float min, const float max, float (*calc_x)(tipsy*, int type, int particle)){
    /* Creates a new profile based on an input tipsy. Allocating all relevant
        memory based on the number of bins stated, as well as calculates the
        binned averages for all default variables. Structure for the profile is
        essentially exactly the same as a tipsy, using the same gas, dark, and
        star particle structs used, but each "particle" would represent the
        average for that bin.

        Currently allows creation of an empty profile. Should I throw an error
        for this?

        Parameters:
            const tipsy* tipsyIn    - tipsy snapshot to make a profile of
            const int nbins         - number of bins to subdivide the profile into
            const float min         - minimum value for the lower bin extreme
            const float max         - maximum value for the uppen bin extreme
            float (*calc_x(tipsy*)) - a pointer to a function that accepts the
                                        tipsy snapshot and returns the value to
                                        be binned

        ToDo:
            - throw error for empty tipsyIn case?
    */

    // indexing vars
    int i, j;

    // Create object (pointer to a struct)
    profile* profileOut = (profile*)malloc(sizeof(profile));
    profileOut->nbins = nbins;
    profileOut->binwidth = (max - min)/((float)nbins);
    profileOut->sim = tipsyIn;

    // Allocate space
    if (tipsyIn->head->nbodies != 0){
        profileOut->bin = (bin_particle*)malloc(nbins*sizeof(bin_particle));
    } else profileOut->bin = NULL;

    // initialize bin values
    for(i=0; i<nbins; i++){
        profileOut->bin[i].min = min + i*profileOut->binwidth;
        profileOut->bin[i].max = min + (i+1)*profileOut->binwidth;
        profileOut->bin[i].ngas = 0;
        profileOut->bin[i].ndark = 0;
        profileOut->bin[i].nstar = 0;
    }

    // calculate bin values
        // (run through all particles and update the bin values accordingly)

    for (i=0; i < tipsyIn->head->nsph; i++){
        j = (int)floor(((*calc_x)(tipsyIn, TYPE_GAS, i) - min)/((float)nbins));

    }
    for (i=0; i < tipsyIn->head->ndark; i++){
        j = (int)floor(((*calc_x)(tipsyIn, TYPE_DARK, i) - min)/((float)nbins));

    }
    for (i=0; i < tipsyIn->head->nstar; i++){
        j = (int)floor(((*calc_x)(tipsyIn, TYPE_STAR, i) - min)/((float)nbins));

    }


    return profileOut;
}
