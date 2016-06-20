#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

profile* profileCreate(tipsy* tipsyIn, const int nbins, const float min, const float max, calc_x xs){
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
    profileOut->sim = tipsyIn;
    profileOut->nbins = nbins;
    profileOut->binwidth = (max - min)/((float)nbins);

    // Allocate and initialize bins
    profileOut->bin = (bin_attributes*)malloc(nbins*sizeof(bin_attributes));
    profileOut->gas = (gas_particle*)malloc(nbins*sizeof(gas_particle));
    profileOut->dark = (dark_particle*)malloc(nbins*sizeof(dark_particle));
    profileOut->star = (star_particle*)malloc(nbins*sizeof(star_particle));
    for (i=0; i<nbins; i++){
        profileOut->bin[i].min = min + i*profileOut->binwidth;
        profileOut->bin[i].max = min + (i+1)*profileOut->binwidth;
        profileOut->bin[i].ngas = 0;
        profileOut->bin[i].ndark = 0;
        profileOut->bin[i].nstar = 0;
        pFlopGas(&profileOut->gas[i], NULL, NULL, flopSetZero);
        pFlopDark(&profileOut->dark[i], NULL, NULL, flopSetZero);
        pFlopStar(&profileOut->star[i], NULL, NULL, flopSetZero);
    }

    // Allocate space for binned particle data and calculate the binned average
        // run through all the particles, finding the bin the particle
        // belongs to and adding it to the running total, then divide by the
        // total number of particles found for that bin
    if (tipsyIn->head->nsph != 0){
        profileOut->gas = (gas_particle*)malloc(nbins*sizeof(gas_particle));
        for (i=0; i < tipsyIn->head->nsph; i++){
            j = (int)floor((xs(tipsyIn, TYPE_GAS, i) - min)/((float)nbins));
            pFlopGas(&profileOut->gas[j], &profileOut->gas[j], &tipsyIn->gas[i], flopAdd);
            profileOut->bin[j].ngas ++;
        }
        vFlopGas(&profileOut->gas[j], &profileOut->gas[j], (float)profileOut->bin[j].ngas, flopDivide);
    } else profileOut->gas = NULL;
    if (tipsyIn->head->ndark != 0){
        profileOut->dark = (dark_particle*)malloc(nbins*sizeof(dark_particle));
        for (i=0; i < tipsyIn->head->ndark; i++){
            j = (int)floor((xs(tipsyIn, TYPE_DARK, i) - min)/((float)nbins));
            pFlopDark(&profileOut->dark[j], &profileOut->dark[j], &tipsyIn->dark[i], flopAdd);
            profileOut->bin[j].ndark ++;
        }
        vFlopDark(&profileOut->dark[j], &profileOut->dark[j], (float)profileOut->bin[j].ndark, flopDivide);
    } else profileOut->dark = NULL;
    if (tipsyIn->head->nstar != 0){
        profileOut->star = (star_particle*)malloc(nbins*sizeof(star_particle));
        for (i=0; i < tipsyIn->head->nstar; i++){
            j = (int)floor((xs(tipsyIn, TYPE_STAR, i) - min)/((float)nbins));
            pFlopStar(&profileOut->star[j], &profileOut->star[j], &tipsyIn->star[i], flopAdd);
            profileOut->bin[j].nstar ++;
        }
        vFlopStar(&profileOut->star[j], &profileOut->star[j], (float)profileOut->bin[j].nstar, flopDivide);
    } else profileOut->star = NULL;

    return profileOut;
}
