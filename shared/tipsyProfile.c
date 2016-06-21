#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

void initializeDerivedVar(derivedvar* variable, const char label[], const char title[], calc_var equation){
    variable->label = label;
    variable->title = title;
    variable->equation = equation;
}

void calculateDerivedVar(derivedvar* variable, profile* profileIn){
    /* Calculates an array of an arbitrary variable, using the equation
        provided in the plottingvar struct, and fills in the values in the
        array pointed to by variable->derived_array. The number of bins is read
        from the input profile and inherited by the derivedvar

        Parameters:
            derivedvar* variable    - pointer to the new variable to be filled in
            profile* profileIn      - profile of nonderived variables whose
                                        values will be used to calculate the
                                        new derived variable
    */
    int i;
    variable->nbins = profileIn->nbins;                                         // inherit nbins
    for (i=0; i<variable->nbins; i++)                                           // calculate value of new variable for all bins in profile
        variable->derived_array[i] = variable->equation(&(profileIn->bin[i]));
}

profile* profileCreate(tipsy* tipsyIn, const int nbins, const float min, const float max, calc_bin xs){
    /* Creates a new profile based on an input tipsy. Allocating all relevant
        memory based on the number of bins stated, as well as calculates the
        binned averages for all default variables. Structure for the profile is
        essentially exactly the same as a tipsy, using the same gas, dark, and
        star particle structs used, but each "particle" would represent the
        average for that bin.

        Currently allows creation of an empty profile. Should I throw an error
        for this?
            - in empty particle case, the bins are still allocated and filled
                with zero, but calculations are not performed to avoid division
                by zero error

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

    // Allocate and initialize bins and particles
    profileOut->bin = (bin_particle*)malloc(nbins*sizeof(bin_particle));
    for (i=0; i<nbins; i++){
        profileOut->bin[i].min = min + i*profileOut->binwidth;
        profileOut->bin[i].max = min + (i+1)*profileOut->binwidth;
        profileOut->bin[i].ngas = 0;
        profileOut->bin[i].ndark = 0;
        profileOut->bin[i].nstar = 0;
        pFlopGas(&(profileOut->bin[i].gas), NULL, NULL, flopSetZero);
        pFlopDark(&(profileOut->bin[i].dark), NULL, NULL, flopSetZero);
        pFlopStar(&(profileOut->bin[i].star), NULL, NULL, flopSetZero);
    }

    // Allocate space for binned particle data and calculate the binned average
        // run through all the particles, finding the bin the particle
        // belongs to and adding it to the running total, then divide by the
        // total number of particles found for that bin
    if (tipsyIn->head->nsph != 0){
        for (i=0; i < tipsyIn->head->nsph; i++){
            j = (int)floor((xs(tipsyIn, TYPE_GAS, i) - min)/((float)nbins));
            pFlopGas(&(profileOut->bin[j].gas), &(profileOut->bin[j].gas), &(tipsyIn->gas[i]), flopAdd);
            profileOut->bin[j].ngas ++;
        }
        vFlopGas(&(profileOut->bin[j].gas), &(profileOut->bin[j].gas), (float)profileOut->bin[j].ngas, flopDivide);
    }
    if (tipsyIn->head->ndark != 0){
        for (i=0; i < tipsyIn->head->ndark; i++){
            j = (int)floor((xs(tipsyIn, TYPE_DARK, i) - min)/((float)nbins));
            pFlopDark(&(profileOut->bin[j].dark), &(profileOut->bin[j].dark), &(tipsyIn->dark[i]), flopAdd);
            profileOut->bin[j].ndark ++;
        }
        vFlopDark(&(profileOut->bin[j].dark), &(profileOut->bin[j].dark), (float)profileOut->bin[j].ndark, flopDivide);
    }
    if (tipsyIn->head->nstar != 0){
        for (i=0; i < tipsyIn->head->nstar; i++){
            j = (int)floor((xs(tipsyIn, TYPE_STAR, i) - min)/((float)nbins));
            pFlopStar(&(profileOut->bin[j].star), &(profileOut->bin[j].star), &(tipsyIn->star[i]), flopAdd);
            profileOut->bin[j].nstar ++;
        }
        vFlopStar(&(profileOut->bin[j].star), &(profileOut->bin[j].star), (float)profileOut->bin[j].nstar, flopDivide);
    }

    return profileOut;
}

void profileDestroy(profile* profileIn){
    free(profileIn->bin);
    free(profileIn);
}
