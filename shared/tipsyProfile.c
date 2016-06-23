#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

void initializeDerivedVar(derivedvar* variable, const char label[], const char title[], const char shortname[], calc_var equation, int type){
    /* Initalizes the values of the derived var with the given parameters.
        Does not allocate any memory, but NULL initializes the pointer to the
        data array and sets nbins to zero to be able to detect array
        non-existence and dynamically allocate later (particularly in
        calculateDerivedVar()).

        Parameters:
            derivedvar* variable    - pointer to the new variable to be initialized
            const char label[]      - short name or variable symbol for graph axes
            const char title[]      - formal name to be shown as graph title
            calc_var equation       - pointer to function representing equation
                                        to calculate the variable's value
    */
    variable->title = title;
    variable->label = label;
    variable->shortname = shortname;
    variable->equation = equation;
    variable->type = type;
    variable->ymin = nanf("uninitialized");
    variable->ymax = nanf("uninitialized");

    // NULL initialize to be able to check if the profile_ys needs to be
        // allocated or not upon calculation
    variable->nbins = 0;
    variable->npoints = 0;
    variable->profile_xs = NULL;
    variable->profile_ys = NULL;
    variable->points_xs = NULL;
    variable->points_ys = NULL;
}

void calculateDerivedVar(derivedvar* variable, profile* profileIn, int type){
    /* Calculates an array of an arbitrary variable, using the equation
        provided in the plottingvar struct, and fills in the values in the
        array pointed to by variable->profile_ys. The number of bins is read
        from the input profile and inherited by the derivedvar.

        Parameters:
            derivedvar* variable    - pointer to the new variable to be filled in
            profile* profileIn      - profile of nonderived variables whose
                                        values will be used to calculate the
                                        new derived variable
    */
    int i;
    variable->type = type;
    // Profile
    // Allocate memory for the variable's x- and y-values to be recorded into
        // realloc should be able to tell if the pointer needs to be freed or is NULL-initialized
    variable->profile_xs = (double*)realloc(variable->profile_xs, profileIn->nbins*sizeof(double));
    variable->profile_ys = (double*)realloc(variable->profile_ys, profileIn->nbins*sizeof(double));

    variable->nbins = profileIn->nbins;                                         // inherit nbins
    // don't combine for loops to write straight down the arrays instead of moving back and forth
    for (i=0; i<variable->nbins; i++) {                                         // calculate xs
        variable->profile_xs[i] = (double)profileIn->bin[i].xval;
    }
    for (i=0; i<variable->nbins; i++) {                                         // calculate value of new variable for all bins in profile
        switch (type) {
            case TYPE_GAS:
                variable->profile_ys[i] = (double)variable->equation(&(profileIn->bin[i].gas), type);
                break;
            case TYPE_DARK:
                variable->profile_ys[i] = (double)variable->equation(&(profileIn->bin[i].dark), type);
                break;
            case TYPE_STAR:
                variable->profile_ys[i] = (double)variable->equation(&(profileIn->bin[i].star), type);
                break;
            default:
                errorCase(ERR_UNKNOWN_PARTICLE);
        }
    }
}

void calculateDerivedVarPoints(derivedvar* variable, profile* profileIn, int type){
    /* Calculates an array of an arbitrary variable, using the equation
        provided in the plottingvar struct, and fills in the values in the
        array pointed to by variable->profile_ys. The number of bins is read
        from the input profile and inherited by the derivedvar.

        Parameters:
            derivedvar* variable    - pointer to the new variable to be filled in
            profile* profileIn      - profile of nonderived variables whose
                                        values will be used to calculate the
                                        new derived variable
    */
    int i;
    variable->type = type;
    // Profile
    // Allocate memory for the variable's x- and y-values to be recorded into
        // realloc should be able to tell if the pointer needs to be freed or is NULL-initialized
    variable->profile_xs = (double*)realloc(variable->profile_xs, profileIn->nbins*sizeof(double));
    variable->profile_ys = (double*)realloc(variable->profile_ys, profileIn->nbins*sizeof(double));

    variable->nbins = profileIn->nbins;                                         // inherit nbins
    // don't combine for loops to write straight down the arrays instead of moving back and forth
    for (i=0; i<variable->nbins; i++) {                                         // calculate xs
        variable->profile_xs[i] = (double)profileIn->bin[i].xval;
    }
    for (i=0; i<variable->nbins; i++) {                                         // calculate value of new variable for all bins in profile
        switch (type) {
            case TYPE_GAS:
                variable->profile_ys[i] = (double)variable->equation(&(profileIn->bin[i].gas), type);
                break;
            case TYPE_DARK:
                variable->profile_ys[i] = (double)variable->equation(&(profileIn->bin[i].dark), type);
                break;
            case TYPE_STAR:
                variable->profile_ys[i] = (double)variable->equation(&(profileIn->bin[i].star), type);
                break;
            default:
                errorCase(ERR_UNKNOWN_PARTICLE);
        }
    }

    // Points
    // Allocate memory for x- and y-values and fill in values
    switch (type){
        case TYPE_GAS:
            variable->points_xs = (double*)realloc(variable->points_xs, profileIn->sim->head->nsph*sizeof(double));
            variable->npoints = profileIn->sim->head->nsph;
            for (i=0; i<profileIn->sim->head->nsph; i++)
                variable->points_xs[i] = (double)(profileIn->eqbin(profileIn->sim, type, i));
            for (i=0; i<profileIn->sim->head->nsph; i++)
                variable->points_ys[i] = (double)(variable->equation(&(profileIn->sim->gas[i]), type));
            break;
        case TYPE_DARK:
            variable->points_xs = (double*)realloc(variable->points_xs, profileIn->sim->head->ndark*sizeof(double));
            variable->npoints = profileIn->sim->head->ndark;
            for (i=0; i<profileIn->sim->head->ndark; i++)
                variable->points_xs[i] = (double)(profileIn->eqbin(profileIn->sim, type, i));
            for (i=0; i<profileIn->sim->head->ndark; i++)
                variable->points_ys[i] = (double)(variable->equation(&(profileIn->sim->dark[i]), type));
            break;
        case TYPE_STAR:
            variable->points_xs = (double*)realloc(variable->points_xs, profileIn->sim->head->nstar*sizeof(double));
            variable->npoints = profileIn->sim->head->nstar;
            for (i=0; i<profileIn->sim->head->nstar; i++)
                variable->points_xs[i] = (double)(profileIn->eqbin(profileIn->sim, type, i));
            for (i=0; i<profileIn->sim->head->nstar; i++)
                variable->points_ys[i] = (double)(variable->equation(&(profileIn->sim->star[i]), type));
            break;
        default:
            errorCase(ERR_UNKNOWN_PARTICLE);
    }
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
    profileOut->eqbin = xs;
    profileOut->binwidth = (max - min)/((float)nbins);

    // Allocate and initialize bins and particles
    profileOut->bin = (bin_particle*)malloc(nbins*sizeof(bin_particle));
    for (i=0; i<nbins; i++){
        profileOut->bin[i].xmin = min + ((float)i)*profileOut->binwidth;
        profileOut->bin[i].xmax = min + (((float)i)+1.0)*profileOut->binwidth;
        profileOut->bin[i].xval = min + (((float)i)+0.5)*profileOut->binwidth;
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
            j = (int)floor((xs(tipsyIn, TYPE_GAS, i) - min)/(profileOut->binwidth));
            if (j < 0) j = 0; if (j >= nbins) j = nbins-1;
            pFlopGas(&(profileOut->bin[j].gas), &(profileOut->bin[j].gas), &(tipsyIn->gas[i]), flopAdd);
            profileOut->bin[j].ngas ++;
        }
        for (j=0; j < nbins; j++)
            vFlopGas(&(profileOut->bin[j].gas), &(profileOut->bin[j].gas), (float)profileOut->bin[j].ngas, flopDivide);
    }
    if (tipsyIn->head->ndark != 0){
        for (i=0; i < tipsyIn->head->ndark; i++){
            j = (int)floor((xs(tipsyIn, TYPE_DARK, i) - min)/(profileOut->binwidth));
            if (j < 0) j = 0; if (j >= nbins) j = nbins-1;
            pFlopDark(&(profileOut->bin[j].dark), &(profileOut->bin[j].dark), &(tipsyIn->dark[i]), flopAdd);
            profileOut->bin[j].ndark ++;
        }
        for (j=0; j < nbins; j++)
            vFlopDark(&(profileOut->bin[j].dark), &(profileOut->bin[j].dark), (float)profileOut->bin[j].ndark, flopDivide);
    }
    if (tipsyIn->head->nstar != 0){
        for (i=0; i < tipsyIn->head->nstar; i++){
            j = (int)floor((xs(tipsyIn, TYPE_STAR, i) - min)/(profileOut->binwidth));
            if (j < 0) j = 0; if (j >= nbins) j = nbins-1;
            pFlopStar(&(profileOut->bin[j].star), &(profileOut->bin[j].star), &(tipsyIn->star[i]), flopAdd);
            profileOut->bin[j].nstar ++;
        }
        for (j=0; j < nbins; j++)
            vFlopStar(&(profileOut->bin[j].star), &(profileOut->bin[j].star), (float)profileOut->bin[j].nstar, flopDivide);
    }

    return profileOut;
}

void profileDestroy(profile* profileIn){
    free(profileIn->bin);
    free(profileIn);
}
