#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "tipsyPlot.h"

float xpos(tipsy* tipsyIn, int type, int particle);

int main() {
    printf("Hello World\n");

    /*
    ########     ###    ########     ###    ##     ##  ######
    ##     ##   ## ##   ##     ##   ## ##   ###   ### ##    ##
    ##     ##  ##   ##  ##     ##  ##   ##  #### #### ##
    ########  ##     ## ########  ##     ## ## ### ##  ######
    ##        ######### ##   ##   ######### ##     ##       ##
    ##        ##     ## ##    ##  ##     ## ##     ## ##    ##
    ##        ##     ## ##     ## ##     ## ##     ##  ######
    */
    // FileIO params
    const char genericfilename[] = "shocktube";
    const char genericTitle[] = "Shocktube\n";
    const int nsteps = 50, interval = 5;
    const int nout = (int)nsteps/interval;

    const float minx = -7.0, maxx = 7.0;
    const int nbins = 200;

    // Physical Properties
    const float gamma = 1.4;


    /*
    ##     ##    ###    ########   ######
    ##     ##   ## ##   ##     ## ##    ##
    ##     ##  ##   ##  ##     ## ##
    ##     ## ##     ## ########   ######
     ##   ##  ######### ##   ##         ##
      ## ##   ##     ## ##    ##  ##    ##
       ###    ##     ## ##     ##  ######
    */
    // indexes
    char simname[100];
    int i, j, k;

    /*
    DERIVED ARRAYS
    */
    int numvars = 2;
    plottingvar* plotvars = (plottingvar*)malloc(numvars*sizeof(plottingvar));
    initializeDerivedVar(&plotvars[0], "rho", "Density",  calc_rho);








    /*
    ##     ##    ###    #### ##    ##
    ###   ###   ## ##    ##  ###   ##
    #### ####  ##   ##   ##  ####  ##
    ## ### ## ##     ##  ##  ## ## ##
    ##     ## #########  ##  ##  ####
    ##     ## ##     ##  ##  ##   ###
    ##     ## ##     ## #### ##    ##
    */

    // Find min-max bounds
    // first timestep, sets the min max bounds
    int filetime = (i+1)*interval;                                              // calculate current filename time
    sprintf(simname, "%s.%05d", genericfilename, filetime);                     // generate filename
    tipsy* snap = readTipsyStd(simname);
    profile* pfile = profileCreate(snap, nbins, minx, maxx, xpos);
    for (i=0; i<numvars; i++){
        plotvars[i].max = findMaxVal(float *arrayIn, int len)
    }


    // the rest, updates the min max bounds
    for (i=0; i<nout; i++){
        int filetime = (i+1)*interval;                                          // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        tipsy* snap = readTipsyStd(simname);
        profile* pfile = profileCreate(snap, nbins, minx, maxx, xpos);


    }
}


// calc_var
float calc_rho(bin_particle* bin){
    return bin->gas.rho;
}

// calc_bin
float xpos(tipsy* tipsyIn, int type, int p){
    switch (type){
        case TYPE_GAS:
            return tipsyIn->gas[p].pos[AXIS_X];
        case TYPE_DARK:
            return tipsyIn->dark[p].pos[AXIS_X];
        case TYPE_STAR:
            return tipsyIn->star[p].pos[AXIS_X];
    }
}
