#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "tipsyPlot.h"
#include <plplot/plplot.h>

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
    int filetime;

    /*
    DERIVED ARRAYS
    */
    int numvars = 2;
    derivedvar* plotvars = (derivedvar*)malloc(numvars*sizeof(derivedvar));
    initializeDerivedVar(&plotvars[0], "rho", "Density",  calc_rho);
    initializeDerivedVar(&plotvars[1], "v_x", "Flow Velocity", calc_velx);
    initializeDerivedVar(&plotvars[2], "T", "Temperature", calc_temp);








    /*
    ########   #######  ##     ## ##    ## ########   ######
    ##     ## ##     ## ##     ## ###   ## ##     ## ##    ##
    ##     ## ##     ## ##     ## ####  ## ##     ## ##
    ########  ##     ## ##     ## ## ## ## ##     ##  ######
    ##     ## ##     ## ##     ## ##  #### ##     ##       ##
    ##     ## ##     ## ##     ## ##   ### ##     ## ##    ##
    ########   #######   #######  ##    ## ########   ######
    */

    // Find min-max bounds
    // first timestep, sets the min max bounds
    filetime = (0+1)*interval;                                              // calculate current filename time
    sprintf(simname, "%s.%05d", genericfilename, filetime);                     // generate filename
    tipsy* snap = readTipsyStd(simname);
    profile* pfile = profileCreate(snap, nbins, minx, maxx, xpos);
    for (j=0; j<numvars; j++){
        calculateDerivedVar(&plotvars[j], pfile);
        plotvars[j].max = findMaxVal(plotvars[j].derived_array, plotvars[j].nbins);
        plotvars[j].min = findMinVal(plotvars[j].derived_array, plotvars[j].nbins);
    }
    // the rest, updates the min max bounds by checking if they are lower and higher
    for (i=1; i<nout; i++){
        filetime = (i+1)*interval;                                          // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        tipsy* snap = readTipsyStd(simname);
        profile* pfile = profileCreate(snap, nbins, minx, maxx, xpos);
        for (j=0; j<numvars; j++){
            calculateDerivedVar(&plotvars[j], pfile);
            if (findMaxVal(plotvars[j].derived_array, plotvars[j].nbins) > plotvars[j].max)
                plotvars[j].max = findMaxVal(plotvars[j].derived_array, plotvars[j].nbins);
            else if (findMinVal(plotvars[j].derived_array, plotvars[j].nbins) < plotvars[j].min)
                plotvars[j].min = findMinVal(plotvars[j].derived_array, plotvars[j].nbins);
        }
        if (i%10 == 0) printf("done t=%d\n", i);
    }
    // extend boundaries
    for (j=0; j<numvars; j++){
        plotvars[j].max += 0.2*(plotvars[j].max - plotvars[j].min);
        plotvars[j].min -= 0.2*(plotvars[j].max - plotvars[j].min);
    }

    /*
    ########  ##        #######  ########
    ##     ## ##       ##     ##    ##
    ##     ## ##       ##     ##    ##
    ########  ##       ##     ##    ##
    ##        ##       ##     ##    ##
    ##        ##       ##     ##    ##
    ##        ########  #######     ##
    */
    printf("found bounds, creating plots now\n");
    // Initialize PLplot
    plsdev("png");
    plinit();

    for (i=0; i<nout; i++){
        filetime = (i+1)*interval;                                          // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        tipsy* snap = readTipsyStd(simname);
        profile* pfile = profileCreate(snap, nbins, minx, maxx, xpos);



    }



    return 0;
}


// calc_var
float calc_rho(bin_particle* bin){
    return bin->gas.rho;
}
float calc_velx(bin_particle* bin){
    return bin->gas.vel[AXIS_X];
}
float calc_temp(bin_particle* bin){
    return bin->gas.temp;
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
    return 0;
}
