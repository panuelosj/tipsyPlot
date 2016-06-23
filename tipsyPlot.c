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
    const char genericfilename[] = "sampledata/scsShock";
    const char genericTitle[] = "Shocktube\n";
    const int nsteps = 100, interval = 5;
    const int nout = (int)nsteps/interval;

    const float xmin = -7.0, xmax = 7.0;
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
    char simname[100];
    char filenameout[100];
    char title[100];
    int i, j, k;
    float minlocalt, maxlocalt;
    int filetime;

    /*
    DERIVED ARRAYS
    */
    int numvars = 3;
    derivedvar* plotvars = (derivedvar*)malloc(numvars*sizeof(derivedvar));
    initializeDerivedVar(&plotvars[0], "rho", "Density", "rho",  calc_rho, TYPE_GAS);
    initializeDerivedVar(&plotvars[1], "v_x", "Flow Velocity", "velx", calc_velx, TYPE_GAS);
    initializeDerivedVar(&plotvars[2], "T", "Temperature", "temp", calc_temp, TYPE_GAS);






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
    filetime = (0+1)*interval;                                                  // calculate current filename time
    sprintf(simname, "%s.%05d", genericfilename, filetime);                     // generate filename
    printf("reading: %s\n", simname);
    tipsy* snap = readTipsyStd(simname);
    profile* pfile = profileCreate(snap, nbins, xmin, xmax, xpos);
    for (j=0; j<numvars; j++){
        printf("\tcalculating: %s\n", plotvars[j].title);
        calculateDerivedVar(&plotvars[j], pfile, TYPE_GAS);
        plotvars[j].ymin = findMinVal(plotvars[j].profile_ys, plotvars[j].nbins);
        plotvars[j].ymax = findMaxVal(plotvars[j].profile_ys, plotvars[j].nbins);
    }
    tipsyDestroy(snap);
    // the rest, updates the min max bounds by checking if they are lower and higher
    for (i=1; i<nout; i++){
        filetime = (i+1)*interval;                                              // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        printf("reading: %s\n", simname);
        snap = readTipsyStd(simname);
        pfile = profileCreate(snap, nbins, xmin, xmax, xpos);
        for (j=0; j<numvars; j++){
            printf("\tcalculating: %s\n", plotvars[j].title);
            calculateDerivedVar(&plotvars[j], pfile, TYPE_GAS);
            minlocalt = findMinVal(plotvars[j].profile_ys, plotvars[j].nbins);
            maxlocalt = findMaxVal(plotvars[j].profile_ys, plotvars[j].nbins);
            if (minlocalt < plotvars[j].ymin)
                plotvars[j].ymin = minlocalt;
            if (maxlocalt > plotvars[j].ymax)
                plotvars[j].ymax = maxlocalt;
        }
        tipsyDestroy(snap);
        if (i%10 == 0) printf("done t=%d\n", i);
    }
    printf("\n\n\n");


    // extend boundaries
    for (j=0; j<numvars; j++){
        printf("adjusting boundaries: %s\n", plotvars[j].title);
        plotvars[j].ymax += 0.2*(plotvars[j].ymax - plotvars[j].ymin);
        plotvars[j].ymin -= 0.2*(plotvars[j].ymax - plotvars[j].ymin);
    }
    printf("\n\n\n");

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

    for (i=0; i<nout; i++){
        filetime = (i+1)*interval;                                          // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        printf("reading: %s\n", simname);
        snap = readTipsyStd(simname);
        pfile = profileCreate(snap, nbins, xmin, xmax, xpos);

        // Individual plots
        for (j=0; j<numvars; j++){
            printf("\tplotting: %s\n", plotvars[j].title);
            calculateDerivedVar(&plotvars[j], pfile, TYPE_GAS);
            sprintf(filenameout, "./test/%s/%s.scsShock%05d.png", plotvars[j].shortname, plotvars[j].shortname, filetime);
            // setup plotting grid
                // plenv(xmin, xmax, ymin, ymax, just, axis);
                // just (axis scaling) - 0 = scaled independently
                // axis - 0 = draw box, ticks, and numeric tick labels
            plsdev("pngcairo");
            plsetopt("geometry", "1080x810");
            plsetopt("-o", filenameout);
            plinit();
            plenv(xmin, xmax, plotvars[j].ymin, plotvars[j].ymax, 0, 0);
            plline(plotvars[j].nbins, plotvars[j].profile_xs, plotvars[j].profile_ys);
            plstring(plotvars[j].npoints, plotvars[j].points_xs, plotvars[j].points_ys, "#0x002e");
            sprintf(title, "%s t=%05d", plotvars[j].title, filetime);
            pllab("x", plotvars[j].label, title);
            plend();
        }


    }



    return 0;
}


// calc_var
float calc_rho(void* particle, int type){
    switch (type) {
        case TYPE_GAS:
            return ((gas_particle*)particle)->rho;
            break;
        default:
            errorCase(ERR_INVALID_ATTRIBUTE);
    }
}
float calc_velx(void* particle, int type){
    switch (type) {
        case TYPE_GAS:
            return ((gas_particle*)particle)->vel[AXIS_X];
            break;
        case TYPE_DARK:
            return ((dark_particle*)particle)->vel[AXIS_X];
            break;
        case TYPE_STAR:
            return ((star_particle*)particle)->vel[AXIS_X];
            break;
    }
}
float calc_temp(void* particle, int type){
    switch (type) {
        case TYPE_GAS:
            return ((gas_particle*)particle)->temp;
            break;
        default:
            errorCase(ERR_INVALID_ATTRIBUTE);
    }
}

// calc_bin
float xpos(tipsy* tipsyIn, int type, int p){
    switch (type){
        case TYPE_GAS:
            return tipsyIn->gas[p].pos[AXIS_X];
            break;
        case TYPE_DARK:
            return tipsyIn->dark[p].pos[AXIS_X];
            break;
        case TYPE_STAR:
            return tipsyIn->star[p].pos[AXIS_X];
            break;
        default:
            errorCase(ERR_UNKNOWN_PARTICLE);
    }
}
