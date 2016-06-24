#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "tipsyPlot.h"
#include <plplot/plplot.h>

#define GAMMA 1.4
#define GAS_CONST 0.4


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
    const char genericTitle[] = "Shocktube";
    const int nsteps = 500, interval = 5;
    const int nout = (int)nsteps/interval;

    const float xmin = -7.0, xmax = 7.0;
    const int nbins = 200;
    const float dDelta = 0.01;

    // Physical Properties
    //const float gamma = 1.4;


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
    char command[500];
    int i, j;
    float minlocalt, maxlocalt;
    int filetime;

    /*
    DERIVED ARRAYS
    */
    int numvars = 6;
    derivedvar* plotvars = (derivedvar*)malloc(numvars*sizeof(derivedvar));
    initializeDerivedVar(&plotvars[0], "rho", "Density", "rho",  calc_rho, TYPE_GAS);
    initializeDerivedVar(&plotvars[1], "T", "Temperature", "temp", calc_temp, TYPE_GAS);
    initializeDerivedVar(&plotvars[2], "P", "Pressure", "p", calc_pressure, TYPE_GAS);
    initializeDerivedVar(&plotvars[3], "v_x", "Flow Velocity", "velx", calc_velx, TYPE_GAS);
    initializeDerivedVar(&plotvars[4], "p/rho^gamma", "Entropy(ish)", "entropy", calc_entropy, TYPE_GAS);
    initializeDerivedVar(&plotvars[5], "C", "Sound Speed", "cSound", calc_entropy, TYPE_GAS);

    /*
       ###    ##    ##    ###    ##       ##    ## ######## ####  ######
      ## ##   ###   ##   ## ##   ##        ##  ##     ##     ##  ##    ##
     ##   ##  ####  ##  ##   ##  ##         ####      ##     ##  ##
    ##     ## ## ## ## ##     ## ##          ##       ##     ##  ##
    ######### ##  #### ######### ##          ##       ##     ##  ##
    ##     ## ##   ### ##     ## ##          ##       ##     ##  ##    ##
    ##     ## ##    ## ##     ## ########    ##       ##    ####  ######

     ######   #######  ##       ##     ## ######## ####  #######  ##    ##
    ##    ## ##     ## ##       ##     ##    ##     ##  ##     ## ###   ##
    ##       ##     ## ##       ##     ##    ##     ##  ##     ## ####  ##
     ######  ##     ## ##       ##     ##    ##     ##  ##     ## ## ## ##
          ## ##     ## ##       ##     ##    ##     ##  ##     ## ##  ####
    ##    ## ##     ## ##       ##     ##    ##     ##  ##     ## ##   ###
     ######   #######  ########  #######     ##    ####  #######  ##    ##
    */
    double truerho[4] = {0.125, 3.485471e-1, 3.48571e-1, 0.125};
    double truevelx[4] = {1.0, 0.0, 0.0, -1.0};
    double truep[4] = {0.05, 2.448958e-1, 2.448958e-1, 0.05};
    double truetemp[4];
    double trueentropy[4];
    double truecsound[4];
    double arrrho[8], arrvelx[8], arrp[8], arrtemp[8], arrentropy[8], arrcsound[8];
    for (i = 0; i<4; i++){
        if (truerho[i] == 0.0) {
            truetemp[i] = 0.0;
            truecsound[i] = 0.0;
            trueentropy[i] = 0.0;
        }
        else {
            truetemp[i] = truep[i]/(truerho[i]*GAS_CONST);
            truecsound[i] = pow(GAMMA*truep[i]/truerho[i], 0.5);
            trueentropy[i] = truep[i]/(pow(truerho[i], GAMMA));
        }
        arrrho[2*i] = truerho[i];
        arrrho[2*i+1] = truerho[i];
        arrvelx[2*i] = truevelx[i];
        arrvelx[2*i+1] = truevelx[i];
        arrp[2*i] = truep[i];
        arrp[2*i+1] = truep[i];
        arrtemp[2*i] = truetemp[i];
        arrtemp[2*i+1] = truetemp[i];
        arrcsound[2*i] = truecsound[i];
        arrcsound[2*i+1] = truecsound[i];
        arrentropy[2*i] = trueentropy[i];
        arrentropy[2*i+1] = trueentropy[i];
    }
    double wavec[6] = {-5.591663e-1, -5.591663e-1, 0.0, 0.0, 5.591663e-1, 5.591663e-1};
    double analytic_xs[8] = {xmin, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, xmax};
    double** analytic_ys = (double**)malloc(numvars*sizeof(double));
    analytic_ys[0] = arrrho;
    analytic_ys[1] = arrtemp;
    analytic_ys[2] = arrp;
    analytic_ys[3] = arrvelx;
    analytic_ys[4] = arrentropy;
    analytic_ys[5] =  arrcsound;



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

    for (i=0; i<nout; i++){
        filetime = (i+1)*interval;                                              // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        printf("reading: %s\n", simname);
        snap = readTipsyStd(simname);
        pfile = profileCreate(snap, nbins, xmin, xmax, xpos);

        for (j=0; j<6; j++){
            analytic_xs[j+1] = wavec[j]*((float)filetime)*((float)dDelta);
        }

        // Individual plots
        for (j=0; j<numvars; j++){
            printf("\tplotting: %s\n", plotvars[j].title);
            calculateDerivedVarPoints(&plotvars[j], pfile, TYPE_GAS);
            sprintf(filenameout, "./test/%s/%s.scsShock.%05d.png", plotvars[j].shortname, plotvars[j].shortname, filetime);

            sprintf(title, "%s: %s t=%05d", genericTitle, plotvars[j].title, filetime);
            // setup plotting grid
                // plenv(xmin, xmax, ymin, ymax, just, axis);
                // just (axis scaling) - 0 = scaled independently
                // axis - 0 = draw box, ticks, and numeric tick labels

            //plsetopt("-o", filenameout);
            plsdev("pngcairo");
            plfontld(1);
            plsetopt("geometry", "810x670");
            plscolbg(255, 255, 255);
            plscol0(1, 0, 0, 0);
            plsfnam(filenameout);
            plinit();
            plenv(xmin, xmax, plotvars[j].ymin, plotvars[j].ymax, 0, 0);
            pllab("x", plotvars[j].label, title);
            plcol0(9);
            //plsym(plotvars[j].npoints, plotvars[j].points_xs, plotvars[j].points_ys, 229);
            plpoin(plotvars[j].npoints, plotvars[j].points_xs, plotvars[j].points_ys, -1);
            plcol0(1);
            plline(plotvars[j].nbins, plotvars[j].profile_xs, plotvars[j].profile_ys);
            plcol0(3);
            plline(8, analytic_xs, analytic_ys[j]);
            plend();
        }


    }

    // make gifs
    float delay = (float)300/(float)nout;
    for (i=0; i<numvars; i++){
        sprintf(command, "convert -layers optimize-transparency -delay %f ./test/%s/%s.*.png ./test/%s.animation.gif", delay, plotvars[i].shortname, plotvars[i].shortname, plotvars[i].shortname);
        printf("%s\n", command);
        system(command);
    }

    return 0;
}


// calc_var
float calc_rho(void* particle, int type){
    if (type == TYPE_GAS)
        return ((gas_particle*)particle)->rho;
    else
        errorCase(ERR_INVALID_ATTRIBUTE);
}
float calc_temp(void* particle, int type){
    if (type == TYPE_GAS)
        return ((gas_particle*)particle)->temp;
    else
        errorCase(ERR_INVALID_ATTRIBUTE);
}
float calc_pressure(void* particle, int type){
    if (type == TYPE_GAS)
        return (((gas_particle*)particle)->rho)*GAS_CONST*(((gas_particle*)particle)->temp);
    else
        errorCase(ERR_INVALID_ATTRIBUTE);
}
float calc_entropy(void* particle, int type){
    if (type == TYPE_GAS)
        return (GAS_CONST*((gas_particle*)particle)->temp/(pow(((gas_particle*)particle)->rho, GAMMA-1.0)));
    else
        errorCase(ERR_INVALID_ATTRIBUTE);
}
float calc_csound(void* particle, int type){
    if (type == TYPE_GAS)
        return pow(GAMMA*GAS_CONST*((gas_particle*)particle)->temp, 0.5);       //return pow(GAMMA*calc_pressure(particle, TYPE_GAS)/((gas_particle*)particle)->rho,0.5);
    else
        errorCase(ERR_INVALID_ATTRIBUTE);
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
