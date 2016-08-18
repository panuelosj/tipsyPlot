#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "tipsyPlot.h"
#include <plplot/plplot.h>
#include <libconfig.h>

#define GAMMA 1.4
#define GAS_CONST 0.4
#define SLICE_MIN -0.5
#define SLICE_MAX 0.5

int main(int argc, char* argv[]){
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

    const char* genericfilename = NULL;
    const char* genericfileout = NULL;
    const char* genericTitle = NULL;
    int nsteps, interval, nbins;
    double dxmin, dxmax, dymin, dymax, dmaxrho, ddDelta;
    float xmin, xmax, ymin, ymax, maxrho, dDelta;


    config_t cfg;
    config_init(&cfg);
    if (argv[1] == NULL || argc == 1){
        system("cat nagato");
        printf("Error: No config file was indicated\n\tPlease enter a config file.\n");
        printf("\tThe following config files were found in your working directory:\n");
        system("ls *.config");
        exit(-1);
    }
    if (config_read_file(&cfg, argv[1]) == 0){
        system("cat nagato");
        printf("Error: cannot find file %s\n", argv[1]);
        exit(-1);
    }
    else {
        printf("reading config file: %s\n", argv[1]);
        config_lookup_string(&cfg, "genericfilename", &genericfilename);
        printf("\tgenericfilename = %s\n", genericfilename);
        config_lookup_string(&cfg, "genericfileout", &genericfileout);
        printf("\tgenericfileout = %s\n", genericfileout);
        config_lookup_string(&cfg, "genericTitle", &genericTitle);
        printf("\tgenericTitle = %s\n", genericTitle);
        config_lookup_int(&cfg, "nsteps", &nsteps);
        config_lookup_int(&cfg, "interval", &interval);
        printf("\tnsteps = %d\tinterval = %d\n", nsteps, interval);
        config_lookup_int(&cfg, "nbins", &nbins);
        printf("\tnbins = %d\n", nbins);
        config_lookup_float(&cfg, "xmin", &dxmin);
        config_lookup_float(&cfg, "xmax", &dxmax);
        config_lookup_float(&cfg, "ymin", &dymin);
        config_lookup_float(&cfg, "ymax", &dymax);
        xmin=(float)dxmin; xmax=(float)dxmax; ymin=(float)dymin; ymax=(float)dymax;
        printf("\txmin = %f\txmax = %f\n\tymin = %f\tymax = %f\n", xmin, xmax, ymin, ymax);
        config_lookup_float(&cfg, "maxrho", &dmaxrho);
        maxrho=(float)dmaxrho;
        printf("\tmaxrho = %f\n", maxrho);
        config_lookup_float(&cfg, "dDelta", &ddDelta);
        dDelta=(float)ddDelta;
        printf("\tdDelta = %f\n", dDelta);
    }

    const int nout = (int)nsteps/interval;

    const double cmapi[2] = {0.0, 1.0};                                        // colour map intensity
    const double cmapr[2] = {1.0, 1.0};                                        // red, left
    const double cmapg[2] = {0.0, 0.0};                                        // green, left
    const double cmapb[2] = {1.0, 0.0};                                        // blue, left
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
    int i, j, k;
    int filetime;

    tipsy* snap;
    profile* pfile;

    /*
    DERIVED ARRAYS
    */
    int numvars = 1;
    derivedvar* plotvars = (derivedvar*)malloc(numvars*sizeof(derivedvar));
    initializeDerivedVar(&plotvars[0], "y", "slice", "slice", calc_y, TYPE_GAS);


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
    printf("creating plots now\n");
    // delete all pictures
    printf("rm test/slicerho/*.png -v\n");
    system("rm test/slicerho/slicerho.*.png -v");

    for (i=0; i<nout; i++){
        filetime = (i+1)*interval;                                              // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        printf("reading: %s\n", simname);

        snap = readTipsyStd(simname);                                           // prepare cropped tipsys
        tipsyCrop(snap, cropThinSliceZ);
        pfile = profileCreate(snap, nbins, xmin, xmax, xpos);

        for (j=0; j<numvars; j++){
            printf("\tplotting: slicerho\n");
            sprintf(filenameout, "./test/slicerho/slicerho.%s.%05d.png", genericfileout, filetime);
            sprintf(title, "%s: %s t=%.4f", genericTitle, plotvars[j].title, (float)filetime*dDelta);

            plsdev("pngcairo");
            plfontld(1);
            plsetopt("geometry", "810x670");
            plsetopt("dpi", "810x670");
            plscolbg(255, 255, 255);
            plscol0(1, 0, 0, 0);
            plsfnam(filenameout);
            plinit();
            plenv(xmin, xmax, ymin, ymax, 0, 0);
            pllab("x", "y", title);
            // plot points
            calculateDerivedVarPoints(&plotvars[j], pfile, TYPE_GAS);
            plscmap1l(1, 2, cmapi, cmapr, cmapg, cmapb, 0);
            for (k=0; k < snap->head->nsph; k++) {
                if ((snap->gas[k].rho/maxrho) > 1)
                    plcol1(1.0);
                else
                    plcol1((snap->gas[k].rho)/maxrho);
                plstring(1, &plotvars[j].points_xs[k], &plotvars[j].points_ys[k], "#(210)");
            }
            plend();
        }

        tipsyDestroy(snap);
        profileDestroy(pfile);
    }

    // make gifs
    float delay = (float)1000/(float)nout;
    sprintf(command, "convert -layers optimize -delay %f ./test/slicerho/slicerho.*.png ./test/slicerho.animation.gif", delay);
    printf("%s\n", command);
    system(command);

    return 0;
}


// slice crops
int cropThinSliceZ(void* particle, int type){
    if (type == TYPE_GAS){
        if (((gas_particle*)particle)->pos[AXIS_Z]>SLICE_MIN && ((gas_particle*)particle)->pos[AXIS_Z]<SLICE_MAX)
            return 1;
        else
            return 0;
    }
    else {
        return 1;
    }
}



// calc_var
float calc_y(void* particle, int type){
    if (type == TYPE_GAS)
        return ((gas_particle*)particle)->pos[AXIS_Y];
    else
        errorCase(ERR_INVALID_ATTRIBUTE);
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
