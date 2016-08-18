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
    int i, j, k, l;

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
    double dxmin, dxmax, dymin, dymax, dmaxdumpval, dmindumpval, ddDelta;
    float xmin, xmax, ymin, ymax, maxdumpval, mindumpval, dDelta;

    // read config file
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
        config_lookup_float(&cfg, "dDelta", &ddDelta);
        dDelta=(float)ddDelta;
        printf("\tdDelta = %f\n", dDelta);
    }

    // read config file arrays
    // find what to plot - scalars
    int numdumps;
    config_lookup_int(&cfg, "numdumps", &numdumps);                             // get number of outputs to dump
    printf("\tnumdumps = %d\n", numdumps);
    char** plotdumps = (char**)malloc(numdumps*sizeof(char*));                 // allocate pointer array
    const config_setting_t* configarraybuff;

    configarraybuff = config_lookup(&cfg, "plotdumps");
    printf("\tplotdumps =");
    for (i=0; i<numdumps; i++){
        plotdumps[i] = (char*)config_setting_get_string_elem(configarraybuff, i);  // read in the output dump strings
        printf("\t%s", plotdumps[i]);
    }
    printf("\n");

    // find what to plot - vectors (each individual axis counts as a separate vector)
    int numvecs;
    config_lookup_int(&cfg, "numvecs", &numvecs);                               // get number of vector outputs to dump
    printf("\tnumvecs = %d\n", numvecs);
    char** plotvecs = (char**)malloc(numvecs*sizeof(char*));                    // allocate pointer array

    configarraybuff = config_lookup(&cfg, "plotvecs");
    printf("\tplotvecs =");
    for (i=0; i<numvecs; i++){
        plotvecs[i] = (char*)config_setting_get_string_elem(configarraybuff, i);  // read in the output dump strings
        printf("\t%s", plotvecs[i]);
    }
    printf("\n");



    const int nout = (int)nsteps/interval;

    const double cmapi[2] = {0.0, 1.0};                                        // colour map intensity
    const double cmapr[2] = {0.0, 1.0};                                        // red
    const double cmapg[2] = {0.0, 0.0};                                        // green
    const double cmapb[2] = {1.0, 0.0};                                        // blue
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
    char dumpname[100];
    // min and max variable name holder are shared with scalar and vectors
    char maxdumpname[100];
    char mindumpname[100];
    char filenameout[100];
    char title[100];
    char command[500];
    int filetime;

    tipsy* snap;
    profile* pfile;
    float* dump;

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
    //printf("rm test/slicerho/*.png -v\n");
    //system("rm test/slicerho/slicerho.*.png -v");

    // colorbar variables
    int cbar_label_opts[] = {PL_COLORBAR_LABEL_BOTTOM};
    const char* cbar_axis_opts[] = {"bcvtm"};
    char** cbar_labels;
    double cbar_axis_ticks[] = {0.0};
    int cbar_axis_subticks[] = {0};
    int cbar_n_values[] = {2};
    double cbar_valuex[2] = {0.0, 1.0};
    double* cbar_values[] = {cbar_valuex};
    double cbar_width;
    double cbar_height;

    for (i=0; i<nout; i++){
        filetime = (i+1)*interval;                                              // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        printf("reading: %s\n", simname);

        snap = readTipsyStd(simname);                                           // prepare cropped tipsys
        //tipsyCrop(snap, cropThinSliceZ);
        pfile = profileCreate(snap, nbins, xmin, xmax, xpos);

        for (j=0; j<numvars; j++){
            for (k=0; k<numdumps+numvecs; k++){
                // output filename and plot title
                // if scalar
                if (k<numdumps){
                    printf("\tplotting: %s\n", plotdumps[k]);
                    sprintf(filenameout, "./test/%s/%s.%s.%05d.png", plotdumps[k], plotdumps[k], genericfileout, filetime);
                    sprintf(title, "%s: %s %s t=%.4f", genericTitle, plotdumps[k], plotvars[j].title, (float)filetime*dDelta);
                } // if vector
                else {
                    printf("\tplotting: %s\n", plotvecs[k-numdumps]);
                    sprintf(filenameout, "./test/%s/%s.%s.%05d.png", plotvecs[k-numdumps], plotvecs[k-numdumps], genericfileout, filetime);
                    sprintf(title, "%s: %s %s t=%.4f", genericTitle, plotvecs[k-numdumps], plotvars[j].title, (float)filetime*dDelta);
                }

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

                // READ IN DUMP DATA
                // if scalar
                if (k<numdumps){
                    sprintf(dumpname, "%s.%s", simname, plotdumps[k]);
                    dump = readTipsyDumpScalar(dumpname);
                    // grab min and max values from config
                    sprintf(maxdumpname, "max%s", plotdumps[k]);
                    config_lookup_float(&cfg, maxdumpname, &dmaxdumpval);
                    maxdumpval=(float)dmaxdumpval;
                    sprintf(mindumpname, "min%s", plotdumps[k]);
                    config_lookup_float(&cfg, mindumpname, &dmindumpval);
                    mindumpval=(float)dmindumpval;
                } // if vector
                else {
                    sprintf(dumpname, "%s.%s", simname, plotvecs[k-numdumps]+sizeof(char));
                    if (plotvecs[k-numdumps][0] == 'x')
                        dump = readTipsyDumpVector(dumpname, AXIS_X);
                    else if (plotvecs[k-numdumps][0] == 'y')
                        dump = readTipsyDumpVector(dumpname, AXIS_Y);
                    else if (plotvecs[k-numdumps][0] == 'z')
                        dump = readTipsyDumpVector(dumpname, AXIS_Z);
                    else if (plotvecs[k-numdumps][0] == 'm')
                        dump = readTipsyDumpVectorMagnitude(dumpname);
                    else {
                        printf("Error: no axis given for vector dump file %s\nExiting program\n", plotvecs[k-numdumps]);
                        exit(-1);
                    }
                    // grab min and max values from config
                    sprintf(maxdumpname, "max%s", plotvecs[k-numdumps]);
                    config_lookup_float(&cfg, maxdumpname, &dmaxdumpval);
                    maxdumpval=(float)dmaxdumpval;
                    sprintf(mindumpname, "min%s", plotvecs[k-numdumps]);
                    config_lookup_float(&cfg, mindumpname, &dmindumpval);
                    mindumpval=(float)dmindumpval;
                }

                // plot points
                calculateDerivedVarPoints(&plotvars[j], pfile, TYPE_GAS);
                plscmap1l(1, 2, cmapi, cmapr, cmapg, cmapb, 0);
                for (l=0; l < snap->head->nsph; l++) {
                    if (((dump[l]-mindumpval)/(maxdumpval-mindumpval)) > 1.0)
                        plcol1(1.0);
                    else if (((dump[l]-mindumpval)/(maxdumpval-mindumpval)) < 0.0)
                        plcol1(0.0);
                    else
                        plcol1((dump[l]-mindumpval)/(maxdumpval-mindumpval));
                    plstring(1, &plotvars[j].points_xs[l], &plotvars[j].points_ys[l], "#(210)");
                }

                // colorbar
                cbar_valuex[0] = mindumpval; cbar_valuex[1] = maxdumpval;
                if (k<numdumps)
                    cbar_labels = &plotdumps[k];
                else
                    cbar_labels = &plotvecs[k-numdumps];
                plcolorbar(&cbar_width, &cbar_height, PL_COLORBAR_GRADIENT | PL_COLORBAR_LABEL_BOTTOM, 0,
                    0.005, 0.0, 0.05, 0.8, 0, 1, 1, 0.0, 0.0,
                    0, 0.0,
                    1, cbar_label_opts, cbar_labels,
                    1, cbar_axis_opts,
                    cbar_axis_ticks, cbar_axis_subticks,
                    cbar_n_values, (const double* const*) cbar_values);

                free(dump);
                plend();
            }
        }

        tipsyDestroy(snap);
        profileDestroy(pfile);
    }

    // make gifs
    float delay = (float)1000/(float)nout;
    for(i=0; i<numdumps; i++) {
        sprintf(command, "convert -layers optimize -delay %f ./test/%s/%s.*.png ./test/%s.animation.gif", delay, plotdumps[i], plotdumps[i], plotdumps[i]);
        printf("%s\n", command);
        system(command);
    }

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
