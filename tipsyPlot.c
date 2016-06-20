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
    ##     ##    ###    #### ##    ##
    ###   ###   ## ##    ##  ###   ##
    #### ####  ##   ##   ##  ####  ##
    ## ### ## ##     ##  ##  ## ## ##
    ##     ## #########  ##  ##  ####
    ##     ## ##     ##  ##  ##   ###
    ##     ## ##     ## #### ##    ##
    */

    // Find min-max bounds
    for(i=0; i<nout; i++){
        int filetime = (i+1)*interval;                                          // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                 // generate filename
        tipsy* snap = readTipsyStd(simname);
        profile* pfile = profileCreate(snap, nbins, minx, maxx, xpos);

        if (i==0){

        } else{

        }
    }
}


// calc_var

// calc_bin
float xpos(tipsy* tipsyIn, int type, int particle){
    switch (type){
        case TYPE_GAS:
            return tipsyIn->gas[particle].pos[AXIS_X];
        case TYPE_DARK:
            return tipsyIn->dark[particle].pos[AXIS_X];
        case TYPE_STAR:
            return tipsyIn->star[particle].pos[AXIS_X];
    }
}
