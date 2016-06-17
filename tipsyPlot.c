#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "tipsyPlot.h"

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

    const float minx = -7, maxx = 7;

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
    for(i=0; i<nout; i++){
        int filetime = (i+1)*interval;                                          // calculate current filename time
        sprintf(simname, "%s.%05d", genericfilename, filetime);                     // generate filename
        tipsy* snap = readTipsyStd(simname);


    }
}
