#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

/*
########   #######  ##     ## ##    ## ########   ######
##     ## ##     ## ##     ## ###   ## ##     ## ##    ##
##     ## ##     ## ##     ## ####  ## ##     ## ##
########  ##     ## ##     ## ## ## ## ##     ##  ######
##     ## ##     ## ##     ## ##  #### ##     ##       ##
##     ## ##     ## ##     ## ##   ### ##     ## ##    ##
########   #######   #######  ##    ## ########   ######
*/
void autoFindBounds(tipsy* tipsyIn){
    /* Finds the positions of point extrema in each axis, and sets the _min/_max
        values in the tipsy attributes to those values.

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy struct whose extrema
                                    will be found and written into tipsyIn->attr
    */
    // Indexing variables
    int i=0;
    // Set initial max and min based on the first particle
    if (tipsyIn->head->nsph != 0) {
        tipsyIn->attr->xmin = tipsyIn->gas[i].pos[AXIS_X];
        tipsyIn->attr->xmax = tipsyIn->gas[i].pos[AXIS_X];
        tipsyIn->attr->ymin = tipsyIn->gas[i].pos[AXIS_Y];
        tipsyIn->attr->ymax = tipsyIn->gas[i].pos[AXIS_Y];
        tipsyIn->attr->zmin = tipsyIn->gas[i].pos[AXIS_Z];
        tipsyIn->attr->zmax = tipsyIn->gas[i].pos[AXIS_Z];
    } else if (tipsyIn->head->ndark != 0) {
        tipsyIn->attr->xmin = tipsyIn->gas[i].pos[AXIS_X];
        tipsyIn->attr->xmax = tipsyIn->gas[i].pos[AXIS_X];
        tipsyIn->attr->ymin = tipsyIn->gas[i].pos[AXIS_Y];
        tipsyIn->attr->ymax = tipsyIn->gas[i].pos[AXIS_Y];
        tipsyIn->attr->zmin = tipsyIn->gas[i].pos[AXIS_Z];
        tipsyIn->attr->zmax = tipsyIn->gas[i].pos[AXIS_Z];
    } else if (tipsyIn->head->nstar != 0) {
        tipsyIn->attr->xmin = tipsyIn->star[i].pos[AXIS_X];
        tipsyIn->attr->xmax = tipsyIn->star[i].pos[AXIS_X];
        tipsyIn->attr->ymin = tipsyIn->star[i].pos[AXIS_Y];
        tipsyIn->attr->ymax = tipsyIn->star[i].pos[AXIS_Y];
        tipsyIn->attr->zmin = tipsyIn->star[i].pos[AXIS_Z];
        tipsyIn->attr->zmax = tipsyIn->star[i].pos[AXIS_Z];
    } else errorCase(ERR_NO_PARTICLES);
    // Find max min that changes
    for (i=1; i<tipsyIn->head->nsph; i++){
        if (tipsyIn->gas[i].pos[AXIS_X] < tipsyIn->attr->xmin)
            tipsyIn->attr->xmin = tipsyIn->gas[i].pos[AXIS_X];
        else if (tipsyIn->gas[i].pos[AXIS_X] > tipsyIn->attr->xmax)
            tipsyIn->attr->xmax = tipsyIn->gas[i].pos[AXIS_X];
        if (tipsyIn->gas[i].pos[AXIS_Y] < tipsyIn->attr->ymin)
            tipsyIn->attr->ymin = tipsyIn->gas[i].pos[AXIS_Y];
        else if (tipsyIn->gas[i].pos[AXIS_Y] > tipsyIn->attr->ymax)
            tipsyIn->attr->ymax = tipsyIn->gas[i].pos[AXIS_Y];
        if (tipsyIn->gas[i].pos[AXIS_Z] < tipsyIn->attr->zmin)
            tipsyIn->attr->zmin = tipsyIn->gas[i].pos[AXIS_Z];
        else if (tipsyIn->gas[i].pos[AXIS_Z] > tipsyIn->attr->zmax)
            tipsyIn->attr->zmax = tipsyIn->gas[i].pos[AXIS_Z];
    }
    for (i=1; i<tipsyIn->head->ndark; i++){
        if (tipsyIn->dark[i].pos[AXIS_X] < tipsyIn->attr->xmin)
            tipsyIn->attr->xmin = tipsyIn->dark[i].pos[AXIS_X];
        else if (tipsyIn->dark[i].pos[AXIS_X] > tipsyIn->attr->xmax)
            tipsyIn->attr->xmax = tipsyIn->dark[i].pos[AXIS_X];
        if (tipsyIn->dark[i].pos[AXIS_Y] < tipsyIn->attr->ymin)
            tipsyIn->attr->ymin = tipsyIn->dark[i].pos[AXIS_Y];
        else if (tipsyIn->dark[i].pos[AXIS_Y] > tipsyIn->attr->ymax)
            tipsyIn->attr->ymax = tipsyIn->dark[i].pos[AXIS_Y];
        if (tipsyIn->dark[i].pos[AXIS_Z] < tipsyIn->attr->zmin)
            tipsyIn->attr->zmin = tipsyIn->dark[i].pos[AXIS_Z];
        else if (tipsyIn->dark[i].pos[AXIS_Z] > tipsyIn->attr->zmax)
            tipsyIn->attr->zmax = tipsyIn->dark[i].pos[AXIS_Z];
    }
    for (i=1; i<tipsyIn->head->nstar; i++){
        if (tipsyIn->star[i].pos[AXIS_X] < tipsyIn->attr->xmin)
            tipsyIn->attr->xmin = tipsyIn->star[i].pos[AXIS_X];
        else if (tipsyIn->star[i].pos[AXIS_X] > tipsyIn->attr->xmax)
            tipsyIn->attr->xmax = tipsyIn->star[i].pos[AXIS_X];
        if (tipsyIn->star[i].pos[AXIS_Y] < tipsyIn->attr->ymin)
            tipsyIn->attr->ymin = tipsyIn->star[i].pos[AXIS_Y];
        else if (tipsyIn->star[i].pos[AXIS_Y] > tipsyIn->attr->ymax)
            tipsyIn->attr->ymax = tipsyIn->star[i].pos[AXIS_Y];
        if (tipsyIn->star[i].pos[AXIS_Z] < tipsyIn->attr->zmin)
            tipsyIn->attr->zmin = tipsyIn->star[i].pos[AXIS_Z];
        else if (tipsyIn->star[i].pos[AXIS_Z] > tipsyIn->attr->zmax)
            tipsyIn->attr->zmax = tipsyIn->star[i].pos[AXIS_Z];
    }
}

/*
########  ######## ########    ###    ##     ## ##       ########  ######
##     ## ##       ##         ## ##   ##     ## ##          ##    ##    ##
##     ## ##       ##        ##   ##  ##     ## ##          ##    ##
##     ## ######   ######   ##     ## ##     ## ##          ##     ######
##     ## ##       ##       ######### ##     ## ##          ##          ##
##     ## ##       ##       ##     ## ##     ## ##          ##    ##    ##
########  ######## ##       ##     ##  #######  ########    ##     ######
*/
void tipsySetDefaults(tipsy* tipsyIn){
    /* Fills in the input struct with some default values in case the struct
        maker couldn't be bothered. Really, I just put this here in case I need
        it again after removing it from createTipsy().

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy to set default values for
    */
    // Indexing variables
    int i;

    // Set gas values
    for (i=0; i < tipsyIn->head->nsph; i++){
        tipsyIn->gas[i].mass = 1.0;
        tipsyIn->gas[i].pos[0] = VAL_NaN;
        tipsyIn->gas[i].pos[1] = VAL_NaN;
        tipsyIn->gas[i].pos[2] = VAL_NaN;
        tipsyIn->gas[i].vel[0] = 0.0;
        tipsyIn->gas[i].vel[1] = 0.0;
        tipsyIn->gas[i].vel[2] = 0.0;
        tipsyIn->gas[i].rho = 1.0;
        tipsyIn->gas[i].temp = 1.0;
        tipsyIn->gas[i].eps = 1.0;
        tipsyIn->gas[i].metals = 0.0;
        tipsyIn->gas[i].phi = 1.0;
    }
    // Set dark values
    for (i=0; i < tipsyIn->head->ndark; i++){
        tipsyIn->dark[i].mass = 1.0;
        tipsyIn->dark[i].pos[AXIS_X] = VAL_NaN;
        tipsyIn->dark[i].pos[AXIS_Y] = VAL_NaN;
        tipsyIn->dark[i].pos[AXIS_Z] = VAL_NaN;
        tipsyIn->dark[i].vel[AXIS_X] = 0.0;
        tipsyIn->dark[i].vel[AXIS_Y] = 0.0;
        tipsyIn->dark[i].vel[AXIS_Z] = 0.0;
        tipsyIn->dark[i].eps = 1.0;
        tipsyIn->dark[i].phi = 1.0;
    }
    // Set star values
    for (i=0; i < tipsyIn->head->nstar; i++){
        tipsyIn->star[i].mass = 1.0;
        tipsyIn->star[i].pos[AXIS_X] = VAL_NaN;
        tipsyIn->star[i].pos[AXIS_Y] = VAL_NaN;
        tipsyIn->star[i].pos[AXIS_Z] = VAL_NaN;
        tipsyIn->star[i].vel[AXIS_X] = 0.0;
        tipsyIn->star[i].vel[AXIS_Y] = 0.0;
        tipsyIn->star[i].vel[AXIS_Z] = 0.0;
        tipsyIn->star[i].metals = 0.0;
        tipsyIn->star[i].tform = 0.0;
        tipsyIn->star[i].eps = 1.0;
        tipsyIn->star[i].phi = 1.0;
    }
}
