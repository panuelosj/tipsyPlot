#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

/*
 ######  ######## ##    ## ######## ######## ########
##    ## ##       ###   ##    ##    ##       ##     ##
##       ##       ####  ##    ##    ##       ##     ##
##       ######   ## ## ##    ##    ######   ########
##       ##       ##  ####    ##    ##       ##   ##
##    ## ##       ##   ###    ##    ##       ##    ##
 ######  ######## ##    ##    ##    ######## ##     ##
*/
void tipsyCenter(tipsy* tipsyIn){
    /* Automatically calculates the center of the tipsy simulation box according
        to _min/_max and shifts the simulation box to have it centered at
        (0,0,0).

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy struct to center
    */

    const float xCenter = (tipsyIn->attr->xmax + tipsyIn->attr->xmin)/2.0;
    const float yCenter = (tipsyIn->attr->ymax + tipsyIn->attr->ymin)/2.0;
    const float zCenter = (tipsyIn->attr->zmax + tipsyIn->attr->zmin)/2.0;
    tipsyTranslate(tipsyIn, -1.0*xCenter, -1.0*yCenter, -1.0*zCenter);
}

/*
 ######  ########   #######  ########
##    ## ##     ## ##     ## ##     ##
##       ##     ## ##     ## ##     ##
##       ########  ##     ## ########
##       ##   ##   ##     ## ##
##    ## ##    ##  ##     ## ##
 ######  ##     ##  #######  ##
*/
void tipsyCrop(tipsy* tipsyIn, crop op){
    /* Crops the simulation to any parameter according to crop op. The function
        uses crop op which returns a binary value telling it to keep or throw
        out the current particle. If crop op returns nonzero, the current
        particle is copied to the last position recorded by nkept_, which should
        either be an unkept particle, or the current particle (only if all
        previous particles have been kept). Once all the kept particles have
        been copied to the beginning of the tipsy struct, the struct is resized
        to include only those particles, effectively deleting the end of the
        struct containing irrelevant particles (only has unkept particles and
        kept particles copied to the beginning of the struct).

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy struct to crop
            crop op             - operation returning binary value to determine
                                    if a particle is to be kept or not
    */

    int i;
    int nkeptgas = 0, nkeptdark = 0, nkeptstar = 0;

    for (i=0; i<tipsyIn->attr->nloadedsph; i++){
        if (op(&tipsyIn->gas[i], TYPE_GAS) != 0) {
            pFlopGas(&tipsyIn->gas[nkeptgas], &tipsyIn->gas[i], NULL, flopCopy);
            nkeptgas++;
        }
    }
    for (i=0; i<tipsyIn->attr->nloadeddark; i++){
        if (op(&tipsyIn->dark[i], TYPE_DARK) != 0){
            pFlopDark(&tipsyIn->dark[nkeptdark], &tipsyIn->dark[i], NULL, flopCopy);
            nkeptdark++;
        }
    }
    for (i=0; i<tipsyIn->attr->nloadedstar; i++){
        if (op(&tipsyIn->star[i], TYPE_STAR) != 0){
            pFlopStar(&tipsyIn->star[nkeptstar], &tipsyIn->star[i], NULL, flopCopy);
        }
    }
    tipsyExtend(tipsyIn, nkeptgas, nkeptdark, nkeptstar);
}


/*
 ######   ######     ###    ##       ########
##    ## ##    ##   ## ##   ##       ##
##       ##        ##   ##  ##       ##
 ######  ##       ##     ## ##       ######
      ## ##       ######### ##       ##
##    ## ##    ## ##     ## ##       ##
 ######   ######  ##     ## ######## ########
*/
void tipsyScaleShrink(tipsy* tipsyIn, const int xShrink, const int yShrink, const int zShrink){
    /* Shrinks the tipsy object in each dimension by the given scaling factors
        by dividing by that value. For consistency with other functions, only
        integer values should be given for the shrinking factors, but are
        casted to floats to give a float result when dividing.

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy struct to be shrunk
            const int xShrink   - x dimension compression factor
            const int yShrink   - y dimension compression factor
            const int zShrink   - z dimension compression factor
        Return:
            void
    */

    int i;                                                                      // indexing variables
    const float xShrinkF = (float) xShrink;                                     // float casts
    const float yShrinkF = (float) yShrink;
    const float zShrinkF = (float) zShrink;
    const float totalShrinkF = xShrinkF * yShrinkF * zShrinkF;

    // Shrink each element's coordinates by the scaling factor given, update rho
    for (i=0; i < tipsyIn->head->nsph; i++){
        tipsyIn->gas[i].pos[AXIS_X] /= xShrinkF;
        tipsyIn->gas[i].pos[AXIS_Y] /= yShrinkF;
        tipsyIn->gas[i].pos[AXIS_Z] /= zShrinkF;
        tipsyIn->gas[i].rho *= totalShrinkF;
    }
    for (i=0; i < tipsyIn->head->ndark; i++){
        tipsyIn->dark[i].pos[AXIS_X] /= xShrinkF;
        tipsyIn->dark[i].pos[AXIS_Y] /= yShrinkF;
        tipsyIn->dark[i].pos[AXIS_Z] /= zShrinkF;
    }
    for (i=0; i < tipsyIn->head->nstar; i++){
        tipsyIn->star[i].pos[AXIS_X] /= xShrinkF;
        tipsyIn->star[i].pos[AXIS_Y] /= yShrinkF;
        tipsyIn->star[i].pos[AXIS_Z] /= zShrinkF;
    }

    // Shrink each boundary by the scaling factor given
    tipsyIn->attr->xmin /= xShrinkF; tipsyIn->attr->xmax /= xShrinkF;
    tipsyIn->attr->ymin /= yShrinkF; tipsyIn->attr->ymax /= yShrinkF;
    tipsyIn->attr->zmin /= zShrinkF; tipsyIn->attr->zmax /= zShrinkF;
}

void tipsyScaleExpand(tipsy* tipsyIn, const float xExpand, const float yExpand, const float zExpand){
    /* Expands the tipsy object in each dimension by the given scaling factors
        by multiplying by that value.

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy struct to be shrunk
            const float xExpand - x dimension expansion factor
            const float yExpand - y dimension expansion factor
            const float zExpand - z dimension expansion factor
        Return:
            void
    */

    int i;                                                                      // indexing variables
    const float xExpandF = (float) xExpand;                                     // float casts
    const float yExpandF = (float) yExpand;
    const float zExpandF = (float) zExpand;
    const float totalExpandF = xExpandF * yExpandF * zExpandF;

    // Expand each element's coordinates by the scaling factor given, update rho
    for (i=0; i < tipsyIn->head->nsph; i++){
        tipsyIn->gas[i].pos[AXIS_X] *= xExpandF;
        tipsyIn->gas[i].pos[AXIS_Y] *= yExpandF;
        tipsyIn->gas[i].pos[AXIS_Z] *= zExpandF;
        tipsyIn->gas[i].rho /= totalExpandF;
    }
    for (i=0; i < tipsyIn->head->ndark; i++){
        tipsyIn->dark[i].pos[AXIS_X] *= xExpandF;
        tipsyIn->dark[i].pos[AXIS_Y] *= yExpandF;
        tipsyIn->dark[i].pos[AXIS_Z] *= zExpandF;
    }
    for (i=0; i < tipsyIn->head->nstar; i++){
        tipsyIn->star[i].pos[AXIS_X] *= xExpandF;
        tipsyIn->star[i].pos[AXIS_Y] *= yExpandF;
        tipsyIn->star[i].pos[AXIS_Z] *= zExpandF;
    }

    // Expand each boundary by the scaling factor given
    tipsyIn->attr->xmin *= xExpandF; tipsyIn->attr->xmax *= xExpandF;
    tipsyIn->attr->ymin *= yExpandF; tipsyIn->attr->ymax *= yExpandF;
    tipsyIn->attr->zmin *= zExpandF; tipsyIn->attr->zmax *= zExpandF;
}

/*
######## ########  ######   ######  ######## ##          ###    ######## ########
   ##    ##       ##    ## ##    ## ##       ##         ## ##      ##    ##
   ##    ##       ##       ##       ##       ##        ##   ##     ##    ##
   ##    ######    ######   ######  ######   ##       ##     ##    ##    ######
   ##    ##             ##       ## ##       ##       #########    ##    ##
   ##    ##       ##    ## ##    ## ##       ##       ##     ##    ##    ##
   ##    ########  ######   ######  ######## ######## ##     ##    ##    ########
*/
void tipsyTesselate(tipsy* tipsyIn, const int xTile, const int yTile, const int zTile){
    /* Tesselates (tile) the tipsy object across each axis, copying over the
        positive axis directions. This will automatically extend (realloc) the
        input tipsy struct to fit the larger tiled tipsy.
       The code currently tesselates by first duplicating the existing particles
        using memcpy, then shifting the relevant copied indices. This is
        repeated using the newly extended tipsy in each direction. Another
        method could be to use a clone of the original as a tiling unit,
        shifting that and copying the shifted clone unit to the original tipsyIn
        thereby tesselating it. This would use more memory for the cloned unit,
        which would also need to be freed later, and would likely not be any
        faster, since the current method copies the finished product of the last
        axis copy (ie, the lengthened x is copied when copying over the y-axis,
        and the xy sheet is copied over the z-axis)

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy struct to be tiled
            const int xTile     - number of units of the original copied over x axis
            const int yTile     - number of units of the original copied over y axis
            const int zTile     - number of units of the original copied over z axis

        ToDo:
            add check to see if tipsyIn is full (ie n_=nloaded_);
                currently assumes tipsyIn is full and uses n_ values
    */

    int i, j;                                                                   // indexing variables
    const int nTile = xTile*yTile*zTile;                                        // total number of tiling units in the final product
    const float xWidth = tipsyIn->attr->xmax - tipsyIn->attr->xmin;             // total width of the tipsy sim box
    const float yWidth = tipsyIn->attr->ymax - tipsyIn->attr->ymin;
    const float zWidth = tipsyIn->attr->zmax - tipsyIn->attr->zmin;

    // Allocate enough memory for the particles to be copied
    tipsyExtend(tipsyIn, tipsyIn->head->nsph*nTile, tipsyIn->head->ndark*nTile, tipsyIn->head->nstar*nTile);

    // Duplicate over x axis
    for (j=1; j < xTile; j++){
        // For each tile repeat stated in xTile, duplicate the memory into the right location
        memcpy(&tipsyIn->gas[j*tipsyIn->attr->nloadedsph], &tipsyIn->gas[0], tipsyIn->attr->nloadedsph*sizeof(gas_particle));
        memcpy(&tipsyIn->dark[j*tipsyIn->attr->nloadeddark], &tipsyIn->dark[0], tipsyIn->attr->nloadeddark*sizeof(dark_particle));
        memcpy(&tipsyIn->star[j*tipsyIn->attr->nloadedstar], &tipsyIn->star[0], tipsyIn->attr->nloadedstar*sizeof(star_particle));
        // Shift all particles over based on the integer multiple of tile width
        for (i=j*tipsyIn->attr->nloadedsph; i < (j+1)*tipsyIn->attr->nloadedsph; i++)
            tipsyIn->gas[i].pos[AXIS_X] += (float)j*xWidth;
        for (i=j*tipsyIn->attr->nloadeddark; i < (j+1)*tipsyIn->attr->nloadeddark; i++)
            tipsyIn->dark[i].pos[AXIS_X] += (float)j*xWidth;
        for (i=j*tipsyIn->attr->nloadedstar; i < (j+1)*tipsyIn->attr->nloadedstar; i++)
            tipsyIn->star[i].pos[AXIS_X] += (float)j*xWidth;
    }
    // Update the number of loaded particles
    tipsyIn->attr->nloadedsph *= xTile;
    tipsyIn->attr->nloadeddark *= xTile;
    tipsyIn->attr->nloadedstar *= xTile;

    // Duplicate over y axis
    for (j=1; j < yTile; j++){
        // For each tile repeat stated in xTile, duplicate the memory into the right location
        memcpy(&tipsyIn->gas[j*tipsyIn->attr->nloadedsph], &tipsyIn->gas[0], tipsyIn->attr->nloadedsph*sizeof(gas_particle));
        memcpy(&tipsyIn->dark[j*tipsyIn->attr->nloadeddark], &tipsyIn->dark[0], tipsyIn->attr->nloadeddark*sizeof(dark_particle));
        memcpy(&tipsyIn->star[j*tipsyIn->attr->nloadedstar], &tipsyIn->star[0], tipsyIn->attr->nloadedstar*sizeof(star_particle));
        // Shift all particles over based on the integer multiple of tile width
        for (i=j*tipsyIn->attr->nloadedsph; i < (j+1)*tipsyIn->attr->nloadedsph; i++)
            tipsyIn->gas[i].pos[AXIS_Y] += (float)j*yWidth;
        for (i=j*tipsyIn->attr->nloadeddark; i < (j+1)*tipsyIn->attr->nloadeddark; i++)
            tipsyIn->dark[i].pos[AXIS_Y] += (float)j*yWidth;
        for (i=j*tipsyIn->attr->nloadedstar; i < (j+1)*tipsyIn->attr->nloadedstar; i++)
            tipsyIn->star[i].pos[AXIS_Y] += (float)j*yWidth;
    }
    // Update the number of loaded particles
    tipsyIn->attr->nloadedsph *= yTile;
    tipsyIn->attr->nloadeddark *= yTile;
    tipsyIn->attr->nloadedstar *= yTile;

    // Duplicate over z axis
    for (j=1; j < zTile; j++){
        // For each tile repeat stated in xTile, duplicate the memory into the right location
        memcpy(&tipsyIn->gas[j*tipsyIn->attr->nloadedsph], &tipsyIn->gas[0], tipsyIn->attr->nloadedsph*sizeof(gas_particle));
        memcpy(&tipsyIn->dark[j*tipsyIn->attr->nloadeddark], &tipsyIn->dark[0], tipsyIn->attr->nloadeddark*sizeof(dark_particle));
        memcpy(&tipsyIn->star[j*tipsyIn->attr->nloadedstar], &tipsyIn->star[0], tipsyIn->attr->nloadedstar*sizeof(star_particle));
        // Shift all particles over based on the integer multiple of tile width
        for (i=j*tipsyIn->attr->nloadedsph; i < (j+1)*tipsyIn->attr->nloadedsph; i++)
            tipsyIn->gas[i].pos[AXIS_Z] += (float)j*zWidth;
        for (i=j*tipsyIn->attr->nloadeddark; i < (j+1)*tipsyIn->attr->nloadeddark; i++)
            tipsyIn->dark[i].pos[AXIS_Z] += (float)j*zWidth;
        for (i=j*tipsyIn->attr->nloadedstar; i < (j+1)*tipsyIn->attr->nloadedstar; i++)
            tipsyIn->star[i].pos[AXIS_Z] += (float)j*zWidth;
    }
    // Update the number of loaded particles
    tipsyIn->attr->nloadedsph *= zTile;
    tipsyIn->attr->nloadeddark *= zTile;
    tipsyIn->attr->nloadedstar *= zTile;

    // Update boundaries
    tipsyIn->attr->xmax += (xTile-1)*xWidth;
    tipsyIn->attr->ymax += (yTile-1)*yWidth;
    tipsyIn->attr->zmax += (zTile-1)*zWidth;
}

/*
######## ########     ###    ##    ##  ######  ##          ###    ######## ########
   ##    ##     ##   ## ##   ###   ## ##    ## ##         ## ##      ##    ##
   ##    ##     ##  ##   ##  ####  ## ##       ##        ##   ##     ##    ##
   ##    ########  ##     ## ## ## ##  ######  ##       ##     ##    ##    ######
   ##    ##   ##   ######### ##  ####       ## ##       #########    ##    ##
   ##    ##    ##  ##     ## ##   ### ##    ## ##       ##     ##    ##    ##
   ##    ##     ## ##     ## ##    ##  ######  ######## ##     ##    ##    ########
*/
void tipsyTranslate(tipsy* tipsyIn, const float xShift, const float yShift, const float zShift){
    /* Translates (shifts) the tipsy object across each axis, moving the
        particle position values according to the input parameters. A positive
        _Shift value will move the particles in the positive axis direction and
        a negative value will move them in the negative axis direction.

        Parameters:
            tipsy* tipsyIn      - pointer to the tipsy struct to be translated
            const float xShift  - amount to translate the particles by in the x axis
            const float yShift  - amount to translate the particles by in the y axis
            const float zShift  - amount to translate the particles by in the z axis

        ToDo:
    */

    int i;                                                                      // indexing variables
    // Shift each particle type
    for (i=0; i < tipsyIn->attr->nloadedsph; i++){
        tipsyIn->gas[i].pos[AXIS_X] += xShift;
        tipsyIn->gas[i].pos[AXIS_Y] += yShift;
        tipsyIn->gas[i].pos[AXIS_Z] += zShift;
    }
    for (i=0; i < tipsyIn->attr->nloadeddark; i++){
        tipsyIn->dark[i].pos[AXIS_X] += xShift;
        tipsyIn->dark[i].pos[AXIS_Y] += yShift;
        tipsyIn->dark[i].pos[AXIS_Z] += zShift;
    }
    for (i=0; i < tipsyIn->attr->nloadedstar; i++){
        tipsyIn->star[i].pos[AXIS_X] += xShift;
        tipsyIn->star[i].pos[AXIS_Y] += yShift;
        tipsyIn->star[i].pos[AXIS_Z] += zShift;
    }
    // Update the box boundaries
    tipsyIn->attr->xmin += xShift; tipsyIn->attr->xmax += xShift;
    tipsyIn->attr->ymin += yShift; tipsyIn->attr->ymax += yShift;
    tipsyIn->attr->zmin += zShift; tipsyIn->attr->zmax += zShift;
}
