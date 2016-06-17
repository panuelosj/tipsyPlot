#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

/*
########  ########    ###    ########
##     ## ##         ## ##   ##     ##
##     ## ##        ##   ##  ##     ##
########  ######   ##     ## ##     ##
##   ##   ##       ######### ##     ##
##    ##  ##       ##     ## ##     ##
##     ## ######## ##     ## ########
*/
tipsy* readTipsyStd(const char filename[]){
    /* Reads in a file in tipsy standard format and writes it into a new tipsy
        struct. createTipsy() is not called since size of particle structs
        are unknown until header is read. Header bit could be read first, then
        the struct created using createTipsy(), but readTipsyStd() and
        createTipsy() are pretty much completed and further changes are
        unlikely, so consistency in changes isn't really a problem. Swaps the
        endianness of the input values when writing to the tipsy struct since
        most modern machines are little endian and tipsy standard is a big
        endian format, but I should check the machine endianness for better
        compatibility.
        (also idk if there's weird stuff with floating point endianness but
        currently doesn't look like there is after checking the values with
        hexdump)

        Parameters:
            const char filename[]   - filename to be opened
        Return:
            pointer to newly created tipsy struct containing the values
                read in from the file

        ToDo:
            Include an endian checker to automatically swap the endianness if
                it is not yet in (big endian) - possibly a compilertime option
                and make two different source files to compile?
            Include checking if the file is in the correct format??
    */

    // Indexing variables
    int i;

    FILE *fp = fopen(filename, "r");                                            // Open file pointer
    // Create object
    tipsy* newTipsy = (tipsy*)malloc(sizeof(tipsy));

    // Create and read in header
    newTipsy->head = (header*)malloc(sizeof(header));                                  // Allocate space
    fread(&newTipsy->head->simtime, sizeof(double), 1, fp);                   // Read in binary
    fread(&newTipsy->head->nbodies, sizeof(int), 6, fp);
    swapEndianBatch(newTipsy, TYPE_HEADER, 0);                                  // Swap endianness
    // Create and read in gas particles
    if (newTipsy->head->nsph != 0){                                           // Check if particles exist
        newTipsy->gas = (gas_particle*)malloc(newTipsy->head->nsph*sizeof(gas_particle));    // Allocate space
        for (i=0; i < newTipsy->head->nsph; i++){
            fread(&newTipsy->gas[i].mass, sizeof(float), 12, fp);               // Read in binary
            swapEndianBatch(newTipsy, TYPE_GAS, i);                             // Swap endianness
        }
    } else newTipsy->gas = NULL;
    // Create and read in dark matter particles
    if (newTipsy->head->ndark != 0){                                          // Check if particles exist
        newTipsy->dark = (dark_particle*)malloc(newTipsy->head->nsph*sizeof(dark_particle));  // Allocate space
        for (i=0; i < newTipsy->head->ndark; i++){
            fread(&newTipsy->dark[i].mass, sizeof(float), 9, fp);               // Read in binary
            swapEndianBatch(newTipsy, TYPE_DARK, i);                            // Swap endianness
        }
    } else newTipsy->dark = NULL;
    // Create and read in star particles
    if (newTipsy->head->nstar != 0){                                          // Check if particles exist
        newTipsy->star = (star_particle*)malloc(newTipsy->head->nsph*sizeof(star_particle));  // Allocate space
        for (i=0; i < newTipsy->head->nstar; i++){
            fread(&newTipsy->star[i].mass, sizeof(float), 11, fp);              // Read in binary
            swapEndianBatch(newTipsy, TYPE_STAR, i);                            // Swap endianness
        }
    } else newTipsy->star = NULL;

    // Create and set attributes
    newTipsy->attr = (attributes*)malloc(sizeof(attributes));
    newTipsy->attr->nloadedsph = newTipsy->head->nsph;
    newTipsy->attr->nloadeddark = newTipsy->head->ndark;
    newTipsy->attr->nloadedstar = newTipsy->head->nstar;
    autoFindBounds(newTipsy);

    // Cleanup
    fclose(fp);
    // Output a pointer to the new tipsy object
    return newTipsy;
}

/*
##      ## ########  #### ######## ########
##  ##  ## ##     ##  ##     ##    ##
##  ##  ## ##     ##  ##     ##    ##
##  ##  ## ########   ##     ##    ######
##  ##  ## ##   ##    ##     ##    ##
##  ##  ## ##    ##   ##     ##    ##
 ###  ###  ##     ## ####    ##    ########
*/
int writeTipsyStd(const char filename[], tipsy* tipsyOut){
    /* Writes an input tipsy struct into a file with the tipsy standard format.
        Function is optimized for minimizing space usage, keeping the endian-
        swapped values in the same tipsy struct as inputted, and swapping the
        endianness again after the values are copied to file.

        Parameters:
            const char filename[]   - filename to write to
            tipsy* tipsyOut         - tipsy struct to write to file
        Return:
            Confirmation of write or error code?

        ToDo:
            Make return confirmation value or error code.
    */
    int i;

    FILE *fp = fopen(filename, "w");
    // Write header
    swapEndianBatch(tipsyOut, TYPE_HEADER, 0);                                  // Swap to .std endianness
    fwrite(&tipsyOut->head->simtime, sizeof(double), 1, fp);                  // Write to binary
    fwrite(&tipsyOut->head->nbodies, sizeof(int), 6, fp);
    swapEndianBatch(tipsyOut, TYPE_HEADER, 0);                                  // Return to working in machine endianness
    // Write gas particles
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            swapEndianBatch(tipsyOut, TYPE_GAS, i);                             // Swap to .std endianness
            fwrite(&tipsyOut->gas[i].mass, sizeof(float), 12, fp);              // Write to binary
            swapEndianBatch(tipsyOut, TYPE_GAS, i);                             // Return to working in machine endianness
        }
    }
    // Write dark matter particles
    if (tipsyOut->head->ndark != 0){
        for (i=0; i < tipsyOut->head->ndark; i++){
            swapEndianBatch(tipsyOut, TYPE_DARK, i);                             // Swap to .std endianness
            fwrite(&tipsyOut->dark[i].mass, sizeof(float), 9, fp);              // Write to binary
            swapEndianBatch(tipsyOut, TYPE_DARK, i);                             // Return to working in machine endianness
        }
    }
    // Write gas particles
    if (tipsyOut->head->nstar != 0){
        for (i=0; i < tipsyOut->head->nstar; i++){
            swapEndianBatch(tipsyOut, TYPE_STAR, i);                             // Swap to .std endianness
            fwrite(&tipsyOut->star[i].mass, sizeof(float), 12, fp);              // Write to binary
            swapEndianBatch(tipsyOut, TYPE_STAR, i);                             // Return to working in machine endianness
        }
    }

    // Cleanup
    fclose(fp);
}
