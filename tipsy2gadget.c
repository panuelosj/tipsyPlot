#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "tipsyPlot.h"

int main(int argc, char* argv[]){
    if (argc == 1){
        system("cat nagato");
        printf("Error: No input file provided\n\tPlease enter an input tipsy to convert\n");
        printf("Usage: tipsy2gadget [input tipsy] [xmin] [xmax] [ymin] [ymax] [zmin] [zmax]\n");
        printf("\tIf no positional extrema are provided, they will be automatically set to the particle extrema found in the input tipsy\n");
        exit(-1);
    }
    if (argc != 2 && argc != 3 && argc != 8 && argc != 9){
        system("cat nagato");
        printf("Error: Incorrect number of parameters\n");
        printf("Usage: tipsy2gadget [input tipsy filename] [output GADGET filename]\n\t\t [xmin] [xmax] [ymin] [ymax] [zmin] [zmax]\n");
        printf("\tPlease ensure you are providing all position extrema, or none at all\n");
        printf("\tIf no positional extrema are provided, they will be automatically set to the particle extrema found in the input tipsy\n");
        printf("\tIf no filename is provided, it will be automatically set to [input tipsy].gadget");
        exit(-1);
    }

    // Read in the input tipsy
    printf("Reading: %s\n", argv[1]);
    tipsy* tipsyIn = readTipsyStd(argv[1]);
    // Set attributes if they are provided, else readTipsyStd() would already have automatically set the found extrema
    if (argc == 8){
        tipsyIn->attr->xmin = (float)atof(argv[2]); tipsyIn->attr->xmax = (float)atof(argv[3]);
        tipsyIn->attr->ymin = (float)atof(argv[4]); tipsyIn->attr->ymax = (float)atof(argv[5]);
        tipsyIn->attr->zmin = (float)atof(argv[6]); tipsyIn->attr->zmax = (float)atof(argv[7]);
    }
    else if (argc == 9){
        tipsyIn->attr->xmin = (float)atof(argv[3]); tipsyIn->attr->xmax = (float)atof(argv[4]);
        tipsyIn->attr->ymin = (float)atof(argv[5]); tipsyIn->attr->ymax = (float)atof(argv[6]);
        tipsyIn->attr->zmin = (float)atof(argv[7]); tipsyIn->attr->zmax = (float)atof(argv[8]);
    }

    printf("Input ");
    printHeader(tipsyIn->head);
    printAttr(tipsyIn->attr);
    printf("=================================================\n");

    if (argc == 2 || argc == 8){
        char defaultfileout[100];
        sprintf(defaultfileout, "%s.gadget", argv[1]);
        writeGadgetFromTipsy(defaultfileout, tipsyIn);
    }
    else if(argc == 3 || argc == 9){
        writeGadgetFromTipsy(argv[2], tipsyIn);
    }

    // Cleanup
    tipsyDestroy(tipsyIn);
}
