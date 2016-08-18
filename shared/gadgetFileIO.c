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
tipsy* readGadgetToTipsy(const char filename[]){
    /* Reads in a file in GADGET format and writes it into a new tipsy struct.
        createTipsy() is not called since size of particle structs are unknown
        until header is read.

        Parameters:
            const char filename[]   - filename to be opened
        Return:
            pointer to newly created tipsy struct containing the values
                read in from the file

        ToDo:
            - ignore mass block if mass table is nonzero
    */

    // Indexing variables
    int i;
    int blocksize;
    int* pID;

    FILE *fp = fopen(filename, "r");                                            // Open the file pointer
    // Create object
    tipsy* newTipsy = (tipsy*)malloc(sizeof(tipsy));

    // Create and read in header
    newTipsy->head = (header*)malloc(sizeof(header));
    fread(&blocksize, sizeof(int), 1, fp);                                      // Read in block size
    if (blocksize != 256){
        printf("Error: Header block size [%d bytes] is incorrect\n", blocksize);
        exit(-1);
    }
    fread(&newTipsy->head->nsph, sizeof(int), 1, fp);                           // Number of particles (inside this file)
    printf("\tnsph = %d\n", newTipsy->head->nsph);
    fread(&newTipsy->head->ndark, sizeof(int), 1, fp);
    printf("\tndark = %d\n", newTipsy->head->ndark);
    fseek(fp, sizeof(int)*2, SEEK_CUR);
    fread(&newTipsy->head->nstar, sizeof(int), 1, fp);
    printf("\tnstar = %d\n", newTipsy->head->nstar);
    fseek(fp, sizeof(int), SEEK_CUR);

    fseek(fp, sizeof(double)*6, SEEK_CUR);                                      // Mass table (ignore this at the moment since masses will be read from file)
    fread(&newTipsy->head->simtime, sizeof(double), 1, fp);                     // Time
    printf("\tsimtime = %f\n", newTipsy->head->simtime);
    fseek(fp, sizeof(double), SEEK_CUR);                                        // Redshift, ignore this atm since I'm not dealing with cosmological stuff
    fseek(fp, sizeof(int)*10, SEEK_CUR);                                        // Ignore all this
    fseek(fp, sizeof(double)*4, SEEK_CUR);                                      // Boxsize, ignore this since extrema will automatically be found from particle positions
    fseek(fp, sizeof(int)*9, SEEK_CUR);                                         // Ignore
    fseek(fp, 256-(sizeof(int)*25+sizeof(double)*12), SEEK_CUR);                // Skip to the end of the 256 byte long header

    fread(&blocksize, sizeof(int), 1, fp);                                      // Read in block size, check if it is consistent
    if (blocksize != 256){
        printf("Error: Header ending block size [%d bytes] is inconsistent\n", blocksize);
        exit(-1);
    }
    newTipsy->head->nbodies = newTipsy->head->nsph + newTipsy->head->ndark + newTipsy->head->nstar;

    // Create particles
    if (newTipsy->head->nbodies == 0){
        printf("Error: No particles found\n");
        exit(-1);
    }
    pID = (int*)malloc(newTipsy->head->nbodies*sizeof(int));
    if (newTipsy->head->nsph != 0)
        newTipsy->gas = (gas_particle*)malloc(newTipsy->head->nsph*sizeof(gas_particle));
    if (newTipsy->head->ndark != 0)
        newTipsy->dark = (dark_particle*)malloc(newTipsy->head->ndark*sizeof(dark_particle));
    if (newTipsy->head->nstar != 0)
        newTipsy->star = (star_particle*)malloc(newTipsy->head->nstar*sizeof(star_particle));

    // Skip to Particle IDs to for writing the particles in the correct order
    fseek(fp, newTipsy->head->nbodies*sizeof(float)*6 + sizeof(int)*4, SEEK_CUR);
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(int)){
        printf("Error: ParticleID block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(int)), newTipsy->head->nbodies);
        exit(-1);
    }
    for (i=0; i<newTipsy->head->nbodies; i++){
        fread(&pID[i], sizeof(int), 1, fp);
        pID[i]--;
        if (pID[i] >= newTipsy->head->nbodies){
            printf("Warning: pID [%d] is greater than the number of particles [%d] allocated\n", pID[i], newTipsy->head->nbodies);
            printf("\tSetting pID to 0\n");
            pID[i] = 0;
        }
    }
    //fread(pID, sizeof(int), newTipsy->head->nbodies, fp);                      // copy particle IDs
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(int)){
        printf("Error: ParticleID ending block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(int)), newTipsy->head->nbodies);
        exit(-1);
    }
    fseek(fp, 256 + sizeof(int)*2, SEEK_SET);                                   // return to the end of the header

    // Read in Positions
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(float)*3){
        printf("Error: Position block size [%lu particles] is inconsistent with the the header [%d particles]\n", blocksize/(sizeof(float)*3), newTipsy->head->nbodies);
        exit(-1);
    }
    if (newTipsy->head->nsph != 0)
        for (i=0; i < newTipsy->head->nsph; i++)
            fread(&newTipsy->gas[pID[i]].pos, sizeof(float), 3, fp);
    if (newTipsy->head->ndark != 0)
        for (i=0; i < newTipsy->head->ndark; i++)
            fread(&newTipsy->dark[pID[i + newTipsy->head->nsph]].pos, sizeof(float), 3, fp);
    if (newTipsy->head->nstar != 0)
        for (i=0; i < newTipsy->head->nstar; i++)
            fread(&newTipsy->star[pID[i  + newTipsy->head->nsph + newTipsy->head->ndark]].pos, sizeof(float), 3, fp);
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(float)*3){
        printf("Error: Position ending block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(float)*3), newTipsy->head->nbodies);
        exit(-1);
    }

    // Read in Velocities
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(float)*3){
        printf("Error: Velocity block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(float)*3), newTipsy->head->nbodies);
        exit(-1);
    }
    if (newTipsy->head->nsph != 0)
        for (i=0; i < newTipsy->head->nsph; i++)
            fread(&newTipsy->gas[pID[i]].vel, sizeof(float), 3, fp);
    if (newTipsy->head->ndark != 0)
        for (i=0; i < newTipsy->head->ndark; i++)
            fread(&newTipsy->dark[pID[i  + newTipsy->head->nsph]].vel, sizeof(float), 3, fp);
    if (newTipsy->head->nstar != 0)
        for (i=0; i < newTipsy->head->nstar; i++)
            fread(&newTipsy->star[pID[i  + newTipsy->head->nsph + newTipsy->head->ndark]].vel, sizeof(float), 3, fp);
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(float)*3){
        printf("Error: Velocity ending block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(float)*3), newTipsy->head->nbodies);
        exit(-1);
    }

    // Skip Particle IDs
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(int)){
        printf("Error: ParticleID block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(int)), newTipsy->head->nbodies);
        exit(-1);
    }
    fseek(fp, newTipsy->head->nbodies*sizeof(int), SEEK_CUR);
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(int)){
        printf("Error: ParticleID ending block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(int)), newTipsy->head->nbodies);
        exit(-1);
    }

    // The hexdump of GIZMO outputs show a full size and a half size block (properly bracketed by the size in bytes) containing all zeroes
    // followed by the mass, energy, density and smoothing length blocks
    // skip the weird enpty blocks block. idk what it's for
    fread(&blocksize, sizeof(int), 1, fp);
    fseek(fp, blocksize, SEEK_CUR);
    fread(&blocksize, sizeof(int), 1, fp);
    fread(&blocksize, sizeof(int), 1, fp);
    fseek(fp, blocksize, SEEK_CUR);
    fread(&blocksize, sizeof(int), 1, fp);

    // Read in Masses
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(float)){
        printf("Error: Mass block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(float)), newTipsy->head->nbodies);
        exit(-1);
    }
    if (newTipsy->head->nsph != 0)
        for (i=0; i < newTipsy->head->nsph; i++)
            fread(&newTipsy->gas[pID[i]].mass, sizeof(float), 1, fp);
    if (newTipsy->head->ndark != 0)
        for (i=0; i < newTipsy->head->ndark; i++)
            fread(&newTipsy->dark[pID[i  + newTipsy->head->nsph]].mass, sizeof(float), 1, fp);
    if (newTipsy->head->nstar != 0)
        for (i=0; i < newTipsy->head->nstar; i++)
            fread(&newTipsy->star[pID[i  + newTipsy->head->nsph + newTipsy->head->ndark]].mass, sizeof(float), 1, fp);
    fread(&blocksize, sizeof(int), 1, fp);
    if ((unsigned int)blocksize != newTipsy->head->nbodies*sizeof(float)){
        printf("Error: Mass ending block size [%lu particles] is inconsistent with the header [%d particles]\n", blocksize/(sizeof(float)), newTipsy->head->nbodies);
        exit(-1);
    }

    // Read in Energy
    if (newTipsy->head->nsph != 0){
        fread(&blocksize, sizeof(int), 1, fp);
        if ((unsigned int)blocksize != newTipsy->head->nsph*sizeof(float)){
            printf("Error: Energy block size [%lu gas particles] is inconsistent with the header [%d gas particles]\n", blocksize/(sizeof(float)), newTipsy->head->nsph);
            exit(-1);
        }
        for (i=0; i < newTipsy->head->nsph; i++)
            fread(&newTipsy->gas[pID[i]].temp, sizeof(float), 1, fp);
        fread(&blocksize, sizeof(int), 1, fp);
        if ((unsigned int)blocksize != newTipsy->head->nsph*sizeof(float)){
            printf("Error: Energy ending block size [%lu gas particles] is inconsistent with the header [%d gas particles]\n", blocksize/(sizeof(float)), newTipsy->head->nsph);
            exit(-1);
        }
    }

    // Read in Density
    if (newTipsy->head->nsph != 0){
        fread(&blocksize, sizeof(int), 1, fp);
        if ((unsigned int)blocksize != newTipsy->head->nsph*sizeof(float)){
            printf("Error: Density block size [%lu gas particles] is inconsistent with the header [%d gas particles]\n", blocksize/(sizeof(float)), newTipsy->head->nsph);
            exit(-1);
        }
        for (i=0; i < newTipsy->head->nsph; i++)
            fread(&newTipsy->gas[pID[i]].rho, sizeof(float), 1, fp);
        fread(&blocksize, sizeof(int), 1, fp);
        if ((unsigned int)blocksize != newTipsy->head->nsph*sizeof(float)){
            printf("Error: Density ending block size [%lu gas particles] is inconsistent with the header [%d gas particles]\n", blocksize/(sizeof(float)), newTipsy->head->nsph);
            exit(-1);
        }
    }

    // Read in Smoothing Length
    if (newTipsy->head->nsph != 0){
        fread(&blocksize, sizeof(int), 1, fp);
        if ((unsigned int)blocksize != newTipsy->head->nsph*sizeof(float)){
            printf("Error: Smoothing Length block size [%lu gas particles] is inconsistent with the header [%d gas particles]\n", blocksize/(sizeof(float)), newTipsy->head->nsph);
            exit(-1);
        }
        for (i=0; i < newTipsy->head->nsph; i++)
            fread(&newTipsy->gas[pID[i]].eps, sizeof(float), 1, fp);
        fread(&blocksize, sizeof(int), 1, fp);
        if ((unsigned int)blocksize != newTipsy->head->nsph*sizeof(float)){
            printf("Error: Smoothing Length ending block size [%lu gas particles] is inconsistent with the header [%d gas particles]\n", blocksize/(sizeof(float)), newTipsy->head->nsph);
            exit(-1);
        }
    }

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
int writeGadgetFromTipsy(const char filename[], tipsy* tipsyOut){
    /* Writes an input tipsy struct into a file with the GADGET format.

        Parameters:
            const char filename[]   - filename to write to
            tipsy* tipsyOut         - tipsy struct to write to file
        Return:
            Confirmation of write or error code?

        ToDo:
            Make return confirmation value or error code.
    */
    int i;
    unsigned int izero = 0;
    double dzero = 0;
    unsigned char bzero = 0;
    int totalsize;
    int idummy;
    double ddummy;


    FILE *fp = fopen(filename, "w");

    // Write header
    /* Header format:
        uint [6]    : number of particles of each type in the current file
        double [6]  : mass table for each particle; function ignores it
                        and just writes masses into mass block
        double      : snapshot time
        double      : redshift (cosmology only, I don't need this rn)
        int [2]     : unused
        int [6]     : number of particles of each type in the full simulation
                        (only used if there are multiple simulation files,
                        I don't need this right now so I'll just set this
                        as the same as above)
        int         : unused
        int         : number of files per snapshot (I assume 1)
        double      : boxsize (use largest of the three sizes I have stored)
        double      : Omega0 (cosmology only, I don't need this rn)
        double      : OmegaLambda (cosmology only, I don't need this rn)
        double      : HubbleParam (cosmology only, I don't need this rn)
        int[2]      : unused
        int[6]      : NumPart_Total_HW (only for 2^64 or more particles,
                        I don't need this rn)
        int         : flag for u block (0=temp, 1=entropy)
        unused      : to total 256 bytes

       Particle types:
           0 = gas
           1 = halo (dark matter in GIZMO)
           2 = disk (dummy collisionless particles in GIZMO)
           3 = bulge (dummy collisionless particles in GIZMO)
           4 = stars
           5 = boundary (black holes in GIZMO)
    */
    // HEADER BRACKET
    idummy = 256;
    fwrite(&idummy, sizeof(int), 1, fp);
    fwrite(&tipsyOut->head->nsph, sizeof(int), 1, fp);                          // uint nparticles
    fwrite(&tipsyOut->head->ndark, sizeof(int), 1, fp);
    fwrite(&izero, sizeof(int), 1, fp);
    fwrite(&izero, sizeof(int), 1, fp);
    fwrite(&tipsyOut->head->nstar, sizeof(int), 1, fp);
    fwrite(&izero, sizeof(int), 1, fp);

    for (i=0; i<6; i++){
        fwrite(&dzero, sizeof(double), 1, fp);                                  // masstable, ignore this since masses will be written explicitly in massblock
    }
    fwrite(&tipsyOut->head->simtime, sizeof(double), 1, fp);                    // double simtime
    fwrite(&dzero, sizeof(double), 1, fp);                                      // redshift
    fwrite(&izero, sizeof(int), 1, fp);
    fwrite(&izero, sizeof(int), 1, fp);

    fwrite(&tipsyOut->head->nsph, sizeof(int), 1, fp);                          // nparticles
    fwrite(&tipsyOut->head->ndark, sizeof(int), 1, fp);
    fwrite(&izero, sizeof(int), 1, fp);
    fwrite(&izero, sizeof(int), 1, fp);
    fwrite(&tipsyOut->head->nstar, sizeof(int), 1, fp);
    fwrite(&izero, sizeof(int), 1, fp);

    fwrite(&izero, sizeof(int), 1, fp);
    idummy = 1;
    fwrite(&idummy, sizeof(int), 1, fp);

    ddummy = (double)abs(tipsyOut->attr->xmax);                                 // find max boundary
    if ((double)abs(tipsyOut->attr->xmin) > ddummy)
        ddummy = (double)abs(tipsyOut->attr->xmin);
    if ((double)abs(tipsyOut->attr->ymax) > ddummy)
        ddummy = (double)abs(tipsyOut->attr->ymax);
    if ((double)abs(tipsyOut->attr->ymin) > ddummy)
        ddummy = (double)abs(tipsyOut->attr->ymin);
    if ((double)abs(tipsyOut->attr->zmax) > ddummy)
        ddummy = (double)abs(tipsyOut->attr->zmax);
    if ((double)abs(tipsyOut->attr->zmin) > ddummy)
        ddummy = (double)abs(tipsyOut->attr->zmin);
    fwrite(&ddummy, sizeof(double), 1, fp);                                     // boxsize
    for (i=0; i<3; i++){
        fwrite(&dzero, sizeof(double), 1, fp);
    }
    for (i=0; i<8; i++){
        fwrite(&izero, sizeof(int), 1, fp);
    }
    fwrite(&izero, sizeof(int), 1, fp);                                         // u flag
    totalsize = sizeof(int)*25 + sizeof(double)*12;
    printf("header excess = %f", ((float)(256-totalsize)/(float)sizeof(unsigned char)));
    for (i=0; i<(int)((256-totalsize)/sizeof(unsigned char)); i++){                   // buffer
        fwrite(&bzero, sizeof(unsigned char), 1, fp);
    }
    idummy = 256;
    fwrite(&idummy, sizeof(int), 1, fp);
    // END HEADER BRACKET

    // POSITION BRACKET
    idummy = (tipsyOut->head->nbodies*sizeof(float)*3);
    fwrite(&idummy, sizeof(int), 1, fp);
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            fwrite(&tipsyOut->gas[i].pos, sizeof(float), 3, fp);
        }
    }
    if (tipsyOut->head->ndark != 0){
        for (i=0; i < tipsyOut->head->ndark; i++){
            fwrite(&tipsyOut->dark[i].pos, sizeof(float), 3, fp);
        }
    }
    if (tipsyOut->head->nstar != 0){
        for (i=0; i < tipsyOut->head->nstar; i++){
            fwrite(&tipsyOut->star[i].pos, sizeof(float), 3, fp);
        }
    }
    idummy = (tipsyOut->head->nbodies*sizeof(float)*3);
    fwrite(&idummy, sizeof(int), 1, fp);
    // END POSITION BRACKET

    // VELOCITY BRACKET
    idummy = (tipsyOut->head->nbodies*sizeof(float)*3);
    fwrite(&idummy, sizeof(int), 1, fp);
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            fwrite(&tipsyOut->gas[i].vel, sizeof(float), 3, fp);
        }
    }
    if (tipsyOut->head->ndark != 0){
        for (i=0; i < tipsyOut->head->ndark; i++){
            fwrite(&tipsyOut->dark[i].vel, sizeof(float), 3, fp);
        }
    }
    if (tipsyOut->head->nstar != 0){
        for (i=0; i < tipsyOut->head->nstar; i++){
            fwrite(&tipsyOut->star[i].vel, sizeof(float), 3, fp);
        }
    }
    idummy = (tipsyOut->head->nbodies*sizeof(float)*3);
    fwrite(&idummy, sizeof(int), 1, fp);
    // END VELOCITY BRACKET

    // PARTICLE ID BRACKET
    idummy = (tipsyOut->head->nbodies*sizeof(int));
    fwrite(&idummy, sizeof(int), 1, fp);
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            idummy = i+1;
            fwrite(&idummy, sizeof(int), 1, fp);
        }
    }
    if (tipsyOut->head->ndark != 0){
        for (i=0; i < tipsyOut->head->ndark; i++){
            idummy = i+1+tipsyOut->head->nsph;
            fwrite(&idummy, sizeof(int), 1, fp);
        }
    }
    if (tipsyOut->head->nstar != 0){
        for (i=0; i < tipsyOut->head->nstar; i++){
            idummy = i+1+tipsyOut->head->nsph+tipsyOut->head->ndark;
            fwrite(&idummy, sizeof(int), 1, fp);
        }
    }
    idummy = (tipsyOut->head->nbodies*sizeof(int));
    fwrite(&idummy, sizeof(int), 1, fp);
    // END PARTICLE ID BRACKET

    // MASS BRACKET
    idummy = (tipsyOut->head->nbodies*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            fwrite(&tipsyOut->gas[i].mass, sizeof(float), 1, fp);
        }
    }
    if (tipsyOut->head->ndark != 0){
        for (i=0; i < tipsyOut->head->ndark; i++){
            fwrite(&tipsyOut->dark[i].mass, sizeof(float), 1, fp);
        }
    }
    if (tipsyOut->head->nstar != 0){
        for (i=0; i < tipsyOut->head->nstar; i++){
            fwrite(&tipsyOut->star[i].mass, sizeof(float), 1, fp);
        }
    }
    idummy = (tipsyOut->head->nbodies*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    // END MASS BRACKET

    // ENERGY BRACKET
    idummy = (tipsyOut->head->nsph*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            fwrite(&tipsyOut->gas[i].temp, sizeof(float), 1, fp);
        }
    }
    idummy = (tipsyOut->head->nsph*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    // END ENERGY BRACKET

    // DENSITY BRACKET
    idummy = (tipsyOut->head->nsph*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            fwrite(&tipsyOut->gas[i].rho, sizeof(float), 1, fp);
        }
    }
    idummy = (tipsyOut->head->nsph*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    // END DENSITY BRACKET

    // SMOOTHING LENGTH BRACKET
    // (maybe eps? idk it's probably just an output variable not a read-in)
    idummy = (tipsyOut->head->nsph*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    if (tipsyOut->head->nsph != 0){
        for (i=0; i < tipsyOut->head->nsph; i++){
            fwrite(&tipsyOut->gas[i].eps, sizeof(float), 1, fp);
        }
    }
    idummy = (tipsyOut->head->nbodies*sizeof(float));
    fwrite(&idummy, sizeof(int), 1, fp);
    // END SMOOTHING LENGTH BRACKET

    fclose(fp);
    return 1;
}
