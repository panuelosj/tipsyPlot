#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../tipsyPlot.h"

void particleFlop(void* dest, void* src1, void* src2, int type, float (*flop)(float val1, float val2)){
    /* Performs the operation stated in the *flop function pointer on entire
        particles (performing the operation on all values enclosed by the
        particle struct), taking in particles src1 and src2, and saving into
        the particle specified by dest.

        Note:
            - It is important to make sure that all three are the same particle
                type specified by the type param.
            - The (*flop) function can be completely arbitrary and even ignore
                the given input. Say the operation that needs to be performed
                is [[dest = src1 /4.0]], then src2 can be set to NULL, and the
                (*flop) can be written so as to ignore val2 and simply use a
                division of 4 on val1.

        Parameters:
            void* dest          - pointer to particle the results will be
                                    written to
            void* src1          - first input
            void* src2          - second input
            int type            - type of particle
            float (*flop)       - operation to be performed

        ToDo:
            - add check to see if all input pointers are the correct type
    */

    gas_particle* gdest; gas_particle* g1; gas_particle* g2;
    dark_particle* ddest; dark_particle* d1; dark_particle* d2;
    star_particle* sdest; star_particle* s1; star_particle* s2;
    switch (type) {
        case TYPE_GAS:
            gdest   = (gas_particle*)dest;
            g1      = (gas_particle*)src1;
            g2      = (gas_particle*)src2;
            gdest->mass         = (*flop)(g1->mass, g2->mass);
            gdest->pos[AXIS_X]  = (*flop)(g1->pos[AXIS_X], g2->pos[AXIS_X]);
            gdest->pos[AXIS_Y]  = (*flop)(g1->pos[AXIS_Y], g2->pos[AXIS_Y]);
            gdest->pos[AXIS_Z]  = (*flop)(g1->pos[AXIS_Z], g2->pos[AXIS_Z]);
            gdest->vel[AXIS_X]  = (*flop)(g1->vel[AXIS_X], g2->vel[AXIS_X]);
            gdest->vel[AXIS_Y]  = (*flop)(g1->vel[AXIS_Y], g2->vel[AXIS_Y]);
            gdest->vel[AXIS_Z]  = (*flop)(g1->vel[AXIS_Z], g2->vel[AXIS_Z]);
            gdest->rho          = (*flop)(g1->rho, g2->rho);
            gdest->temp         = (*flop)(g1->temp, g2->temp);
            gdest->eps          = (*flop)(g1->eps, g2->eps);
            gdest->metals       = (*flop)(g1->metals, g2->metals);
            gdest->phi          = (*flop)(g1->phi, g2->phi);
            break;
        case TYPE_DARK:
            ddest   = (dark_particle*)dest;
            d1      = (dark_particle*)src1;
            d2      = (dark_particle*)src2;
            ddest->mass         = (*flop)(d1->mass, d2->mass);
            ddest->pos[AXIS_X]  = (*flop)(d1->pos[AXIS_X], d2->pos[AXIS_X]);
            ddest->pos[AXIS_Y]  = (*flop)(d1->pos[AXIS_Y], d2->pos[AXIS_Y]);
            ddest->pos[AXIS_Z]  = (*flop)(d1->pos[AXIS_Z], d2->pos[AXIS_Z]);
            ddest->vel[AXIS_X]  = (*flop)(d1->vel[AXIS_X], d2->vel[AXIS_X]);
            ddest->vel[AXIS_Y]  = (*flop)(d1->vel[AXIS_Y], d2->vel[AXIS_Y]);
            ddest->vel[AXIS_Z]  = (*flop)(d1->vel[AXIS_Z], d2->vel[AXIS_Z]);
            ddest->eps          = (*flop)(d1->eps, d2->eps);
            ddest->phi          = (*flop)(d1->phi, d2->phi);
            break;
        case TYPE_STAR:
            sdest   = (star_particle*)dest;
            s1      = (star_particle*)src1;
            s2      = (star_particle*)src2;
            sdest->mass         = (*flop)(s1->mass, s2->mass);
            sdest->pos[AXIS_X]  = (*flop)(s1->pos[AXIS_X], s2->pos[AXIS_X]);
            sdest->pos[AXIS_Y]  = (*flop)(s1->pos[AXIS_Y], s2->pos[AXIS_Y]);
            sdest->pos[AXIS_Z]  = (*flop)(s1->pos[AXIS_Z], s2->pos[AXIS_Z]);
            sdest->vel[AXIS_X]  = (*flop)(s1->vel[AXIS_X], s2->vel[AXIS_X]);
            sdest->vel[AXIS_Y]  = (*flop)(s1->vel[AXIS_Y], s2->vel[AXIS_Y]);
            sdest->vel[AXIS_Z]  = (*flop)(s1->vel[AXIS_Z], s2->vel[AXIS_Z]);
            sdest->metals       = (*flop)(s1->metals, s2->metals);
            sdest->tform        = (*flop)(s1->tform, s2->tform);
            sdest->eps          = (*flop)(s1->eps, s2->eps);
            sdest->phi          = (*flop)(s1->phi, s2->phi);
            break;
    }
}

void particleSetZero(void* particle, int type){
    /* Fills in all values of the input particle with zero. The type must be
        explicitly stated for the function to know how to treat the void pointer

        Parameters:
            void* particle      - pointer to the particle to zero out
            int type            - type of particle
    */
    gas_particle* g; dark_particle* d; star_particle* s;
    switch (type) {
        case TYPE_GAS:
            g = (gas_particle*)particle;
            g->mass         = 0.0;
            g->pos[AXIS_X]  = 0.0;
            g->pos[AXIS_Y]  = 0.0;
            g->pos[AXIS_Z]  = 0.0;
            g->vel[AXIS_X]  = 0.0;
            g->vel[AXIS_Y]  = 0.0;
            g->vel[AXIS_Z]  = 0.0;
            g->rho          = 0.0;
            g->temp         = 0.0;
            g->eps          = 0.0;
            g->metals       = 0.0;
            g->phi          = 0.0;
            break;
        case TYPE_DARK:
            d = (dark_particle*)particle;
            d->mass         = 0.0;
            d->pos[AXIS_X]  = 0.0;
            d->pos[AXIS_Y]  = 0.0;
            d->pos[AXIS_Z]  = 0.0;
            d->vel[AXIS_X]  = 0.0;
            d->vel[AXIS_Y]  = 0.0;
            d->vel[AXIS_Z]  = 0.0;
            d->eps          = 0.0;
            d->phi          = 0.0;
            break;
        case TYPE_STAR:
            s = (star_particle*)particle;
            s->mass         = 0.0;
            s->pos[AXIS_X]  = 0.0;
            s->pos[AXIS_Y]  = 0.0;
            s->pos[AXIS_Z]  = 0.0;
            s->vel[AXIS_X]  = 0.0;
            s->vel[AXIS_Y]  = 0.0;
            s->vel[AXIS_Z]  = 0.0;
            s->metals       = 0.0;
            s->tform        = 0.0;
            s->eps          = 0.0;
            s->phi          = 0.0;
            break;
    }
}

void particleAdd(void* dest, void* src1, void* src2, int type){
    /* Adds the values of the partices in src1 and src2 and saved them to dest.
            [[dest = src1 + src2]]
        src1 and src2 remains unchanged while values at dest are updated.
        dest and src1 can be the same, then the operation simply becomes
            [[dest += src2]]
        same for src 2. All three values can be the same, then the operation becomes
            [[dest *= 2.0]]
        No new memory allocation is performed.

        Note:
            - It is important to make sure that all three are the same particle
                type specified by the type param.

        Parameters:
            void* dest          - pointer to particle the results will be
                                    written to
            void* src1          - first addend
            void* src2          - second addend
            int type            - type of particle

        ToDo:
            - add check to see if all input pointers are the correct type
    */

    gas_particle* gdest; gas_particle* g1; gas_particle* g2;
    dark_particle* ddest; dark_particle* d1; dark_particle* d2;
    star_particle* sdest; star_particle* s1; star_particle* s2;
    switch (type) {
        case TYPE_GAS:
            gdest   = (gas_particle*)dest;
            g1      = (gas_particle*)src1;
            g2      = (gas_particle*)src2;
            gdest->mass         = g1->mass + g2->mass;
            gdest->pos[AXIS_X]  = g1->pos[AXIS_X] + g2->pos[AXIS_X];
            gdest->pos[AXIS_Y]  = g1->pos[AXIS_Y] + g2->pos[AXIS_Y];
            gdest->pos[AXIS_Z]  = g1->pos[AXIS_Z] + g2->pos[AXIS_Z];
            gdest->vel[AXIS_X]  = g1->vel[AXIS_X] + g2->vel[AXIS_X];
            gdest->vel[AXIS_Y]  = g1->vel[AXIS_Y] + g2->vel[AXIS_Y];
            gdest->vel[AXIS_Z]  = g1->vel[AXIS_Z] + g2->vel[AXIS_Z];
            gdest->rho          = g1->rho + g2->rho;
            gdest->temp         = g1->temp + g2->temp;
            gdest->eps          = g1->eps + g2->eps;
            gdest->metals       = g1->metals + g2->metals;
            gdest->phi          = g1->phi + g2->phi;
            break;
        case TYPE_DARK:
            ddest   = (dark_particle*)dest;
            d1      = (dark_particle*)src1;
            d2      = (dark_particle*)src2;
            ddest->mass         = d1->mass + d2->mass;
            ddest->pos[AXIS_X]  = d1->pos[AXIS_X] + d2->pos[AXIS_X];
            ddest->pos[AXIS_Y]  = d1->pos[AXIS_Y] + d2->pos[AXIS_Y];
            ddest->pos[AXIS_Z]  = d1->pos[AXIS_Z] + d2->pos[AXIS_Z];
            ddest->vel[AXIS_X]  = d1->vel[AXIS_X] + d2->vel[AXIS_X];
            ddest->vel[AXIS_Y]  = d1->vel[AXIS_Y] + d2->vel[AXIS_Y];
            ddest->vel[AXIS_Z]  = d1->vel[AXIS_Z] + d2->vel[AXIS_Z];
            ddest->eps          = d1->eps + d2->eps;
            ddest->phi          = d1->phi + d2->phi;
            break;
        case TYPE_STAR:
            sdest   = (star_particle*)dest;
            s1      = (star_particle*)src1;
            s2      = (star_particle*)src2;
            sdest->mass         = s1->mass + s2->mass;
            sdest->pos[AXIS_X]  = s1->pos[AXIS_X] + s2->pos[AXIS_X];
            sdest->pos[AXIS_Y]  = s1->pos[AXIS_Y] + s2->pos[AXIS_Y];
            sdest->pos[AXIS_Z]  = s1->pos[AXIS_Z] + s2->pos[AXIS_Z];
            sdest->vel[AXIS_X]  = s1->vel[AXIS_X] + s2->vel[AXIS_X];
            sdest->vel[AXIS_Y]  = s1->vel[AXIS_Y] + s2->vel[AXIS_Y];
            sdest->vel[AXIS_Z]  = s1->vel[AXIS_Z] + s2->vel[AXIS_Z];
            sdest->metals       = s1->metals + s2->metals;
            sdest->tform        = s1->tform + s2->tform;
            sdest->eps          = s1->eps + s2->eps;
            sdest->phi          = s1->phi + s2->phi;
            break;
    }
}
