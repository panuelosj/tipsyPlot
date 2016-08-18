#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../tipsyPlot.h"


float flopSetZero(float val1, float val2){
    return 0.0;
}
float flopAdd(float val1, float val2){
    return val1 + val2;
}
float flopDivide(float val1, float val2){
    return val1/val2;
}
float flopCopy(float val1, float val2){
    return val1;
}

/*
##     ##    ######## ##        #######  ########   ######
##     ##    ##       ##       ##     ## ##     ## ##    ##
##     ##    ##       ##       ##     ## ##     ## ##
##     ##    ######   ##       ##     ## ########   ######
 ##   ##     ##       ##       ##     ## ##              ##
  ## ##      ##       ##       ##     ## ##        ##    ##
   ###       ##       ########  #######  ##         ######
   Value FLoating-point OPerationS
       Arbitrary floating point operations on a particle and a value
   Performs the operation stated in the *flop function pointer on an entire
       particle, taking in particle src1 and a float src2, and saving into
       the struct specified by dest. Identical to the corresponding functions
       under pFlops (particle floating point operations), but takes in a
       particle and a float instead of two particles. The float is applied as
       the value over all attributes stored by the particle struct

       Note:
           - It is important to make sure that all three are the same particle
               type specified by the type param.
           - The (*flop) function (typedef'd) can be completely arbitrary and
               even ignore the given input. Say the operation that needs to be
               performed is [[dest = src1 /4.0]], then src2 can be set to NULL,
               and the (*flop) can be written so as to ignore src2 and simply
               use a division of 4 on src1.
           - src1 and src2 are local variables (pointers) on the stack, so
               changing them only changes where they point to and don't actually
               change values, are automatically deallocated when the function
               finishes

       Parameters:
           void* dest          - pointer to particle the results will be
                                   written to
           void* src1          - first input
           float* src2          - second input
           float (*flop)       - operation to be performed

       ToDo:
           - add check to see if all input pointers are the correct type
           - add return to check if operation was performed
*/
void vFlopGas(gas_particle* dest, gas_particle* src1, float src2, flop op){
    gas_particle dummy;
    if (dest == NULL) return;                                                   // destination cannot be NULL, unsafe
    if (src1 == NULL) src1 = &dummy;
    assert(dest != NULL && src1 != NULL);

    dest->mass         = op(src1->mass, src2);
    dest->pos[AXIS_X]  = op(src1->pos[AXIS_X], src2);
    dest->pos[AXIS_Y]  = op(src1->pos[AXIS_Y], src2);
    dest->pos[AXIS_Z]  = op(src1->pos[AXIS_Z], src2);
    dest->vel[AXIS_X]  = op(src1->vel[AXIS_X], src2);
    dest->vel[AXIS_Y]  = op(src1->vel[AXIS_Y], src2);
    dest->vel[AXIS_Z]  = op(src1->vel[AXIS_Z], src2);
    dest->rho          = op(src1->rho, src2);
    dest->temp         = op(src1->temp, src2);
    dest->eps          = op(src1->eps, src2);
    dest->metals       = op(src1->metals, src2);
    dest->phi          = op(src1->phi, src2);
}
void vFlopDark(dark_particle* dest, dark_particle* src1, float src2, flop op){
    dark_particle dummy;
    if (dest == NULL) return;                                                   // destination cannot be NULL, unsafe
    if (src1 == NULL) src1 = &dummy;
    assert(dest != NULL && src1 != NULL);

    dest->mass         = op(src1->mass, src2);
    dest->pos[AXIS_X]  = op(src1->pos[AXIS_X], src2);
    dest->pos[AXIS_Y]  = op(src1->pos[AXIS_Y], src2);
    dest->pos[AXIS_Z]  = op(src1->pos[AXIS_Z], src2);
    dest->vel[AXIS_X]  = op(src1->vel[AXIS_X], src2);
    dest->vel[AXIS_Y]  = op(src1->vel[AXIS_Y], src2);
    dest->vel[AXIS_Z]  = op(src1->vel[AXIS_Z], src2);
    dest->eps          = op(src1->eps, src2);
    dest->phi          = op(src1->phi, src2);
}
void vFlopStar(star_particle* dest, star_particle* src1, float src2, flop op){
    star_particle dummy;
    if (dest == NULL) return;                                                   // destination cannot be NULL, unsafe
    if (src1 == NULL) src1 = &dummy;
    assert(dest != NULL && src1 != NULL);

    dest->mass         = op(src1->mass, src2);
    dest->pos[AXIS_X]  = op(src1->pos[AXIS_X], src2);
    dest->pos[AXIS_Y]  = op(src1->pos[AXIS_Y], src2);
    dest->pos[AXIS_Z]  = op(src1->pos[AXIS_Z], src2);
    dest->vel[AXIS_X]  = op(src1->vel[AXIS_X], src2);
    dest->vel[AXIS_Y]  = op(src1->vel[AXIS_Y], src2);
    dest->vel[AXIS_Z]  = op(src1->vel[AXIS_Z], src2);
    dest->metals       = op(src1->metals, src2);
    dest->tform        = op(src1->tform, src2);
    dest->eps          = op(src1->eps, src2);
    dest->phi          = op(src1->phi, src2);
}


/*
########     ######## ##        #######  ########   ######
##     ##    ##       ##       ##     ## ##     ## ##    ##
##     ##    ##       ##       ##     ## ##     ## ##
########     ######   ##       ##     ## ########   ######
##           ##       ##       ##     ## ##              ##
##           ##       ##       ##     ## ##        ##    ##
##           ##       ########  #######  ##         ######
    Particle FLoating-point OPerationS
        Arbitrary floating point operations on two particles
    Performs the operation stated in the *flop function pointer on entire
        particles (performing the operation on all values enclosed by the
        particle struct), taking in particles src1 and src2, and saving into
        the struct specified by dest.

        Note:
            - It is important to make sure that all three are the same particle
                type specified by the type param.
            - The (*flop) function (typedef'd) can be completely arbitrary and
                even ignore the given input. Say the operation that needs to be
                performed is [[dest = src1 /4.0]], then src2 can be set to NULL,
                and the (*flop) can be written so as to ignore src2 and simply
                use a division of 4 on src1.
            - src1 and src2 are local variables (pointers) on the stack, so
                changing them only changes where they point to and don't actually
                change values, are automatically deallocated when the function
                finishes

        Parameters:
            void* dest          - pointer to particle the results will be
                                    written to
            void* src1          - first input
            void* src2          - second input
            float (*flop)       - operation to be performed

        ToDo:
            - add check to see if all input pointers are the correct type
            - add return to check if operation was performed
*/
void pFlopGas(gas_particle* dest, gas_particle* src1, gas_particle* src2, flop op){
    gas_particle dummy;
    if (dest == NULL) return;                                                   // destination cannot be NULL, unsafe
    if (src1 == NULL) src1 = &dummy;
    if (src2 == NULL) src2 = &dummy;
    assert(dest != NULL && src1 != NULL && src2 != NULL);

    dest->mass         = op(src1->mass, src2->mass);
    dest->pos[AXIS_X]  = op(src1->pos[AXIS_X], src2->pos[AXIS_X]);
    dest->pos[AXIS_Y]  = op(src1->pos[AXIS_Y], src2->pos[AXIS_Y]);
    dest->pos[AXIS_Z]  = op(src1->pos[AXIS_Z], src2->pos[AXIS_Z]);
    dest->vel[AXIS_X]  = op(src1->vel[AXIS_X], src2->vel[AXIS_X]);
    dest->vel[AXIS_Y]  = op(src1->vel[AXIS_Y], src2->vel[AXIS_Y]);
    dest->vel[AXIS_Z]  = op(src1->vel[AXIS_Z], src2->vel[AXIS_Z]);
    dest->rho          = op(src1->rho, src2->rho);
    dest->temp         = op(src1->temp, src2->temp);
    dest->eps          = op(src1->eps, src2->eps);
    dest->metals       = op(src1->metals, src2->metals);
    dest->phi          = op(src1->phi, src2->phi);
}
void pFlopDark(dark_particle* dest, dark_particle* src1, dark_particle* src2, flop op){
    dark_particle dummy;
    if (dest == NULL) return;                                                   // destination cannot be NULL, unsafe
    if (src1 == NULL) src1 = &dummy;
    if (src2 == NULL) src2 = &dummy;
    assert(dest != NULL && src1 != NULL && src2 != NULL);

    dest->mass         = op(src1->mass, src2->mass);
    dest->pos[AXIS_X]  = op(src1->pos[AXIS_X], src2->pos[AXIS_X]);
    dest->pos[AXIS_Y]  = op(src1->pos[AXIS_Y], src2->pos[AXIS_Y]);
    dest->pos[AXIS_Z]  = op(src1->pos[AXIS_Z], src2->pos[AXIS_Z]);
    dest->vel[AXIS_X]  = op(src1->vel[AXIS_X], src2->vel[AXIS_X]);
    dest->vel[AXIS_Y]  = op(src1->vel[AXIS_Y], src2->vel[AXIS_Y]);
    dest->vel[AXIS_Z]  = op(src1->vel[AXIS_Z], src2->vel[AXIS_Z]);
    dest->eps          = op(src1->eps, src2->eps);
    dest->phi          = op(src1->phi, src2->phi);
}
void pFlopStar(star_particle* dest, star_particle* src1, star_particle* src2, flop op){
    star_particle dummy;
    if (dest == NULL) return;                                                   // destination cannot be NULL, unsafe
    if (src1 == NULL) src1 = &dummy;
    if (src2 == NULL) src2 = &dummy;
    assert(dest != NULL && src1 != NULL && src2 != NULL);

    dest->mass         = op(src1->mass, src2->mass);
    dest->pos[AXIS_X]  = op(src1->pos[AXIS_X], src2->pos[AXIS_X]);
    dest->pos[AXIS_Y]  = op(src1->pos[AXIS_Y], src2->pos[AXIS_Y]);
    dest->pos[AXIS_Z]  = op(src1->pos[AXIS_Z], src2->pos[AXIS_Z]);
    dest->vel[AXIS_X]  = op(src1->vel[AXIS_X], src2->vel[AXIS_X]);
    dest->vel[AXIS_Y]  = op(src1->vel[AXIS_Y], src2->vel[AXIS_Y]);
    dest->vel[AXIS_Z]  = op(src1->vel[AXIS_Z], src2->vel[AXIS_Z]);
    dest->metals       = op(src1->metals, src2->metals);
    dest->tform        = op(src1->tform, src2->tform);
    dest->eps          = op(src1->eps, src2->eps);
    dest->phi          = op(src1->phi, src2->phi);
}


void pFlop(void* dest, void* src1, void* src2, int type, flop op){
    // empty pointers to be used
    gas_particle* gdest; gas_particle* g1; gas_particle* g2;
    dark_particle* ddest; dark_particle* d1; dark_particle* d2;
    star_particle* sdest; star_particle* s1; star_particle* s2;
    // dummy structs
    gas_particle gdummy; dark_particle ddummy; star_particle sdummy;

    switch (type) {
        case TYPE_GAS:
            if (dest != NULL)
                gdest = (gas_particle*)dest;
            else
                break;                                                          // destination cannot be NULL, unsafe
            if (src1 != NULL)
                g1 = (gas_particle*)src1;
            else
                g1 = &gdummy;
            if (src2 != NULL)
                g2 = (gas_particle*)src2;
            else
                g2 = &gdummy;
            assert(gdest != NULL && g1 != NULL && g2 != NULL);

            gdest->mass         = op(g1->mass, g2->mass);
            gdest->pos[AXIS_X]  = op(g1->pos[AXIS_X], g2->pos[AXIS_X]);
            gdest->pos[AXIS_Y]  = op(g1->pos[AXIS_Y], g2->pos[AXIS_Y]);
            gdest->pos[AXIS_Z]  = op(g1->pos[AXIS_Z], g2->pos[AXIS_Z]);
            gdest->vel[AXIS_X]  = op(g1->vel[AXIS_X], g2->vel[AXIS_X]);
            gdest->vel[AXIS_Y]  = op(g1->vel[AXIS_Y], g2->vel[AXIS_Y]);
            gdest->vel[AXIS_Z]  = op(g1->vel[AXIS_Z], g2->vel[AXIS_Z]);
            gdest->rho          = op(g1->rho, g2->rho);
            gdest->temp         = op(g1->temp, g2->temp);
            gdest->eps          = op(g1->eps, g2->eps);
            gdest->metals       = op(g1->metals, g2->metals);
            gdest->phi          = op(g1->phi, g2->phi);
            break;
        case TYPE_DARK:
            if (dest != NULL)
                ddest = (dark_particle*)dest;
            else
                break;                                                          // destination cannot be NULL, unsafe
            if (src1 != NULL)
                d1 = (dark_particle*)src1;
            else
                d1 = &ddummy;
            if (src2 != NULL)
                d2 = (dark_particle*)src2;
            else
                d2 = &ddummy;
            assert(ddest != NULL && d1 != NULL && d2 != NULL);

            ddest->mass         = op(d1->mass, d2->mass);
            ddest->pos[AXIS_X]  = op(d1->pos[AXIS_X], d2->pos[AXIS_X]);
            ddest->pos[AXIS_Y]  = op(d1->pos[AXIS_Y], d2->pos[AXIS_Y]);
            ddest->pos[AXIS_Z]  = op(d1->pos[AXIS_Z], d2->pos[AXIS_Z]);
            ddest->vel[AXIS_X]  = op(d1->vel[AXIS_X], d2->vel[AXIS_X]);
            ddest->vel[AXIS_Y]  = op(d1->vel[AXIS_Y], d2->vel[AXIS_Y]);
            ddest->vel[AXIS_Z]  = op(d1->vel[AXIS_Z], d2->vel[AXIS_Z]);
            ddest->eps          = op(d1->eps, d2->eps);
            ddest->phi          = op(d1->phi, d2->phi);
            break;
        case TYPE_STAR:
            if (dest != NULL)
                sdest = (star_particle*)dest;
            else
                break;                                                          // destination cannot be NULL, unsafe
            if (src1 != NULL)
                s1 = (star_particle*)src1;
            else
                s1 = &sdummy;
            if (src2 != NULL)
                s2 = (star_particle*)src2;
            else
                s2 = &sdummy;
            assert(sdest != NULL && s1 != NULL && s2 != NULL);

            sdest->mass         = op(s1->mass, s2->mass);
            sdest->pos[AXIS_X]  = op(s1->pos[AXIS_X], s2->pos[AXIS_X]);
            sdest->pos[AXIS_Y]  = op(s1->pos[AXIS_Y], s2->pos[AXIS_Y]);
            sdest->pos[AXIS_Z]  = op(s1->pos[AXIS_Z], s2->pos[AXIS_Z]);
            sdest->vel[AXIS_X]  = op(s1->vel[AXIS_X], s2->vel[AXIS_X]);
            sdest->vel[AXIS_Y]  = op(s1->vel[AXIS_Y], s2->vel[AXIS_Y]);
            sdest->vel[AXIS_Z]  = op(s1->vel[AXIS_Z], s2->vel[AXIS_Z]);
            sdest->metals       = op(s1->metals, s2->metals);
            sdest->tform        = op(s1->tform, s2->tform);
            sdest->eps          = op(s1->eps, s2->eps);
            sdest->phi          = op(s1->phi, s2->phi);
            break;
    }
}
