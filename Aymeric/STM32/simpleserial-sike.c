/*
    This file is part of the ChipWhisperer Example Targets
    Copyright (C) 2012-2017 NewAE Technology Inc.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "api.h"
#include "P434_internal.h"

/* Initial nbits, can be changed by sending 'n' followed by one byte */
#define SIKE_BOBSK3_P434_BYTES SECRETKEY_B_BYTES /* ((218 - 1 + 7) / 8) */
#define SIKE_ALICEPK_P434_BYTES 346 /* 2*3*434 bits (R0->x[0], R0->x[1], R->x[0], R->x[1], R2->x[0], R2->x[1]) + 16 bytes (MSG) */
#define SIKE_MP_ADDFAST_BYTES 112 /* 2*56 bytes */
#define INIT_NBITS OBOB_BITS - 1 /* Actually 217 since last bit is zero */

#define UINT8_RADIX 8
#define UINT8_LOG2RADIX 3

/* Experimental durations, roughly according to (Y-4.08)/0.95 for Y microseconds */
#define SLEEP_DURATION_10US 6
#define SLEEP_DURATION_100US 100
#define SLEEP_DURATION_200US 207
#define SLEEP_DURATION_1MS 1050
#define SLEEP_DURATION_10MS 10543
#define SLEEP_DURATION_100MS 105472
#define SLEEP_DURATION_1S 1054744

#define SLEEP_SEQ SLEEP_DURATION_1S
#define SLEEP_BETWEEN SLEEP_DURATION_1MS

/* Number of bits used in LADDER3PT */
static uint8_t nbits = INIT_NBITS;

/* Z = 0x0000ECEEA7BD2EDAE93254545F77410CD801A4FB559FACD4B90FF404FC00000000000000000000000000000000000000000000000000742C */
const digit_t custom_Montgomery_one[NWORDS_FIELD] = {
   0x0000742C, 0x00000000,
   0x00000000, 0x00000000,
   0x00000000, 0x00000000,
   0xFC000000, 0xB90FF404,
   0x559FACD4, 0xD801A4FB, 
   0x5F77410C, 0xE9325454,
   0xA7BD2EDA, 0x0000ECEE
};

const digit_t one[NWORDS_FIELD] = {
    0x00000001, 0x00000000,
    0x00000000, 0x00000000,
    0x00000000, 0x00000000,
    0x00000000, 0x00000000,
    0x00000000, 0x00000000,
    0x00000000, 0x00000000,
    0x00000000, 0x00000000
};

uint8_t sk[SIKE_BOBSK3_P434_BYTES] = { 0x00 };


static void custom_print(char * str)
{
    int i = 0;
    while(str[i]!= '\0')
    {
        putch(str[i++]);
    }
}

/*
 * Sleeps 1.36*duration + 4.46 microseconds.
 */
static inline void sleep(volatile uint32_t duration)
{
    while (duration-- != 0);
}

void __attribute__ ((noinline, naked))  custom_mp_addfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision addition, c = a+b.
    /*
     * Add 440-bit values in chunks of 32 bits (and propagate carry).
     *
     * Pattern (same):
     *     1) load in r2-r5 the next four values of "a" (hence, 128 bits)
     *     2) load in r6-r9 the next four values of "b" (hence, 128 bits)
     *     3) add (with carry) resp. r2-r5 with r6-r9
     *     4) stores results so far in "c"
     *
     * Instructions
     *      PUSH (ENCAPS ONTO THE STACK)
     *      POP (DECAPS FROM THE STACK)
     *      MOV (MOVES VALUES)
     *      ADDS (ADD: ADDITION, S: SET FLAGS)
     *      ADCS (AD: ADDITION, C: CARRY, S: SET FLAGS)
     *      STMIA (ST: STORE, M: MULTIPLE, IA: INCREMENT AFTER)
     *      LDMIA (LD: LOAD, M: MULTIPLE, IA: INCREMENT AFTER)
     */
    asm(

            "push  {r4-r9,lr}           \n\t"
            "mov r14, r2                \n\t" // (r2 = c) r14 <- c (pointer)

            /* pattern */
            // a is a pointer, b is a pointer
            // a[i] goes to position i pointed by a
            "ldmia r0!, {r2-r5}             \n\t" // (r0 = a) r2, r3, r4, r5 <- a[0:4]
            "ldmia r1!, {r6-r9}             \n\t" // (r1 = b) r6, r7, r8, r9 <- b[0:4] 

            // (a + b) {r2-r5} + {r6-r9}, register by register (with carry)
            "adds r2, r2, r6                \n\t"
            "adcs r3, r3, r7                \n\t"
            "adcs r4, r4, r8                \n\t"
            "adcs r5, r5, r9                \n\t"

            // stores (and updates) in (memory poitned by) r14 results so far
            "stmia r14!, {r2-r5}            \n\t"

            // repeat
            "ldmia r0!, {r2-r5}             \n\t"
            "ldmia r1!, {r6-r9}             \n\t"

            "adcs r2, r2, r6                \n\t"
            "adcs r3, r3, r7                \n\t"
            "adcs r4, r4, r8                \n\t"
            "adcs r5, r5, r9                \n\t"

            "stmia r14!, {r2-r5}            \n\t"

            // repeat
            "ldmia r0!, {r2-r5}             \n\t"
            "ldmia r1!, {r6-r9}             \n\t"

            "adcs r2, r2, r6                \n\t"
            "adcs r3, r3, r7                \n\t"
            "adcs r4, r4, r8                \n\t"
            "adcs r5, r5, r9                \n\t"

            "stmia r14!, {r2-r5}            \n\t"

            // repeat
            "ldmia r0!, {r2-r3}             \n\t"
            "ldmia r1!, {r6-r7}             \n\t"

            "adcs r2, r2, r6                \n\t"
            "adcs r3, r3, r7                \n\t"

            "stmia r14!, {r2-r3}            \n\t"

        
            "pop  {r4-r9,pc}                \n\t"
    :
    :
    :
    );

}

static void fp2_decode(const unsigned char *enc, f2elm_t x)
{ // Parse byte sequence back into GF(p^2) element, and conversion to Montgomery representation
    unsigned int i;

    for (i = 0; i < 2*(MAXBITS_FIELD / 8); i++) ((unsigned char *)x)[i] = 0;
    for (i = 0; i < FP2_ENCODED_BYTES / 2; i++) {
        ((unsigned char*)x)[i] = enc[i];
        ((unsigned char*)x)[i + MAXBITS_FIELD / 8] = enc[i + FP2_ENCODED_BYTES / 2];
    }
    to_fp2mont(x, x);
}


static void swap_points(point_proj_t P, point_proj_t Q, const digit_t option)
{ // Swap points.
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->X[0][i] ^ Q->X[0][i]);
        P->X[0][i] = temp ^ P->X[0][i];
        Q->X[0][i] = temp ^ Q->X[0][i];
        temp = option & (P->Z[0][i] ^ Q->Z[0][i]);
        P->Z[0][i] = temp ^ P->Z[0][i];
        Q->Z[0][i] = temp ^ Q->Z[0][i];
        temp = option & (P->X[1][i] ^ Q->X[1][i]);
        P->X[1][i] = temp ^ P->X[1][i];
        Q->X[1][i] = temp ^ Q->X[1][i];
        temp = option & (P->Z[1][i] ^ Q->Z[1][i]);
        P->Z[1][i] = temp ^ P->Z[1][i];
        Q->Z[1][i] = temp ^ Q->Z[1][i];
    }
}

void custom_fp2mul434_mont(const f2elm_t a, const f2elm_t b, f2elm_t c)
{ // GF(p^2) multiplication using Montgomery arithmetic, c = a*b in GF(p^2).
  // Inputs: a = a0+a1*i and b = b0+b1*i, where a0, a1, b0, b1 are in [0, 2*p-1]
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1]
    felm_t t1, t2;
    dfelm_t tt2;             //#fixme
    //dfelm_t tt1, tt2, tt3; #fixme
    //digit_t mask;          #fixme
    //unsigned int i;        #fixme

    sleep(SLEEP_BETWEEN);
    trigger_high();
    custom_mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1
    trigger_low();
    sleep(SLEEP_BETWEEN);
    sleep(SLEEP_BETWEEN);
    trigger_high();
    custom_mp_addfast(b[0], b[1], t2);                      // t2 = b0+b1
    trigger_low();
    sleep(SLEEP_BETWEEN);

    fpmul434_mont(a[0], b[0], c[0]);
    fpmul434_mont(a[1], b[1], tt2);
    fpmul434_mont(t1, t2, c[1]);

    fpsub434(c[1],c[0],c[1]);
    fpsub434(c[1],tt2,c[1]);

    fpsub434(c[0],tt2,c[0]);
}

void custom_fp2sqr_mont(const f2elm_t a, f2elm_t c)
{ // GF(p^2) squaring using Montgomery arithmetic, c = a^2 in GF(p^2).
  // Inputs: a = a0+a1*i, where a0, a1 are in [0, 2*p-1]
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1]
    felm_t t1, t2, t3;

    sleep(SLEEP_BETWEEN);
    trigger_high();
    custom_mp_addfast(a[0], a[1], t1);                  // t1 = a0+a1
    trigger_low();
    sleep(SLEEP_BETWEEN);
    fpsub434(a[0], a[1], t2);                           // t2 = a0-a1
    sleep(SLEEP_BETWEEN);
    trigger_high();
    custom_mp_addfast(a[0], a[0], t3);                  // t3 = 2a0
    trigger_low();
    sleep(SLEEP_BETWEEN);
    fpmul434_mont(t1, t2, c[0]);                        // c0 = (a0+a1)(a0-a1)
    fpmul434_mont(t3, a[1], c[1]);                      // c1 = 2a0*a1
}

void custom_xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP.
    f2elm_t t0, t1, t2;

    fpadd434(P->X, P->Z, t0);                         //  1. t0 = XP+ZP
    fpsub434(P->X, P->Z, t1);                         //  2. t1 = XP-ZP
    custom_fp2sqr_mont(t0, P->X);                     //  3. XP = (XP+ZP)^2
    fpsub434(Q->X, Q->Z, t2); fp2correction434(t2);   //  4. t2 = XQ-ZQ
    fpadd434(Q->X, Q->Z, Q->X);                       //  5. XQ = XQ+ZQ
    custom_fp2mul434_mont(t0, t2, t0);                //  6. t0 = (XP+ZP)*(XQ-ZQ)
    custom_fp2sqr_mont(t1, P->Z);                     //  7. ZP = (XP-ZP)^2
    custom_fp2mul434_mont(t1, Q->X, t1);              //  8. t1 = (XP-ZP)*(XQ+ZQ)
    fpsub434(P->X, P->Z, t2);                         //  9. t2 = (XP+ZP)^2-(XP-ZP)^2
    custom_fp2mul434_mont(P->X, P->Z, P->X);          // 10. XP = (XP+ZP)^2*(XP-ZP)^2
    custom_fp2mul434_mont(t2, A24, Q->X);             // 11. XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    fpsub434(t0, t1, Q->Z);                           // 12. ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    fpadd434(Q->X, P->Z, P->Z);                       // 13. ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    fpadd434(t0, t1, Q->X);                           // 14. XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    custom_fp2mul434_mont(P->Z, t2, P->Z);            // 15. ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    custom_fp2sqr_mont(Q->Z, Q->Z);                   // 16. ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    custom_fp2sqr_mont(Q->X, Q->X);                   // 17. XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    custom_fp2mul434_mont(Q->Z, xPQ, Q->Z);           // 18. ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}

uint8_t get_nbits(uint8_t* n)
{
    nbits =  n[0];
    return 0x00;
}

uint8_t get_key(uint8_t* k)
{
    for (int i=0; i < SIKE_BOBSK3_P434_BYTES; ++i)
    {
        sk[i] = k[i];
    }
    return 0x00;
}

uint8_t do_mp_addfast(uint8_t* pt)
{
    felm_t tmp;

    trigger_high();
    custom_mp_addfast((digit_t*) (pt), (digit_t*) (pt + 56), tmp);
    trigger_low();

    return 0x00;
}

uint8_t get_pt(uint8_t* pt)
{
    int i = 0, bit = 0, swap = 0, prevbit = 0;
    point_proj_t R = {0};
    point_proj_t R2 = {0};
    point_proj_t R0 = {0};
    digit_t mask;
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0};
    f2elm_t A24 = {0};

    // Initialize images of Alice's basis
    fp2_decode(pt, R->X); /* 0:110 */
    fp2_decode(pt + FP2_ENCODED_BYTES, R0->X); /* 110:220 */
    fp2_decode(pt + 2*FP2_ENCODED_BYTES, R2->X); /* 220:330 */

    // Initialize constants: A24plus = A+2C, A24minus = A-2C, where C=1
    get_A(R->X, R0->X, R2->X, A);
    fpadd434((digit_t*)&custom_Montgomery_one, (digit_t*)&custom_Montgomery_one, A24minus[0]);
    fp2add434(A, A24minus, A24plus);
    fp2sub434(A, A24minus, A24minus);

    fpcopy434((digit_t*)&custom_Montgomery_one, A24[0]);
    fp2add434(A24, A24, A24);
    fp2add434(A, A24, A24);
    fp2div2_434(A24, A24);
    fp2div2_434(A24, A24); 

    fpcopy434((digit_t*)&custom_Montgomery_one, (digit_t*)R2->Z);
    fpcopy434((digit_t*)&custom_Montgomery_one, (digit_t*)R->Z);
    fpcopy434((digit_t*)&custom_Montgomery_one, (digit_t*)R0->Z);
    fpzero434((digit_t*)(R->Z)[1]);

    for (i = 0; i < nbits; i++) {
        bit = (sk[i >> UINT8_LOG2RADIX] >> (i & (UINT8_RADIX - 1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(R, R2, mask);
        custom_xDBLADD(R0, R2, R->X, A24);
        custom_fp2mul434_mont(R2->X, R->Z, R2->X);
        sleep(SLEEP_SEQ);
    }

    return 0x00;
}

uint8_t test_trig(uint8_t* x)
{
    trigger_high();
    sleep(((x[3] << 24) | (x[2] << 16) | (x[1] << 8) | x[0]));
    trigger_low();

    return 0x00;
}

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    /* Prints hello */
    custom_print("Hello\n");

    simpleserial_init();

    simpleserial_addcmd('k', SIKE_BOBSK3_P434_BYTES, get_key);
    simpleserial_addcmd('p', SIKE_ALICEPK_P434_BYTES, get_pt);
    simpleserial_addcmd('q', SIKE_MP_ADDFAST_BYTES, do_mp_addfast);
    simpleserial_addcmd('n', 1, get_nbits);
    simpleserial_addcmd('t', 4, test_trig);

    while(1)
        simpleserial_get();
}
