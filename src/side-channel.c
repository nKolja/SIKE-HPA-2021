/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny key encapsulation (SIKE) protocol
*********************************************************************************************/ 

#include <string.h>
#include "sha3/fips202.h"



//ADDED FOR TESTING PURPOSES
#include "gmp.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>



unsigned char FromHex(char c)
    {
    switch(c)
        {
        case '0': return 0;
        case '1': return 1;
        case '2': return 2;
        case '3': return 3;
        case '4': return 4;
        case '5': return 5;
        case '6': return 6;
        case '7': return 7;
        case '8': return 8;
        case '9': return 9;
        case 'a': return 10;
        case 'b': return 11;
        case 'c': return 12;
        case 'd': return 13;
        case 'e': return 14;
        case 'f': return 15;
        }
    // Report a problem here!
    return -1;
    }




static int hamming_weight(uint32_t i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (int)((((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24);
}




void read_bob_pk(char* bob_pk_path, unsigned char* bob_pk)
{
    FILE *fp = fopen(bob_pk_path, "rb");
    assert((fp != NULL) && "Error opening Bob's public key file");

    assert( fread(bob_pk, CRYPTO_PUBLICKEYBYTES, 1, fp) != 0 );
    fclose(fp);
}

void read_alice_pk(char* alice_pk_path, unsigned char* alice_pk)
{
    FILE *fp = fopen(alice_pk_path, "rb");
    assert((fp != NULL) && "Error opening Alice's public key file");

    assert( fread(alice_pk, CRYPTO_CIPHERTEXTBYTES, 1, fp) != 0 );
    fclose(fp);
}

void write_bob_pk(char *bob_pk_path, unsigned char* bob_pk)
{
    FILE *fp = fopen(bob_pk_path, "wb");
    assert((fp != NULL) && "Error opening Bob's public key file");

    fwrite((unsigned char *)bob_pk, CRYPTO_PUBLICKEYBYTES, 1, fp);
    fclose(fp);
}

void write_alice_pk(char *alice_pk_path, unsigned char* alice_pk)
{
    FILE *fp = fopen(alice_pk_path, "wb");
    assert((fp != NULL) && "Error opening Alice's public key file");

    fwrite((unsigned char *)alice_pk, CRYPTO_CIPHERTEXTBYTES, 1, fp);
    fclose(fp);
}


void write_bob_sk(char *bob_sk_path, unsigned char* bob_sk)
{
    FILE *fp = fopen(bob_sk_path, "wb");
    assert((fp != NULL) && "Error opening Bob's secret key file");

    fwrite((unsigned char *)bob_sk, CRYPTO_SECRETKEYBYTES, 1, fp);
    fclose(fp);
}


void write_only_bob_sk(char *bob_sk_path, unsigned char* bob_sk)
{
    FILE *fp = fopen(bob_sk_path, "wb");
    assert((fp != NULL) && "Error opening Bob's secret key file");

    fwrite((unsigned char *)bob_sk, SECRETKEY_B_BYTES, 1, fp);
    fclose(fp);
}


int create_empty_bob_sk(unsigned char *bob_pk, unsigned char *bob_sk)
{ 
    // SIKE's key generation
    // Outputs: secret key bob_sk (CRYPTO_SECRETKEYBYTES = MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes)
    //          public key bob_pk (CRYPTO_PUBLICKEYBYTES bytes) 

    // Generate lower portion of secret key bob_sk <- s||SK
    randombytes(bob_sk, MSG_BYTES);                        //WRITE RANDOMNESS IN THE BEGGINING
    memset(bob_sk+MSG_BYTES, 0x00, SECRETKEY_B_BYTES);     //WRITE ALL ZEROS -> THIS IS WHAT WE'RE SEARCHING FOR

    // Append public key bob_pk to secret key bob_sk
    memcpy(&bob_sk[MSG_BYTES + SECRETKEY_B_BYTES], bob_pk, CRYPTO_PUBLICKEYBYTES);    //APPEND THE PUBLIC KEY THAT WAS READ FROM FILE

    return 0;
}




void write_points(char *point_path, point_proj_t P0, point_proj_t P1, point_proj_t P2)
{
    FILE *fp = fopen(point_path, "ab");
    assert((fp != NULL) && "Error writing point file");

    fwrite((point_proj *)P0, sizeof(point_proj), 1, fp);
    fwrite((point_proj *)P1, sizeof(point_proj), 1, fp);
    fwrite((point_proj *)P2, sizeof(point_proj), 1, fp);
    fclose(fp);
}


void read_points(char *point_path, point_proj_t P0, point_proj_t P1, point_proj_t P2)
{
    FILE *fp = fopen(point_path, "rb");
    assert((fp != NULL) && "Error reading point file");

    assert( fread((point_proj *)P0, sizeof(point_proj), 1, fp) != 0 );
    assert( fread((point_proj *)P1, sizeof(point_proj), 1, fp) != 0 );
    assert( fread((point_proj *)P2, sizeof(point_proj), 1, fp) != 0 );
    fclose(fp);
}






void generate_alice_sk_pk(unsigned char* bob_pk, unsigned char* temp, unsigned char* alice_sk, unsigned char* alice_pk, unsigned char* shared_secret, unsigned char* jinvariant, unsigned char* h)
{
    // SIKE's encapsulation
    // Input:   public key          bob_pk          (CRYPTO_PUBLICKEYBYTES bytes)
    // Outputs: shared secret       shared_secret      (CRYPTO_BYTES bytes)
    //          ciphertext message  alice_pk       (CRYPTO_CIPHERTEXTBYTES = CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes)

    // Generate alice_sk <- G(temp||bob_pk) mod oA 
    randombytes(temp, MSG_BYTES);                                                       //Generate random bytes
    memcpy(&temp[MSG_BYTES], bob_pk, CRYPTO_PUBLICKEYBYTES);                            //Concatenate random bytes with BOB's public key
    shake256(alice_sk, SECRETKEY_A_BYTES, temp, CRYPTO_PUBLICKEYBYTES+MSG_BYTES);       //Hash the above
    alice_sk[SECRETKEY_A_BYTES - 1] &= MASK_ALICE;                                      //Pad correctly

    // Encrypt
    EphemeralKeyGeneration_A(alice_sk, alice_pk);                                       //Generate Alice's public key
    EphemeralSecretAgreement_A(alice_sk, bob_pk, jinvariant);                           //Generate the common j invariant
    shake256(h, MSG_BYTES, jinvariant, FP2_ENCODED_BYTES);                              //Hash the j invariant
    for (int i = 0; i < MSG_BYTES; i++) 
    {
        alice_pk[i + CRYPTO_PUBLICKEYBYTES] = temp[i] ^ h[i];                           //XOR hashed j inv with the random bytes and write it to the public key
    }

    // Generate shared secret shared_secret <- H(temp||alice_pk)
    memcpy(&temp[MSG_BYTES], alice_pk, CRYPTO_CIPHERTEXTBYTES);                         //Concatenate random bytes with Alice's public key
    shake256(shared_secret, CRYPTO_BYTES, temp, CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);      //Hash to obtain the common shared secret K
}







void print_felm(const felm_t a)
{
    for (int i = 0; i < 4*14; i++)
        printf("%02x", ((uint8_t *)(a))[i]);
    printf("\n");
}



void extract_hamming_weight(const felm_t a, const felm_t b, const felm_t c,  int hamming_weights[3][14])
{
    //Compute hamming weights of each 32bit word of a field element (tot 14)
    // To be more precise the Hamming distance from the previous word in the pipeline register is computed.
    uint32_t hwa0, hwa1, hwa2, hwa3, hwb0, hwb1, hwb2, hwb3, hwc0, hwc1, hwc2, hwc3;
    uint32_t backup=0;



    for(int i = 0; i < 12; i+= 4)
    {
        hwa0 = ((uint32_t *) (a))[i+0];
        hwa1 = ((uint32_t *) (a))[i+1];
        hwa2 = ((uint32_t *) (a))[i+2];
        hwa3 = ((uint32_t *) (a))[i+3];

        hwb0 = ((uint32_t *) (b))[i+0];
        hwb1 = ((uint32_t *) (b))[i+1];
        hwb2 = ((uint32_t *) (b))[i+2];
        hwb3 = ((uint32_t *) (b))[i+3];

        hwc0 = ((uint32_t *) (c))[i+0];
        hwc1 = ((uint32_t *) (c))[i+1];
        hwc2 = ((uint32_t *) (c))[i+2];
        hwc3 = ((uint32_t *) (c))[i+3];


        hamming_weights[0][i+0] = hamming_weight(hwa0^backup);
        hamming_weights[0][i+1] = hamming_weight(hwa1^hwa0);
        hamming_weights[0][i+2] = hamming_weight(hwa2^hwa1);
        hamming_weights[0][i+3] = hamming_weight(hwa3^hwa2);

        hamming_weights[1][i+0] = hamming_weight(hwb0^hwa3);
        hamming_weights[1][i+1] = hamming_weight(hwb1^hwb0);
        hamming_weights[1][i+2] = hamming_weight(hwb2^hwb1);
        hamming_weights[1][i+3] = hamming_weight(hwb3^hwb2);

        hamming_weights[2][i+0] = hamming_weight(hwc0^hwb3);
        hamming_weights[2][i+1] = hamming_weight(hwc1^hwc0);
        hamming_weights[2][i+2] = hamming_weight(hwc2^hwc1);
        hamming_weights[2][i+3] = hamming_weight(hwc3^hwc2);

        backup = hwc3;
    }



    hwa0 = ((uint32_t *) (a))[12];
    hwa1 = ((uint32_t *) (a))[13];

    hwb0 = ((uint32_t *) (b))[12];
    hwb1 = ((uint32_t *) (b))[13];

    hwc0 = ((uint32_t *) (c))[12];
    hwc1 = ((uint32_t *) (c))[13];


    hamming_weights[0][12] = hamming_weight(hwa0^backup);
    hamming_weights[0][13] = hamming_weight(hwa1^hwa0);

    hamming_weights[1][12] = hamming_weight(hwb0^hwa1);
    hamming_weights[1][13] = hamming_weight(hwb1^hwb0);

    hamming_weights[2][12] = hamming_weight(hwc0);
    hamming_weights[2][13] = hamming_weight(hwc1^hwc0);


}




void fp2mul_mont_weights(const f2elm_t a, const f2elm_t b, f2elm_t c, int hamming_weights_0[3][14], int hamming_weights_1[3][14])
{ // GF(p^2) multiplication using Montgomery arithmetic, c = a*b in GF(p^2).
  // Inputs: a = a0+a1*i and b = b0+b1*i, where a0, a1, b0, b1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2;
    dfelm_t tt1, tt2, tt3; 


    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1
    extract_hamming_weight(a[0], a[1], t1, hamming_weights_0);

    mp_addfast(b[0], b[1], t2);                      // t2 = b0+b1
    extract_hamming_weight(b[0], b[1], t2, hamming_weights_1);

    fpmul_mont(a[0], b[0], c[0]);                     // tt1 = a0*b0
    fpmul_mont(a[1], b[1], tt2);                     // tt2 = a1*b1
    fpmul_mont(t1, t2, c[1]);                         // tt3 = (a0+a1)*(b0+b1)


    fpsub(c[1], c[0], c[1]);
    fpsub(c[1], tt2, c[1]);

    fpsub(c[0], tt2, c[0]);


    // mp_mul(a[0], b[0], tt1, NWORDS_FIELD);           // tt1 = a0*b0
    // mp_mul(a[1], b[1], tt2, NWORDS_FIELD);           // tt2 = a1*b1
    // mp_mul(t1, t2, tt3, NWORDS_FIELD);               // tt3 = (a0+a1)*(b0+b1)
    // mp_dblsubfast(tt1, tt2, tt3);                    // tt3 = (a0+a1)*(b0+b1) - a0*b0 - a1*b1
    // mp_subaddfast(tt1, tt2, tt1);                    // tt1 = a0*b0 - a1*b1 + p*2^MAXBITS_FIELD if a0*b0 - a1*b1 < 0, else tt1 = a0*b0 - a1*b1
    // rdc_mont(tt3, c[1]);                             // c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 
    // rdc_mont(tt1, c[0]);                             // c[0] = a0*b0 - a1*b1
}


void fp2sqr_mont_weights(const f2elm_t a, f2elm_t c, int hamming_weights_0[3][14], int hamming_weights_1[3][14])
{ // GF(p^2) squaring using Montgomery arithmetic, c = a^2 in GF(p^2).
  // Inputs: a = a0+a1*i, where a0, a1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2, t3;
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1 
    extract_hamming_weight(a[0], a[1], t1, hamming_weights_0);

    // sub_p4(a[0], a[1], t2);                          // t2 = a0-a1
    fpsub(a[0], a[1], t2);                          // t2 = a0-a1

    mp_addfast(a[0], a[0], t3);                      // t3 = 2a0
    extract_hamming_weight(a[0], a[0], t3, hamming_weights_1);

    fpmul_mont(t1, t2, c[0]);                        // c0 = (a0+a1)(a0-a1)
    fpmul_mont(t3, a[1], c[1]);                      // c1 = 2a0*a1
}


void xDBLADD_weights(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24, int hamming_weights[10][3][14])
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP. 
    f2elm_t t0, t1, t2;

    //32-bit ARM code    
    fp2add(P->X, P->Z, t0);                                                         // t0 = XP+ZP
    fp2sub(P->X, P->Z, t1);                                                         // t1 = XP-ZP
    fp2sqr_mont(t0, P->X);                                                          // XP = (XP+ZP)^2
    fp2sub(Q->X, Q->Z, t2);                                                         // t2 = XQ-ZQ
    fp2correction(t2);
    fp2add(Q->X, Q->Z, Q->X);                                                       // XQ = XQ+ZQ

    //OLD 64-bit code
    // mp2_add(P->X, P->Z, t0);                        // t0 = XP+ZP
    // mp2_sub_p2(P->X, P->Z, t1);                     // t1 = XP-ZP
    // fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2
    // mp2_sub_p2(Q->X, Q->Z, t2);                     // t2 = XQ-ZQ
    // mp2_add(Q->X, Q->Z, Q->X);                      // XQ = XQ+ZQ



    // 0          // old 2-3
    fp2mul_mont_weights(t0, t2, t0, hamming_weights[0], hamming_weights[0]);        // t0 = (XP+ZP)*(XQ-ZQ)  

    fp2sqr_mont(t1, P->Z);                                                          // ZP = (XP-ZP)^2

    // 1          // old 6-7
    fp2mul_mont_weights(t1, Q->X, t1, hamming_weights[1], hamming_weights[1]);      // t1 = (XP-ZP)*(XQ+ZQ)

    //32-bit ARM code
    fp2sub(P->X, P->Z, t2);                                                     // t2 = (XP+ZP)^2-(XP-ZP)^2
    fp2mul_mont(P->X, P->Z, P->X);                                                  // XP = (XP+ZP)^2*(XP-ZP)^2
    fp2mul_mont(A24, t2, Q->X);                                                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sub(t0, t1, Q->Z);                                                       // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    fp2add(Q->X, P->Z, P->Z);                                                      // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    fp2add(t0, t1, Q->X);                                                          // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    fp2mul_mont(P->Z, t2, P->Z);                                                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]

    // OLD 64-bit code
    // mp2_sub_p2(P->X, P->Z, t2);                     // t2 = (XP+ZP)^2-(XP-ZP)^2
    // fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2
    // fp2mul_mont(A24, t2, Q->X);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    // mp2_sub_p2(t0, t1, Q->Z);                       // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    // mp2_add(Q->X, P->Z, P->Z);                      // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    // mp2_add(t0, t1, Q->X);                          // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    // fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]

    // 2-3        // old 14-15
    fp2sqr_mont_weights(Q->Z, Q->Z, hamming_weights[2], hamming_weights[3]);        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2

    // 4-5        // old 16-17
    fp2sqr_mont_weights(Q->X, Q->X, hamming_weights[4], hamming_weights[5]);        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2

    // 6-7        // old 18-19
    fp2mul_mont_weights(Q->Z, xPQ, Q->Z, hamming_weights[6], hamming_weights[7]);   // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}



void compute_hamming_weights(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, int start_step, int end_step, int hamming_weights[2][10][3][14], FILE* readpts, FILE* write_0, FILE *write_1)
{
    // Bob's ephemeral shared secret computation
    // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
    // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1]. 
    //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
    // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t G = {0}, G0 = {0}, G2 = {0};
    point_proj_t G_1 = {0}, G0_1 = {0}, G2_1 = {0};
    f2elm_t PKB[3];
    f2elm_t A = {0}, A24 = {0};
    digit_t mask, SecretKeyB[NWORDS_ORDER] = {0};
    int i, bit, swap, prevbit = 0;


    // Initialize images of Alice's basis
    fp2_decode(PublicKeyA,                          PKB[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES,      PKB[1]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES,    PKB[2]);

    // Initializing constant A
    get_A(PKB[0], PKB[1], PKB[2], A);
    fpcopy((digit_t*)&Montgomery_one, A24[0]);
    mp2_add(A24, A24, A24);
    mp2_add(A, A24, A24);
    fp2div2(A24, A24);  
    fp2div2(A24, A24);

    //Read Bob's sk (obtained so far) in usable format
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);

    // READ STARTING POINT (If 0 STARTING POINT = PUBLIC KEY ALICE)
    // G0 = Q, G2 = Q-P, G = P // R0 R2 R
    assert( fread((point_proj_t *)G0, sizeof(point_proj_t), 1, readpts) != 0 );
    assert( fread((point_proj_t *)G2, sizeof(point_proj_t), 1, readpts) != 0 );
    assert( fread((point_proj_t *)G , sizeof(point_proj_t), 1, readpts) != 0 );


    // We swap G and G2 if and only if (BIT ^ prevbit == 1), BIT being the bit we try to guess
    if (start_step > 0) 
        prevbit = (SecretKeyB[(start_step-1) >> LOG2RADIX] >> ((start_step-1) & (RADIX-1))) & 1;

    // Main loop
    for (i = start_step; i < end_step; i++) {
        bit = (SecretKeyB[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(G, G2, mask);
        xDBLADD(G0, G2, G->X, A24);
        fp2mul_mont(G2->X, G->Z, G2->X);
    }

    
    // NOW WE DO ANOTHER STEP OF THE LOOP WITH i = end_step and bit = 0
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(G, G2, mask);    

    //CASE IF BIT = 0 -> xDBLADD(G0, G2, G)
    //CASE IF BIT = 1 -> xDBLADD(G0, G, G2)
    //IN OTHER WORDS WE SWAP G with G2
    //WE DO THIS BY COPYING POINTS G0, G2, AND G TO NEW VARIABLES TO BE USED
    //IN THE BIT = 1 CASE, THE SWAPPING BEING DONE AUTOMATICALLY IN THE COPY
    fp2copy(G->X, G2_1->X);
    fp2copy(G->Z, G2_1->Z);    

    fp2copy(G2->X, G_1->X);
    fp2copy(G2->Z, G_1->Z);

    fp2copy(G0->X, G0_1->X);
    fp2copy(G0->Z, G0_1->Z);


    //CASE IF BIT = 0
    xDBLADD_weights(G0, G2, G->X, A24, hamming_weights[0]);
    fp2mul_mont_weights(G2->X, G->Z, G2->X, hamming_weights[0][8], hamming_weights[0][9]); 
    fwrite((point_proj_t *)G0, sizeof(point_proj_t), 1, write_0);
    fwrite((point_proj_t *)G2, sizeof(point_proj_t), 1, write_0);
    fwrite((point_proj_t *)G , sizeof(point_proj_t), 1, write_0);


    //CASE IF BIT = 1
    xDBLADD_weights(G0_1, G2_1, G_1->X, A24, hamming_weights[1]);
    fp2mul_mont_weights(G2_1->X, G_1->Z, G2_1->X, hamming_weights[1][8], hamming_weights[1][9]); 
    fwrite((point_proj_t *)G0_1, sizeof(point_proj_t), 1, write_1);
    fwrite((point_proj_t *)G2_1, sizeof(point_proj_t), 1, write_1);
    fwrite((point_proj_t *)G_1 , sizeof(point_proj_t), 1, write_1);

}





int hamming_weight_algorithm(int argc, char* argv[])
{

    // OALICE_BITS               216  
    // OBOB_BITS                 218     
    // MSG_BYTES                  16       
    // SECRETKEY_A_BYTES          27    ((OALICE_BITS + 7) / 8)
    // SECRETKEY_B_BYTES          28    ((OBOB_BITS - 1 + 7) / 8)
    // CRYPTO_PUBLICKEYBYTES     330
    // CRYPTO_SECRETKEYBYTES     374    // MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes
    // CRYPTO_BYTES               16
    // CRYPTO_CIPHERTEXTBYTES    346    // CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes  

    //bob_sk = (randomness || secret key as a value up to 3^... || public key of corresponding secret key)
    //bob_pk = (public key of corresponding secret key) ---> these are three x coordinates encoded with fp2encode (Px || Qx || PQx)

    //alice_sk = HASH_G(randomness temp || bob_pk)   with some padding in the last byte
    //jinvariant = (jinvariant of common curve encoded with fp2encode)
    //h = hash (jinvariant)
    //alice_pk = (public key of corresponding secret key || (randomness temp XOR h) )  ---> these are three x coordinates encoded with fp2encode (Px || Qx || PQx), and then xor some bits
    //shared_secret = HASH( randomness temp || alice_pk )

    // Bob reads alice_pk, takes the public key and computes the common j invariant and the hash
    // from that he obtains the randomness by xoring and generates ephemeral_sk by doing the same procedure with the randomness and his static_pk
    // If this key gives the same alices public key then Bob can trust alice, and he chooses the shared secret to be the hash of randomness temp concatenated with the alices pk
    // Otherwise Bob uses his own randomness instead of the one obtained from Alices public key.


    /*
    Input 
    argv[1] = Number of Alice's (ephermeral) keys 
    argv[2],...,argv[number_keys+1] -> Indices of keys that are used
    argv[number_keys + 2] = Number of bits of Bob's secret key that is already recovered
    argv[number_keys + 3] = String bits of Bob's secret key that is already known 
    */
    

    // Modified get_options(argc, argv)
    long number_keys = strtol(argv[1], NULL, 10);
    int KEYS[number_keys];

    for(int i = 0; i < number_keys; i++)
        KEYS[i] = strtol(argv[i + 2], NULL, 10);


    int start_step = strtol(argv[number_keys + 2], NULL, 10);
    int end_step = strtol(argv[number_keys + 3], NULL, 10);
    char BITS[end_step]; 
    char* keys_directory;

    keys_directory = argv[number_keys + 4];

    for(int j = 0; j < end_step; j++)
        BITS[j] = argv[number_keys + 5][j];


    // DEFINE BOB'S PUBLIC KEY AND SECRET KEY
    unsigned char bob_sk[CRYPTO_SECRETKEYBYTES];
    unsigned char bob_pk[CRYPTO_PUBLICKEYBYTES];
    char bob_pk_path[200];
    sprintf(bob_pk_path, "%sbob_keys/bob_pk", keys_directory);

    // READ BOB'S PUBLIC KEY
    read_bob_pk(bob_pk_path, bob_pk);

    // ZERO OUT BOB'S SECRET KEY
    create_empty_bob_sk(bob_pk, bob_sk);

    // FILL BOB's SECRET KEY
    for(int i = 0; i < end_step; i++)
        if(BITS[i] == '1')
            bob_sk[(MSG_BYTES*8 + i) >> 3] |= ( 1 << ((MSG_BYTES*8 + i) & (8-1)) );


    //DEFINE ALICE'S PUBLIC KEY
    unsigned char alice_pk[CRYPTO_CIPHERTEXTBYTES];
    char alice_pk_path[200];

    //We are trying to guess the secret key of Bob
    //At each step we are guessing one bit of Bob's secret key by computing a secret key of Alice which makes Bob's computation leak a bit of ouptut
    int hamming_weights[2][10][3][14];
    char weights_file_name[200];
    char points_paths[3][200];


    sprintf(weights_file_name, "./weights/bit_%03d", end_step);
    FILE *weights_file  = fopen(weights_file_name, "wb");
    assert((weights_file != NULL) && "Error writing weights file");

    sprintf(points_paths[0], "points/pt_%03d", start_step);
    FILE *readpts = fopen(points_paths[0], "rb");
    assert((readpts != NULL) && "Error reading point file");

    sprintf(points_paths[1], "points/pt_%03d_0", end_step + 1);
    FILE *write_0 = fopen(points_paths[1], "wb");
    assert((write_0 != NULL) && "Error writing point file");

    sprintf(points_paths[2], "points/pt_%03d_1", end_step + 1);
    FILE *write_1 = fopen(points_paths[2], "wb");
    assert((write_1 != NULL) && "Error writing point file");

    //Go through all the keys, and for each one compute the hamming weights of the output of the ADDS function for some calls
    for(int alice_key_counter = 0; alice_key_counter < number_keys; alice_key_counter++)
    {
        point_proj_t G = {0}, G0 = {0}, G2 = {0}; //Guess points in case bit=0 or bit=1
        f2elm_t const_A = {0};

        memset(alice_pk, 0x00, CRYPTO_CIPHERTEXTBYTES);
        sprintf(alice_pk_path, "%salice_keys/alice_pk_%05d", keys_directory, KEYS[alice_key_counter]);
        read_alice_pk(alice_pk_path, alice_pk);

        compute_hamming_weights(bob_sk + MSG_BYTES, alice_pk, start_step, end_step, hamming_weights, readpts, write_0, write_1);

        // PRINT WEIGHTS TO FILE
        fwrite(hamming_weights, sizeof(int), 2*10*3*14, weights_file);
    }

    fclose(weights_file);
    fclose(readpts);
    fclose(write_0);
    fclose(write_1);




  return 0;
}




int make_starting_points(int argc, char* argv[])
{
    //DEFINE ALICE'S PUBLIC KEY
    unsigned char alice_pk[CRYPTO_CIPHERTEXTBYTES];
    char alice_pk_path[200];
    char points_path[200];
    char* keys_directory;
    int number_keys = argc - 2;
    int KEYS[number_keys];

    for(int i = 0 ; i < number_keys; i++)
        KEYS[i] = strtol(argv[i + 1], NULL, 10);

    keys_directory = argv[number_keys + 1];

    sprintf(points_path, "points/pt_%03d", 0);


    FILE *fp = fopen(points_path, "wb");
    assert((fp != NULL) && "Error writing point file");

    for(int alice_key_counter = 0; alice_key_counter < number_keys; alice_key_counter++)
    {
        point_proj_t R = {0}, R0 = {0}, R2 = {0};
        f2elm_t xP, xQ, xPQ;

        memset(alice_pk, 0x00, CRYPTO_CIPHERTEXTBYTES);
        sprintf(alice_pk_path, "%salice_keys/alice_pk_%05d", keys_directory, KEYS[alice_key_counter]);
        read_alice_pk(alice_pk_path, alice_pk);
        
        fp2_decode(alice_pk,                        xP);
        fp2_decode(alice_pk + FP2_ENCODED_BYTES,    xQ);
        fp2_decode(alice_pk + 2*FP2_ENCODED_BYTES,  xPQ);
    
        // Initializing points p q pq
        fp2copy(xQ, R0->X);
        fpcopy((digit_t*)&Montgomery_one, (digit_t*)R0->Z);
        fpzero((digit_t*)(R0->Z)[1]);

        fp2copy(xPQ, R2->X);
        fpcopy((digit_t*)&Montgomery_one, (digit_t*)R2->Z);
        fpzero((digit_t*)(R2->Z)[1]);

        fp2copy(xP, R->X);
        fpcopy((digit_t*)&Montgomery_one, (digit_t*)R->Z);
        fpzero((digit_t*)(R->Z)[1]);

        fwrite((point_proj_t *)R0, sizeof(point_proj_t), 1, fp);
        fwrite((point_proj_t *)R2, sizeof(point_proj_t), 1, fp);
        fwrite((point_proj_t *)R, sizeof(point_proj_t), 1, fp);

    }

    fclose(fp);

    return 0;
}







int bob_key_generation(int argc, char* argv[])
{
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    char static_pk_path[200] = "bob_keys/bob_pk";
    char static_sk_path[200] = "bob_keys/bob_sk";
    char only_sk_path[200] = "bob_keys/only_bob_sk";



    // Create Bob's keys
    randombytes(sk, MSG_BYTES);

    if(argc <= 1)
        random_mod_order_B(sk + MSG_BYTES);

    else
        for(int i = 0; i < SECRETKEY_B_BYTES; i++)
            sk[MSG_BYTES + i] = ((uint8_t)(FromHex(argv[1][2*i])<<4)|(uint8_t)FromHex(argv[1][2*i+1]));
    
    EphemeralKeyGeneration_B(sk + MSG_BYTES, pk);


    // Write Bob's keys
    write_bob_pk(static_pk_path, pk);
    memcpy(&sk[MSG_BYTES + SECRETKEY_B_BYTES], pk, CRYPTO_PUBLICKEYBYTES);
    write_bob_sk(static_sk_path, sk);
    write_only_bob_sk(only_sk_path, &sk[MSG_BYTES]);

    return 0;
}



int alice_key_generation(int argc, char* argv[])
{
    unsigned char temp[CRYPTO_CIPHERTEXTBYTES+MSG_BYTES];
    unsigned char ephemeral_sk[SECRETKEY_A_BYTES];
    unsigned char ephemeral_pk[CRYPTO_CIPHERTEXTBYTES];
    unsigned char shared_secret[CRYPTO_BYTES] = {0};
    unsigned char jinvariant[FP2_ENCODED_BYTES];
    unsigned char h[MSG_BYTES];
    char ephemeral_pk_path[200];

    unsigned char static_pk[CRYPTO_PUBLICKEYBYTES];
    char static_pk_path[200];
    sprintf(static_pk_path, "bob_keys/bob_pk");
    int NUMBER_KEYS = 1;
    if (argc >= 2)
        NUMBER_KEYS = strtol(argv[1], NULL, 10);

    // READ BOB'S PUBLIC KEY
    read_bob_pk(static_pk_path, static_pk);

    for(int key_counter = 0; key_counter < NUMBER_KEYS; key_counter++)
        {

            memset(temp, 0x00, CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);
            memset(ephemeral_sk, 0x00, SECRETKEY_A_BYTES);
            memset(ephemeral_pk, 0x00, CRYPTO_CIPHERTEXTBYTES);
            memset(shared_secret, 0x00, CRYPTO_BYTES);
            memset(jinvariant, 0x00, FP2_ENCODED_BYTES);
            memset(h, 0x00, MSG_BYTES);


            // GENERATE A SECRET AND PUBLIC EPHEMERAL KEY
            generate_alice_sk_pk(static_pk, temp, ephemeral_sk, ephemeral_pk, shared_secret, jinvariant, h);
            sprintf(ephemeral_pk_path, "alice_keys/alice_pk_%05d", key_counter);
            write_alice_pk(ephemeral_pk_path, ephemeral_pk);

        }

    return 0;
}


/*
int alice_key_generation(int argc, char* argv[])
{
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
    unsigned char ss[CRYPTO_BYTES]; 
    char path_alice_pk[200] = "static_keys/alice_pk";
    char path_bob_pk[200] ="static_keys/bob_pk";

    read_bob_pk(path_bob_pk, pk);
    crypto_kem_enc(ct, ss, pk);
    write_alice_pk(path_alice_pk, ct);
    return 0;
}

void complain(char *fmt, ...)
{
    va_list arglist;
    va_start(arglist, fmt);
    vfprintf(stderr, fmt, arglist);
    exit(1);
}

static void usage()
{
    fprintf(stderr, "key_generator -n number_keys [ options ]\n");
    fprintf(stderr, "  -n : number of keys to create\n");
    fprintf(stderr, "  -b : bits already known\n");
    fprintf(stderr, "  -h help\n");
    exit(1);
}

static void get_options(int argc, char **argv)
{
    char c;
    NUMBER_KEYS = 0;
    memset(BITS,0x00,512);
    
    while ((c = getopt(argc, argv, "n:h:b")) != (char)(-1))
    {
        switch (c)
        {
            case 'n':
                if (sscanf(optarg, "%d", &NUMBER_KEYS) != 1)
                    complain("Bad argument to -n!\n");
                break;
            case 'h':
                usage();
                break;
            case 'b':
                if (sscanf(optarg, "%s", BITS) != 1) //to do
                    complain("Bad argument to -b!\n");
                break;
            default:
                fprintf(stderr, "Bad option %c\n", (char)c);
                usage();
        }
    }

    if (NUMBER_KEYS == 0)
    {
        fprintf(stderr, "argument -n number_keys is necessary\n");
        usage();
    }
}
*/




