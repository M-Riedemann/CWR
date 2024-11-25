

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "../my_numerics.h"




// Struktur des inneren Zustandes eines Zufallsgenerators
struct randu_random_t {
    uint32_t state;
};


// Zufallsgenerator des Typs RANDNU
uint32_t randu_random_r(struct randu_random_t *rng) {
    uint32_t a = (1u << 16) + 3;
    uint32_t b = 0;
    uint32_t m =  (1u << 31);
    rng->state = (a * (uint64_t)rng->state + b) % m;
    return rng->state;
}


int main(void) {
    struct randu_random_t PRN_Randu;
    PRN_Randu.state = time(NULL);

    gsl_rng *PRN_MT = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(PRN_MT, time(NULL));

    FILE* randu_file = fopen("A25_Randu.csv", "w");
    FILE* mt_file = fopen("A25_Mersenne_Twister.csv", "w");

    for (int i = 0; i < 2500; i++) {
        double randu1 = (double) randu_random_r(&PRN_Randu) / (1u << 31);
        double randu2 = (double) randu_random_r(&PRN_Randu) / (1u << 31);
        double randu3 = (double) randu_random_r(&PRN_Randu) / (1u << 31);

        double mt1 = gsl_rng_uniform(PRN_MT);
        double mt2 = gsl_rng_uniform(PRN_MT);
        double mt3 = gsl_rng_uniform(PRN_MT);

        fprintf(randu_file, "%g, %g, %g\n", randu1, randu2, randu3);
        fprintf(mt_file, "%g, %g, %g\n", mt1, mt2, mt3);
    }
    fclose(randu_file);
    fclose(mt_file);
    return 0;
}
