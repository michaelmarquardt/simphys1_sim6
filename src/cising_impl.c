#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>

gsl_rng *r = NULL;

/* Initialize the random number generator. This function is called from cising.pyx. */
void c_init() {
  r = gsl_rng_alloc(gsl_rng_taus);
}

/* Generate a random integer between 0 and N-1. */
long gsl_randint(long N) {
  return gsl_rng_uniform_int(r, N);
}

/* Generate a random float between 0.0 and 1.0. */
double gsl_rand() {
  return gsl_rng_uniform(r);
}

/* This is an example for a function that is exported to Cython and
   uses the GSL RNG. */
void c_random_list(long N, double* rs) {
  long i;
  for (i = 0; i < N; i++)
    rs[i] = gsl_rand();
}

/* Spin flip*/
void c_spin_flip(double beta, long L, double* E, double* mu, int* sigma) {
    long V = L*L;
    for (long step=0; step<V; step++) {
        // k = i + j*L   =>   i = k%L   j = k//L
        long k = gsl_randint(V);
        // Test deltaE before flipping spin
        int deltaE = 2*sigma[k]*(sigma[k/L*L+(k+1+L)%L]+sigma[k/L*L+(k-1+L)%L]+sigma[(k+L+V)%V]+sigma[(k-L+V)%V]);
        if (gsl_rand() < exp(-beta*deltaE)) {
            sigma[k] *= -1;     // spin flip
            *E += deltaE;
            *mu += 2*sigma[k];
        }
    }
}


