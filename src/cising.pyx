import cython
import numpy
cimport numpy

# declare the interface to the C code
cdef extern void c_init()
cdef extern void c_random_list(long N, double *rs)
cdef extern void c_spin_flip(double beta, long L, double* E, double* mu, int* sigma)

# Cython interface to C function
def random_list(numpy.ndarray[double, ndim=1, mode='c'] rs not None):
    # pass the array to the C function
    cdef long N = rs.shape[0]
    c_random_list(N, &rs[0])
    return rs

def spin_flip(beta, E_in, mu_in, numpy.ndarray[int, ndim=2, mode='c'] sigma not None):
    cdef long L = sigma.shape[0]
    cdef double E = E_in
    cdef double mu = mu_in
    c_spin_flip(beta, L, &E, &mu, &sigma[0,0])
    return E, mu

# call c_init() when the module is loaded to init the gsl RNG
c_init()
