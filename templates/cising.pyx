import cython
import numpy
cimport numpy

# declare the interface to the C code
cdef extern void c_init()
cdef extern void c_random_list(long N, double *rs)

# Cython interface to C function
def random_list(numpy.ndarray[double, ndim=1, mode='c'] rs not None):
    # pass the array to the C function
    cdef long N = rs.shape[0]
    c_random_list(N, &rs[0])
    return rs

# call c_init() when the module is loaded to init the gsl RNG
c_init()
