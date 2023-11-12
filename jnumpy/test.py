# from jl_fftw import fft as j_fft
from jl_fftw import mc_pi
# import jnumpy as jnp
# import numpy as np
from time import time

n = 100000
# x = np.random.rand(n)
js = time()
pi = mc_pi(n)
print("time:", time()-js)
# ps = time()
# np_f = np.fft.fft(x)
# print("time:", time()-ps)


