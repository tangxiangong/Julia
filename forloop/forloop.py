import numpy as np
from time import time
from numba import njit

@njit
def seriel_numba(x):
    y = np.zeros_like(x)
    for k in range(len(x)):
        y[k] = np.sin(x[k])
    return y

n = 100_000
x = np.random.rand(n)

s = time()
list(map(np.sin, x))
print(f"Python for 循环\n {1000*(time()-s)} ms")

seriel_numba(np.random.rand(10))    # JIT 预热
s = time()
seriel_numba(x)
print(f"Numba for 循环\n {1000*(time()-s)} ms")

s = time()
np.sin(x)
print(f"NumPy 向量化\n {1000*(time()-s)} ms")


