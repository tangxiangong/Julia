import numpy as np
from numpy.random import randn


def bm(T, tau):
    n = int(T/tau) + 1
    t = np.linspace(0, T, n)
    noise = np.sqrt(tau) * randn(n-1)
    noise = np.insert(noise, 0, 0)
    x = np.cumsum(noise)
    return t, x


if __name__ == "__main__":
    import time
    s = time.time()
    bm(1000, 0.01)
    print(f"time = {(time.time()-s) * 1000} ms")
