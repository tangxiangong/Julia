import numpy as np
from numba import njit
from typing import Callable


def solve(T: float, f: Callable[[float, float, float], float], y0, tau):
    t = np.linspace(0, T, int(T/tau)+1)
    y = np.zeros_like(t)
    y[0] = y0
    for n in range(1, len(t)):
        k1 = tau * f(t[n-1], y[n-1], T)
        k2 = tau * f(t[n-1]+tau/2, y[n-1]+k1/2, T)
        k3 = tau * f(t[n-1]+tau/2, y[n-1]+k2/2, T)
        k4 = tau * f(t[n], y[n-1]+k3, T)
        y[n] = y[n-1] + (k1+2*k2+2*k3+k4)/6
    return t, y

if __name__ == "__main__":
    def f(t: float, y: float, T):
        return y + np.exp(-t/T) * t * np.sin(t)
    T, y0 = 10.0, 1.0 
    tau = 1/256
    from time import time
    start = time()
    solve(T, f, y0, tau)
    print((time()-start)*1000)
    #  24ms