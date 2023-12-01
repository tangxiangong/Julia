import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg


def Euler(T, u0, f, tau, M):
    t = np.arange(0, T+tau, step=tau)
    # u = np.zeros((t.size, u0.size))
    x = u0
    A = sp.eye(u0.size) - tau * M
    b = np.zeros_like(u0)
    for n in range(1, t.size):
        b = tau * f(t[n]) + x
        x, _ = cg(A, b)
    return t, x

h = 1/128
x = np.arange(h, 1, step=h)
u0 = np.sin(np.pi * x)
f = lambda t: (np.pi**2-1)*np.exp(-t)*np.sin(np.pi * x)
diagonals = [-2*np.ones(x.size), np.ones(x.size-1), np.ones(x.size-1)]
M = sp.diags(diagonals=diagonals, offsets=[0, 1, -1])/h**2

# %timeit Euler(1, u0, f, 1/256, M)  22 ms
