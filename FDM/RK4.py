import numpy as np
from typing import Callable


class ODE(object):
    def __init__(self, T: float, f:Callable[[float, float], float], y0: float):
        self.T = T
        self.f = f
        self.y0 = y0

class Solver(object):
    def __init__(self, tau: float):
        self.tau = tau

class RKSolver(Solver):
    def __init__(self, tau: float):
        super().__init__(tau)

class RK4(RKSolver):
    def __init__(self, tau: float):
        super().__init__(tau)

def solve(ode: ODE, solver: RK4):
    T, f, y0 = ode.T, ode.f, ode.y0
    tau = solver.tau
    t = np.linspace(0, T, int(T/tau)+1)
    y = np.zeros_like(t)
    y[0] = y0
    for n in range(1, len(t)):
        k1 = tau * f(t[n-1], y[n-1])
        k2 = tau * f(t[n-1]+tau/2, y[n-1]+k1/2)
        k3 = tau * f(t[n-1]+tau/2, y[n-1]+k2/2)
        k4 = tau * f(t[n], y[n-1]+k3)
        y[n] = y[n-1] + (k1+2*k2+2*k3+k4)/6
    return t, y

if __name__ == "__main__":
    def f_(t: float, y: float, T):
        return y + np.exp(-t/T) * t * np.sin(t)
    T, y0 = 10.0, 1.0 
    tau = 1/256
    rk4solver = RK4(tau)
    f = lambda t, x: f_(t, x, T)
    ode = ODE(T, f, y0)
    from time import time
    start = time()
    solve(ode, rk4solver)
    print((time()-start)*1000)
    #  24ms