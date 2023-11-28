import numpy as np
from numpy.random import choice

class Mod(object):
    def __init__(self, buffer, mod):
        assert(isinstance(buffer, int) and isinstance(mod, int))
        self.buffer = buffer
        self.mod = mod
        self.value = buffer % mod
    def __repr__(self):
        return f"{self.value}(mod {self.mod})"
    def __add__(self, rhs):
        assert(self.mod == rhs.mod)
        return Mod(self.value+rhs.value, self.mod)
    def __neg__(self):
        return Mod(self.mod-self.value, self.mod)
    def __sub__(self, rhs):
        return self + (-rhs)
    def __mul__(self, rhs):
        assert(self.mod == rhs.mod)
        return Mod(self.value*rhs.value, self.mod)

def rand(mod, *args):
    nums = [Mod(n, mod) for n in range(mod)]
    return choice(nums, size=args)