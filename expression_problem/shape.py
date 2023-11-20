
# 已有代码
from numpy import pi
from abc import abstractmethod, ABCMeta


class Shape(metaclass=ABCMeta):
    @abstractmethod
    def perimeter(self): pass

class Circle(Shape):
    def __init__(self, radius):
        assert(radius > 0)
        self.__radius = radius

    def get_radius(self):
        return self.__radius

    def perimeter(self):
        return 2 * self.__radius * pi
    
class Rectangle(Shape):
    def __init__(self, height, width):
        assert(height > 0 and width > 0)
        self.__height = height
        self.__width = width

    def get_height_width(self):
        return self.__height, self.__width

    def perimeter(self):
        return 2 * (self.__height + self.__width)
    
# 新增代码
from functools import singledispatch
from math import sqrt

class Triangle(Shape):
    def __init__(self, a, b, c):
        assert(a > 0 and b > 0 and c > 0)
        assert(a + b > c and a + c > b and b + c > a)
        self.__a = a
        self.__b = b
        self.__c = c
    
    def get_edges(self):
        return self.__a, self.__b, self.__c

    def perimeter(self):
        return self.__a + self.__b + self.__c
    
@singledispatch
def area(Shape):
    print("算法未实现")

@area.register
def _(c: Circle):
    r = c.get_radius()
    return pi * r * r

@area.register
def _(r: Rectangle):
    h, w = r.get_height_width()
    return h * w

@area.register
def _(t: Triangle):
    a, b, c = t.get_edges()
    s = t.perimeter()/2
    return sqrt(s * (s - a) * (s - b) * (s - c))

if __name__ == "__main__":
    c = Circle(2)
    r = Rectangle(1.4, 3)
    t = Triangle(3, 4, 5)
    print(f"圆周长 {c.perimeter()}")
    print(f"圆面积 {area(c)}")
    print(f"矩形周长 {r.perimeter()}")
    print(f"矩形面积 {area(r)}")
    print(f"三角形周长 {t.perimeter()}")
    print(f"三角形面积 {area(t)}")