/*
假设已有一个圆类和一个矩形类，都实现了计算其周长的方法。
如何在不修改原有代码的前提条件下，新增一个三角形类、实现其计算周长的算法，并且对这三个类型新增计算面积的方法。
*/

// 已有代码
#include <cmath>
#include <cassert>

const double PI = 3.1415926;

class Shape { public: virtual double perimeter() = 0; };

template <class T>
class Circle:public Shape {
private:
    T m_radius;
public:
    Circle() {};
    Circle(const T&);
    T get_radius();
    double perimeter();
};

template <class T>
Circle<T>::Circle(const T& radius) { assert(radius > 0); this->m_radius = radius; }

template <class T>
T Circle<T>::get_radius() { return this->m_radius; }

template <class T>
double Circle<T>::perimeter() { return 2 * PI * this->m_radius; }

template<class T1, class T2>
class Rectangle:public Shape {
private:
    T1 m_height;
    T2 m_width;
public:
    Rectangle() {};
    Rectangle(const T1&, const T2&);
    T1 get_height();
    T2 get_width();
    double perimeter();
};

template <class T1, class T2> 
Rectangle<T1, T2>::Rectangle(const T1& height, const T2& width) {
    assert(height > 0 && width > 0);
    this->m_height = height;
    this->m_width = width;
}

template <class T1, class T2>
T1 Rectangle<T1, T2>::get_height() { return this->m_height; }

template <class T1, class T2>
T2 Rectangle<T1, T2>::get_width() { return this->m_width; }

template <class T1, class T2>
double Rectangle<T1, T2>::perimeter() { return (double)(this->m_height + this->m_width) * 2; }

// 新增代码
template <class T1, class T2, class T3>
class Triangle:public Shape {
private:
    T1 m_a;
    T2 m_b;
    T3 m_c;
public:
    Triangle() {};
    Triangle(const T1&, const T2&, const T3&);
    auto get_edges();
    double perimeter();
};

template <class T1, class T2, class T3>
Triangle<T1, T2, T3>::Triangle(const T1& a, const T2& b, const T3& c) {
    assert(a > 0 && b > 0 && c > 0);
    assert(a + b > c && a + c > b && b + c > a);
    this->m_a = a;
    this->m_b = b;
    this->m_c = c;
}

template <class T1, class T2, class T3>
double Triangle<T1, T2, T3>::perimeter() {return (double)(this->m_a + this->m_b + this->m_c); }

template <class T1, class T2, class T3>
auto Triangle<T1, T2, T3>::get_edges() {
    using T = decltype(this->m_a + this->m_b + this->m_c);
    T* edges = new T[3];
    edges[0] = this->m_a;
    edges[1] = this->m_b;
    edges[2] = this->m_c;
    return edges;
}

template <class T>
double area(Circle<T> c) { 
    T r = c.get_radius();  
    double a = PI * r * r;
    return a; 
}

template <class T1, class T2>
auto area(Rectangle<T1, T2> r) { 
    T1 height = r.get_height(); 
    T2 width = r.get_width();    
    auto a = height * width;
    return a;
}

template <class T1, class T2, class T3>
double area(Triangle<T1, T2, T3> t) {
    auto edges = t.get_edges();
    double s = t.perimeter()/2;
    double a = sqrt(s * (s - edges[0]) * (s - edges[1]) * (s - edges[2]));
    delete[] edges;
    return a;
}