/*
g++ main.cpp -std=c++17 -I./ -o main && ./main
*/

#include <iostream>
#include <shape.hpp>

using std::cout;
using std::endl;

int main() {
    Circle c {2};
    cout << "圆周长 " << c.perimeter() << endl;
    cout << "圆面积 " << area(c) << endl;
    Rectangle r {1.4, 3};
    cout << "矩形周长 " << r.perimeter() << endl;
    cout << "矩形面积 " << area(r) << endl;
    Triangle t {3, 4, 5.0};
    cout << "三角形周长 " << t.perimeter() << endl;
    cout << "三角形面积 " << area(t) << endl;
    return 0;
}

