#include <iostream>
using std::cout;
using std::endl;

template <typename T>
void swap(T& a, T& b) {
    T c = a;
    a = b;
    b = c;
}

int main() {
    double a = 1.2;
    double b = 2.3;
    swap(a, b);
    cout << a << endl;
    cout << b << endl;
    return 0; 
}