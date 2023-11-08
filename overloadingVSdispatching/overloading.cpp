#include<iostream>
#include<string>
using namespace std;

class NDArray { public: string name; };
class Vector : public NDArray {};
class Matrix : public NDArray {};

string multiply(NDArray a, NDArray b) { return "不确定"; }
string multiply(Vector a, Vector b) { return "向量点乘"; }
string multiply(Matrix a, Vector b) { return "矩阵乘向量"; }
string multiply(Vector a, Matrix b) { return "维数不一致"; }
string multiply(Matrix a, Matrix b) { return "矩阵乘积"; }

void result(NDArray a, NDArray b) {cout << a.name << "和" << b.name << "相乘的结果是" << multiply(a, b) << endl;}

int main() {
    Vector v {"向量"}; Matrix m {"矩阵"};
    result(v, v); result(m, v); result(v, m); result(m, m);
}