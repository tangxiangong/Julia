#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


double f(double, double, double);
void solve(double (*)(double, double, double), double, double, double, double*, double*, int);

double f(double t, double y, double T) {
    double result;
    result = y + exp(-t/T) * t * sin(t);
    return result;
}

void solve(double (*f)(double, double, double), double T, double y0, double tau, double* t, double* y, int n) {
    t[0] = 0.0;
    y[0] = y0;
    double k1, k2, k3, k4;
    for(int k=1; k<=n; k++) {
        t[k] = t[k-1] + tau;
        k1 = tau * f(t[k-1], y[k-1], T);
        k2 = tau * f(t[k-1]+tau/2, y[k-1]+k1/2, T);
        k3 = tau * f(t[k-1]+tau/2, y[k-1]+k2/2, T);
        k4 = tau * f(t[k], y[k-1]+k3, T);
        y[k] = y[k-1] + (k1+2*k2+2*k3+k4)/6;
    }
}

int main() {
    clock_t start, end;
    double T, tau, y0;
    T = 10.0; tau = 1/256; y0 = 1.0;
    int n = 2561;
    
    double* t = (double*)malloc(sizeof(double)*n);
    double* y = (double*)malloc(sizeof(double)*n);
    start = clock();
    solve(f, T, y0, tau, t, y, n);
    end = clock();
    
    printf("time=%f\n", (double)1000*(end-start)/CLOCKS_PER_SEC);
    return 0;
}