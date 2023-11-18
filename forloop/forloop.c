#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double serial(double* x, int n) {
    double result = 0.0;
    for(int k=0; k<n; k++)
        result += sin(x[k]) * sin(x[k]) + cos(x[k]) * cos(x[k]);
    result /= n;
    return result;
}


int main() {
    int n = 1000000;
    double* x = (double*)malloc(sizeof(double)*n);
    srand((unsigned)time(NULL));
    for(int k = 0; k < n; k++)
        x[k] = rand()/(double)RAND_MAX;    
    size_t start = clock();
    double result = serial(x, n);
    size_t end = clock();
    printf("C for 循环\n  结果 %f\n  时间 %f ms\n", result, (double)1000*(end-start)/CLOCKS_PER_SEC);
    return 0;
}