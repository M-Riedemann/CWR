#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// zu untersuchende Funktion
double p1(double y)
{
    return y*y*y - y/2;
}


// num. Integration (Trapez)
double integrate_trapez(double left, double right, int N, double integrand(double)) {
    double sum = 0;
    double delta_y = (right - left) / N;
    for (int k = 0; k < N; k++) {
        sum += (integrand(left + k*delta_y) + integrand(left + (k+1)*delta_y)) * delta_y / 2;
    }
    return sum;
}

// num. Integration (Rechteck)
double integrate(double left, double right, int N, double integrand(double)) {
    double sum = 0;
    double delta_y = (right - left) / N;
    for(int k = 0; k < N; ++k)  {
        double y = left + k * delta_y;
        double f = integrand(y);
        double A = f * delta_y;
        sum += A;
    }
    return sum;
}

// analytische Lösung
double genau(double x)
{
    return (x*x*x*x)/4 -x*x/4;
}


int main(void) {
    double left = 0.0;
    double right = 2.0;
    for(int k = 1; k < 8; ++k) {
        int N = pow(10,k);
        double result = integrate(left, right, N, p1);
        double resultg = genau(right);
        //printf("Integralwert: %f\n", result);
        //printf("genauer Integralwert: %f\n", resultg);
        double diff = fabs(result - resultg);
        printf("absolute Differenz: %f\n", diff);
        FILE* myFile = fopen("data.csv", "a");
        fprintf(myFile, "%g, %d\n", diff, N);
        fclose(myFile);
        }
    return EXIT_SUCCESS;
}

