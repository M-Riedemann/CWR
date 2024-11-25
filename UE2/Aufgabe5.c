#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Funktion
double f(double x){
    return x*x*(x-1);
}

//einfache Ableitung
double diff(double x, double delta, double func(double)){
    return (func(x+delta)-func(x))/delta;
}

double centered_diff(double x, double delta, double func(double)) {
    return ( func(x + delta) - func(x - delta) ) / ( 2 * delta);
}




int main(void) {
    double x = 1, delta, f_diff_numerical, f_centered_diff_numerical;
    
    FILE* file_diff_f = fopen("A05_Differenzierung.csv", "w");
    fprintf(file_diff_f, "Delta, Abweichung1, Abweichung2\n");
    for (int exp = 0; exp < 17; exp++) {
         delta = pow(10, -exp);
         f_diff_numerical = diff(x, delta, f);
         f_centered_diff_numerical = centered_diff(x, delta, f);
         fprintf(file_diff_f, "%.16f, %.10f, %.12f\n", delta, fabs(f_diff_numerical - 1), fabs(f_centered_diff_numerical - 1));}
    fclose(file_diff_f);

    return EXIT_SUCCESS;
}

