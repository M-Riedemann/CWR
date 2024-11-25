#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// zu untersuchende Funktion f
double erf(double y) {
    return 2 * exp(-pow(y, 2)) / sqrt(M_PI);
}

// num. Integration (Mittelpunktsregel)
double erf_midpoint(double x, double delta_x) {
    int N = fabs(x)/delta_x;
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += erf( (i-0.5)*delta_x) * delta_x;
    }
    return sum;
}

// num. Integration von Cosh (Mittelpunktsregel)
double cosh_midpoint(double x, double delta_x) {
    int N = fabs(x)/delta_x;
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += cosh( (i-0.5)*delta_x) * delta_x;
    }
    return sum;
}


// num. Integration (Simpsonregel)
double erf_simpson(double x, double delta_x) {
    int N = fabs(x) / delta_x;
    double sum = 0, x_i, m_i, x_ii;
    for (int i = 0; i < N; i++) {
        x_i = delta_x * i;
        m_i = delta_x * (i+0.5);
        x_ii = delta_x * (i+1);
        sum += ( erf(x_i) + 4*erf(m_i) + erf(x_ii) ) * delta_x / 6;
    }
    return sum;
}

// num. Integration von Cosh (Simpsonregel) 
double cosh_simpson(double x, double delta_x) {
    int N = fabs(x) / delta_x;
    double sum = 0, x_i, m_i, x_ii;
    for (int i = 0; i < N; i++) {
        x_i = delta_x * i;
        m_i = delta_x * (i+0.5);
        x_ii = delta_x * (i+1);
        sum += ( cosh(x_i) + 4*cosh(m_i) + cosh(x_ii) ) * delta_x / 6;
    }
    return sum;
}



int main(void) {
    double x = 2.0, delta_x = 3*pow(10,-4);
    double Msum = erf_midpoint(x,delta_x);
    double Ssum = erf_simpson(x,delta_x);
    printf("Integralwert über Mitellpunktregel: %f\n", Msum);
    printf("Integralwert über Simpsonsregel: %f\n", Ssum);

    FILE* file_erf_simpson = fopen("data(A03).csv", "w");
    fprintf(file_erf_simpson, "Integralgrenze x, Ergebnis S_erf, Ergebnis M_erf\n");
    for (int i = 0; i < 100; i++) {
        x = -2 + i / (99.0 / 4);
        double erf_integrated_simpson = erf_simpson(x, delta_x);
        double erf_integrated_midpoint = erf_midpoint(x, delta_x);
        fprintf(file_erf_simpson, "%.8f, %.8f, %.8f\n", x, erf_integrated_simpson, erf_integrated_midpoint);
    }
    fclose(file_erf_simpson);

    double I = sinh(1), cosh_integrated_midpoint, cosh_integrated_simpson;
    FILE* file_cosh_numint = fopen("data(A03_cosh).csv", "w");
    fprintf(file_cosh_numint, "Intervallbreite delta_x, Fehler_S, Fehler_M\n");
    for (int exp = 0; exp < 10; exp++) {
        delta_x = pow(2, -exp);
        cosh_integrated_midpoint = cosh_midpoint(1, delta_x);
        cosh_integrated_simpson = cosh_simpson(1, delta_x);
        fprintf(file_cosh_numint, "%.8f, %.8f, %.8f\n", delta_x, fabs(I - cosh_integrated_simpson), fabs(I - cosh_integrated_midpoint));
    }
    fclose(file_cosh_numint);

    return EXIT_SUCCESS;
}

