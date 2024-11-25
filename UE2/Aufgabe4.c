#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"

double erf(double y) {
    return 2 * exp(-pow(y, 2)) / sqrt(M_PI);
}

int main(void) {
    double x = 2.0, delta_x = 3*pow(10,-4);
    double Msum = midpoint_integrator(x,delta_x, erf);
    double Ssum = simpson_integrator(x,delta_x, erf);
    printf("Integralwert über Mitellpunktregel: %f\n", Msum);
    printf("Integralwert über Simpsonsregel: %f\n", Ssum);

    FILE* file_erf_simpson = fopen("data(A03).csv", "w");
    fprintf(file_erf_simpson, "Integralgrenze x, Ergebnis S_erf, Ergebnis M_erf\n");
    for (int i = 0; i < 100; i++) {
        x = -2 + i / (99.0 / 4);
        double erf_integrated_simpson = simpson_integrator(x, delta_x, erf);
        double erf_integrated_midpoint = midpoint_integrator(x, delta_x, erf);
        fprintf(file_erf_simpson, "%.8f, %.8f, %.8f\n", x, erf_integrated_simpson, erf_integrated_midpoint);
    }
    fclose(file_erf_simpson);


    return EXIT_SUCCESS;
}
