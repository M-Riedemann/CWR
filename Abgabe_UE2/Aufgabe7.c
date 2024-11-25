#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"

//von der Aufgabe festgelegte Konstante
int w=4;

// Fourierintegrand Cosinus für x und Parameterarray mit omega und k
double fourier_integrand_cosine(double x, double k) {
    return cos(w*x) * cos(k*x);
}


int main(void) {
    //Definition der Benötigten Größen
    double dk = 20.0/1000.0, k, delta_x = pow(10,-3), M[] = {10,30,100,500};
    
    //Ausgabedatei
    FILE* file_f = fopen("data(A07).csv", "w");
    fprintf(file_f, "k, Integralwerte (für M = 10, 30, 100, 500)\n");

    //Iteration über 1000 Werte für zwischen [-10, 10]
    for(int i = 0; i < 1001; i++) {
        double k = -10+i*dk;
        fprintf(file_f, "%g", k);
        
        //Integration über die verschiedenen Integralgrenzen
        for(int i = 0; i < 4; i++) {
            double left = -M[i], right = M[i];
            double sum = integrate_simpson_2_params(left, right, delta_x, fourier_integrand_cosine, k);
            fprintf(file_f, ", %g", sum);}

        fprintf(file_f, "\n");}
    
    fclose(file_f);


    return 0;
}



