#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"


// Quantenzahlen und physikalische Konstanten
const int m = 1;
const int k = 1;
const int h_cross = 1;

//ODE des oszillators
int schroedinger_harmonic_ODE(double x, const double y[], double f[], void *params) {
    double E = *(double*) params;                                                           // Übergabe des E Pointers
    double V = (1.0/2.0) * k * pow(x,2);                                                          // Definition des Potentials
    f[0] = y[1];                                                                            // Übergabe der ersten Ableitung 
    f[1] = -2*m/pow(h_cross,2)*(E-V) * y[0];                                                // Übergabe der zweiten Ableitung
    return 0;
}

double test_energy(double E){
    double Earr[1] = {E}; 
    double psi[2] = {1,0};                                                      // Anfangswert für psi
    double x = -10.;                                                            // Startwert für X
    double delta_x = 1e-3;                                                      // Schrittweite

    while (x < 0){                                                              // Iteration bis x zu >= 0
        rk4_step(x, delta_x, psi, schroedinger_harmonic_ODE, 2, Earr);          // rk4 step
        x += delta_x;
    }

    rk4_step(x, 0 - x, psi, schroedinger_harmonic_ODE, 2, Earr);                 // extra step

    return psi[0]*psi[1];                                                        //test
}


int main(void){

    double E = 0.1;                                                                     // Anfangswert für E
    double E_s = 0.1;                                                   
    double delta_E = 1e-2;                                                              // iterationsschritt
    double rel_tol = 1e-3;                                                              // Relative Toleranz 
    double max_iter = 1e4;                                                              // maximale Iterationsschritte
    int n=0;
    double root_array[2];
   
    FILE* test_energy_file = fopen("A20_test_energy.csv", "w");                         
    fprintf(test_energy_file, "Energiewert, Psi(x1) * Psi_Prime(x1)\n");

    while(E<=10){                                                                        // Variation von E 
        fprintf(test_energy_file, "%g, %g\n", E, test_energy(E));                        // Abspeicherung des Produktes von Ort und seinem differential sowie von E
        E += delta_E;
    }
    fclose(test_energy_file);                                                               

    while(E_s < 2){                                                                                 //Variation von E_s um Startpunkte der Nullstellensuche festzulegen
        double E_root = find_root_newton_raphson(test_energy, E_s, delta_E, rel_tol, max_iter);
        printf("Eigenenergie bei n=%i: %g\n", n, E_root);
        root_array[n] = E_root;
        n += 1;
        E_s += 1;
    }


    double x = -5;                                                                      // Festlegung der Startwerte zur Bestimmung der Wellengleichungen
    double psi_0[2] = {1,0};
    double psi_1[2] = {1,0};
    double delta_x = 1e-3;  


    FILE* psi_0_file = fopen("A20_psi_0.csv", "w");
    fprintf(psi_0_file, "x, Psi_0(x)\n");


    FILE* psi_1_file = fopen("A20_psi_1.csv", "w");
    fprintf(psi_1_file, "x, Psi_1(x)\n");


    while (x < 5){                                                                          // Iteration zu bis x >= 5
        rk4_step(x, delta_x, psi_0, schroedinger_harmonic_ODE, 2, &root_array[0]);          // Bestimmung der Wellengleichungen über den rk4-stepper
        rk4_step(x, delta_x, psi_1, schroedinger_harmonic_ODE, 2, &root_array[1]);                          
        fprintf(psi_0_file, "%g, %g\n",x, psi_0[0]);
        fprintf(psi_1_file, "%g, %g\n",x, psi_1[0]);
        x += delta_x;
    }

    return EXIT_SUCCESS;
}

