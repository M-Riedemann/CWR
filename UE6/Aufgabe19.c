#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"


// Konstanten des Systems
const double GAMMA = 4.0 / 3;
const double KAPPA = 3.85e9;
const double R_SUN = 7e8;
const double G = 6.67384e-11;
const double U_G_CONSTANT = 8.314;
const double VOLUME_MOL = 22.4e-3;
const double M_WASSERSTOFF = 1;
const double R = R_SUN * 1.2;

// Integrationsparameter
const double epsilon = 1e-4;
const int dimension = 2;


// Berechnung der Masse- und Dichteänderung für gegebenen Radius, sowie Masse-, und Dichtezustand
int ODE_star_density(double r, const double y[], double f[], void *params) {
    double rho = y[0];                                                                              // Dichte
    double mass = y[1];                                                                             // Masse
    f[0] = -3.0 / 2.0 / (KAPPA * GAMMA) * pow(rho, 2-GAMMA) * G * mass / (r*r);              // Dichteänderung
    f[1] = 4 * M_PI * rho * r*r;                                                       // Masseänderung
    return 0;
}


double R_tilde_funktion(double rho_0) {

    double m_0 = 4.0 / 3 * M_PI * rho_0 * pow(epsilon, 3);

    double y[2] = {rho_0, m_0}; 

    double r = epsilon;
    double delta_r = 1000;


    while (y[0] >= 1e-1) {
        
        rk4_step(r, delta_r, y, ODE_star_density, dimension, NULL);                                         

        r += delta_r;

        }
    
    return fabs(r-R);
}


int main(void){
    double rho_0 = 100000;
    double delta_rho = 1e3;

    double num_rho_0 = find_root_newton_raphson(R_tilde_funktion, rho_0, delta_rho, 1e5, 100000);
    
    printf("Numerisches rho_0: %g\n", num_rho_0);

    return EXIT_SUCCESS;
}

