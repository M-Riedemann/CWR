#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(void) {

    // Parameter
    const double alpha = 1 / 6.1;
    const double gamma2 = 1 / 3.7;

    const double N = 83e6;
    const double E0 = 30000;
    const double I0 = 9000;


    // R0 = 1.25
    double S = N - (E0 + I0), E = E0, I = I0, R = 0, t = 0, I_max = I0;
    double R0 = 1.25;
    double beta = R0 * gamma2;

    // Daten
    FILE* file_epidemic_modelling_1_25 = fopen("A09_1.25.csv", "w");
    fprintf(file_epidemic_modelling_1_25, "Zeit t in Tagen, Susceptible S, Exposed E, Infectious I, Removed R\n");
    fprintf(file_epidemic_modelling_1_25, "%f, %f, %f, %f, %f\n", t, S, E, I, R);

    // Integrations-Parameter
    double delta_t = 0.01;
    while (t < 365) {
        t = t + delta_t;
        S = S - beta * (I*S) / N * delta_t;
        E = E + (beta * (I*S) / N  - alpha*E) * delta_t;
        I = I + (alpha*E - gamma2*I) * delta_t;
        R = N - (S + E + R);
        fprintf(file_epidemic_modelling_1_25, "%f, %f, %f, %f, %f\n", t, S, E, I, R);
        if (I > I_max) {
            I_max = I;
        }
    }
    
    printf("Maximale Anzahl Infizierter (R0 = 1.25): %f\n", I_max);
    fclose(file_epidemic_modelling_1_25);

    // R0 = 1.5
    S = N - (E0 + I0), E = E0, I = I0, R = 0, t = 0, I_max = I0;
    R0 = 1.5;
    beta = R0 * gamma2;

    // Daten
    FILE* file_epidemic_modelling_1_5 = fopen("A09_1.5.csv", "w");
    fprintf(file_epidemic_modelling_1_5, "Zeit t in Tagen, Susceptible S, Exposed E, Infectious I, Removed R\n");
    fprintf(file_epidemic_modelling_1_5, "%f, %f, %f, %f, %f\n", t, S, E, I, R);

    // Integrations-Parameter
    while (t < 365) {
        t += delta_t;
        S -= beta * (I*S) / N * delta_t;
        E += (beta * (I*S) / N  - alpha*E) * delta_t;
        I += (alpha*E - gamma2*I) * delta_t;
        R = N - (S + E + R);
        fprintf(file_epidemic_modelling_1_5, "%f, %f, %f, %f, %f\n", t, S, E, I, R);
        if (I > I_max) {
            I_max = I;
        }
    }
    printf("Maximale Anzahl Infizierter (R0 = 1.5): %f\n", I_max);
    fclose(file_epidemic_modelling_1_5);


    // R0 = 2
    S = N - (E0 + I0), E = E0, I = I0, R = 0, t = 0, I_max = I0;
    R0 = 2;
    beta = R0 * gamma2;

    // Daten
    FILE* file_epidemic_modelling_2 = fopen("A09_2.0.csv", "w");
    fprintf(file_epidemic_modelling_2, "Zeit t in Tagen, Susceptible S, Exposed E, Infectious I, Removed R\n");
    fprintf(file_epidemic_modelling_2, "%f, %f, %f, %f, %f\n", t, S, E, I, R);

    // Integrations-Parameter
    while (t < 365) {
        t += delta_t;
        S -= beta * (I*S) / N * delta_t;
        E += (beta * (I*S) / N  - alpha*E) * delta_t;
        I += (alpha*E - gamma2*I) * delta_t;
        R = N - (S + E + R);
        fprintf(file_epidemic_modelling_2, "%f, %f, %f, %f, %f\n", t, S, E, I, R);
        if (I > I_max) {
            I_max = I;
        }
    }
    printf("Maximale Anzahl Infizierter (R0 = 2): %f\n", I_max);
    fclose(file_epidemic_modelling_2);

    return 0;
}

