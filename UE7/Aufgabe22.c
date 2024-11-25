#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"


const double L = 1;                                                 // Länge der Domäne
const double delta_x = 0.01;                                        // Diskretisierung der Domäne
const double D = 0.1;                                               // Wärmediffusivität
const double T_LINKS = 1;                                           // Direchlet-Randbedingung links: T = 1
const double T_RECHTS = -1;                                         // Direchlet-Randbedingung rechts: T = -1
const int N = L / delta_x;                                          // Anzahl der Diskretisierungsschritte




int heat_FTCS(double t, double y[], double dt) {
    double *y_temp = (double*) malloc(N * sizeof(double));
    for(size_t i = 0; i < N; i++) {
        double T = y[i % N];
        double TW = y[(i - 1) % N];
        double TO = y[(i + 1) % N];

        if (i == 0) {
            y_temp[i] = y[i] + dt * D * (T_LINKS - 2*T + TO) / pow(delta_x, 2);
        } else if (i == N - 1) {
            y_temp[i] = y[i] + dt * D * (TW - 2*T + T_RECHTS) / pow(delta_x, 2);
        } else {
            y_temp[i] = y[i] + dt * D * (TW - 2*T + TO) / pow(delta_x, 2);
        }
    }

    // Übertragen der neuen Werte in das Ursprungsaray
    for(int i = 0; i < N; i++) {
        y[i] = y_temp[i];
    }

    free(y_temp);
    return 0;
}


int main(void) {
    // Zeitparameter
    double time = 0;
    double time_max = 1;
    double delta_time = pow(delta_x, 2) / (2 * D) / 4;         
    printf("delta_time: %g\n", delta_time);

    double *y = (double*) calloc(N, sizeof(double));                // Initialisieren des Temperaturfeldes

    // Initialisieren des Files
    FILE* heat_file = fopen("A22_Waermeleitungsgleichung_data.csv", "w");
    fprintf(heat_file, "0");
    for (int i = 0; i < N; i++) {
        fprintf(heat_file, ", %g", i*delta_x);
    }
    fprintf(heat_file, "\n%g", time);
    for (int i = 0; i < N; i++) {
        fprintf(heat_file, ", %g", y[i]);
    }

    // Durchlaufen der Zeitschritte bis time_max
    while (time < time_max) {
        heat_FTCS(time, y, delta_time);
        time += delta_time;

        // Beschreiben des Heat-Files
        fprintf(heat_file, "\n%g", time);
        for (int i = 0; i < N; i++) {
            fprintf(heat_file, ", %g", y[i]);
        }
    }

    free(y);
    fclose(heat_file);
    return 0;
}
