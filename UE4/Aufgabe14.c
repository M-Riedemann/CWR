#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"


/* Konstanten*/
const int N = 2;                                                                // Dimension
const double mass = 1;                                                          // Masse 
const double k = 10;                                                            // Federkonstante beider Federn           
const double r_0 = 1;                                                           // Ruhelänge beider Federn
const double r_A[2] = {-1, 0};                                                  // Position der Verankerung der linken Feder
const double r_B[2] = {1, 0};                                                   // Position der Verankerung der rechten Feder

// +++ Aufgabe 1: Bestimmung der DGL +++

// Dgl des Federsystems
int ODE_dual_springs(double t, const double y[], double f[], void *params) {

    // Übertragung der Geschwindigkeiten
    f[0] = y[2];
    f[1] = y[3];

    // Bestimmung der Abstände zu der jeweiligen Verankerung
    double dist_r_A = sqrt((y[0]-r_A[0])*(y[0]-r_A[0])+(y[1]-r_A[1])*(y[1]-r_A[1]));
    double dist_r_B = sqrt((y[0]-r_B[0])*(y[0]-r_B[0])+(y[1]-r_B[1])*(y[1]-r_B[1]));

    // Berechnung der Beschleunigungen
    double a_A = -(k * (dist_r_A - r_0)) / mass;
    double a_B = -(k * (dist_r_B - r_0)) / mass;

    // Trennung der X und Y komponenten der Beschleunigung
    double a_A_x = a_A * (y[0] - r_A[0]) / dist_r_A;
    double a_A_y = a_A * (y[1] - r_A[1]) / dist_r_A;

    double a_B_x = a_B * (y[0] - r_B[0]) / dist_r_B;
    double a_B_y = a_B * (y[1] - r_B[1]) / dist_r_B;

    // Übertragung der Beschleunigung
    f[2] = (a_A_x + a_B_x);
    f[3] = (a_A_y + a_B_y);
    return 0;
}


int main(void) {
    // Simulationsparameter
    int dimension = 2 * N;                                                      // Dimensionalität des Zustandsvektors des Systems
    double T_x = 2 * M_PI / sqrt(2*k);                                          // analytische Periodendauer T_x
    double T_x_Euler = 0;
    double T_x_RK4 = 0;
    double x_0 = -0.5;                                                          // Startkoordinate x_0
    double t = 0;                                                               // Startzeit t = 0
    double delta_t = 1e-4;                                                      // Größe der Zeitschritte delta_t

    // Initialisierung der Startposition der Masse
    double y_euler[4] = {x_0, 0}; 
    double y_rk4[4] = {x_0, 0};


    // +++ Aufgabe 4, 5: Position x in Abhängigkeit von der Zeit +++

    // numerische Nullstellensuche über Vorzeichenwechsel der x-Koordinate
    int euler_count_root = 0;                                        // Anzahl der Nulldurchgänge (Euler)
    int rk4_count_root = 0;                                          // Anzahl der Nulldurchgänge (RK4)
    double *euler_root = (double*) malloc(sizeof(double));            // Definition des Pointers für die Nulldurchgänge (Euler)
    double *rk4_root = (double*) malloc(sizeof(double));              // Definition des Pointers für die Nulldurchgänge (RK4)

    // zu überschreibende csv datei
    FILE* pos_file = fopen("Abgabe_UE4/A14_pos.csv", "w"); 
    fprintf(pos_file, "Zeit, x (Euler), y (Euler), x (RK4), y (RK4)\n");
    fprintf(pos_file, "%g, %g, %g, %g, %g\n", t, y_euler[0], y_euler[1], y_rk4[0], y_rk4[1]);

    // Integration über 10 Periodendauern T_x
    while (t <= 10 * T_x) {
       
        double x_pre_euler = y_euler[0];                                        // x-Koordinate vor euler-step 
        double x_pre_rk4 = y_rk4[0];                                            // x-Koordinate vor rk4-step

        euler_step(t, delta_t, y_euler, ODE_dual_springs, dimension, NULL);
        rk4_step(t, delta_t, y_rk4, ODE_dual_springs, dimension, NULL);

        double x_post_euler = y_euler[0];                                       // x-Koordinate vor euler-step
        double x_post_rk4 = y_rk4[0];                                           // x-Koordinate nach rk4-step 

        t += delta_t;

        // Bestimmung der Zeiten der Nulldurchgänge
        if (x_pre_euler * x_post_euler < 0) {                                   // Prüfung ob Nulldurchgang erfolgt anhand eines Vorzeichenwechsels bei der Euler-Integration  
            euler_count_root += 1;               
            euler_root = (double*) realloc(euler_root, sizeof(double) * euler_count_root);
            euler_root[euler_count_root - 1] = t;          // Eintragen des Nulldurchganges in vergrößerten pointer          

        }
        if (x_pre_rk4 * x_post_rk4 < 0) {                                       // Prüfung ob Nulldurchgang erfolgt anhand eines Vorzeichenwechsels bei der rk4-Integration
            rk4_count_root += 1;
            rk4_root = (double*) realloc(rk4_root, sizeof(double) * rk4_count_root);
            rk4_root[rk4_count_root - 1] = t;              // Eintragen des Nulldurchganges in vergrößerten pointer 

        }
        // Überschreibung der csv datei mit den Positionen
        fprintf(pos_file, "%g, %g, %g, %g, %g\n", t, y_euler[0], y_euler[1], y_rk4[0], y_rk4[1]);
    }

    // +++ Aufgabe 3 (Ausgabe) +++

    // Ausgabe der Periodendauer
    printf("analytisch bestimmte Periodendauer: %g\n", T_x);
    if (euler_count_root < 2) {                                      // Prüfen ob ausreichend Nulldurchgänge (Euler)
        printf("Fehler! Nicht genügent Nullstellen");
    } else {                                                                    // Berechnung der Periodendauer
        T_x_Euler = 2.0 / (euler_count_root - 1)  * (euler_root[euler_count_root - 1] - euler_root[0]);
        printf("Periodendauer T_x bei der Euler-Integration: %g\n", T_x_Euler);                   // Ausgabe der Periodendauer bei der euler-Integration
    }
    if (rk4_count_root < 2) {                                        // Prüfen ob ausreichend Nulldurchgänge (RK4)
        printf("Fehler! Nicht genügent Nullstellen");
    } else {                                                                    // Berechnung der Periodendauer                           
        T_x_RK4 = 2.0 / (rk4_count_root - 1) * (rk4_root[rk4_count_root - 1] - rk4_root[0]);
        printf("Periodendauer T_x bei der RK4-Integration: %g\n", T_x_RK4);                       // Ausgabe der Periodendauer bei der rk4-Integration
    }
    T_x_RK4 = 2.0 / (rk4_count_root -1) * (rk4_root[rk4_count_root - 1] - rk4_root[0]);
    fclose(pos_file), free(euler_root), free(rk4_root);


    // +++ Aufgabe 6: Abweichung der numerischen von der analytischen Lösung für variierte Zeitschritte delta_t +++

    // zu überschreibende csv datei
    FILE* res_file = fopen("Abgabe_UE4/A14_res.csv", "w");
    fprintf(res_file, "delta_t, Residue (Euler), Residue (RK4)\n"); 

    // Iteration über 150 logarithmisch gleichverteilte Zeitschritte 
    for (int i = 0; i < 150; i++) {

        t = 0, delta_t = pow(10, -1-4*i/149.0);                                  

        // initialisierung der Startposition
        double y_euler_res[4] = {x_0, 0}; 
        double y_rk4_res[4] = {x_0, 0};

        // Integration über 10 Periodendauern T_x
        while (t + delta_t < 10 * T_x_RK4) {
            t += delta_t;
            euler_step(t, delta_t, y_euler_res, ODE_dual_springs, dimension, NULL);
            rk4_step(t, delta_t, y_rk4_res, ODE_dual_springs, dimension, NULL);
        }
        
        // zusätzlicher Zeitschritt für Vergleichbarkeit mit analytischer Lösung
        euler_step(t, 10*T_x_RK4 - t, y_euler_res, ODE_dual_springs, dimension, NULL);
        rk4_step(t, 10*T_x_RK4 - t, y_rk4_res, ODE_dual_springs, dimension, NULL);

        // Überschreibung der csv datei mit den Residuen
        fprintf(res_file, "%g, %g, %g\n", delta_t, fabs(x_0 - y_euler_res[0]), fabs(x_0 - y_rk4_res[0]));
        
    }
    fclose(res_file);


    // +++ Aufgabe 7: Numerische Berechnung der Trajektorie mit drei unterschiedlichen Anfangsbedingungen und dem überlegenden Integrationsverfahren RK4 +++
    t = 0, delta_t = 10e-4;

    // Initialisierung der Startpositionen der Masse
    double y1_rk4[4] = {-0.2, 0.5, 0.2};
    double y2_rk4[4] = {-0.7, 0.3, 0};
    double y3_rk4[4] = {0, -0.7, 1.4};

    // zu überschreibende csv datei
    FILE* weird_pos_file = fopen("Abgabe_UE4/A14_weird_pos.csv", "w"); 
    fprintf(weird_pos_file, "Zeit, x1, y1, x2, y2, x3, y3\n");
    fprintf(weird_pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y1_rk4[0], y1_rk4[1], y2_rk4[0], y2_rk4[1], y3_rk4[0], y3_rk4[1]);

    // Integration über 10 Periodendauern T_x
    while (t <= 10 * T_x) {
        t += delta_t;
        rk4_step(t, delta_t, y1_rk4, ODE_dual_springs, dimension, NULL);
        rk4_step(t, delta_t, y2_rk4, ODE_dual_springs, dimension, NULL);
        rk4_step(t, delta_t, y3_rk4, ODE_dual_springs, dimension, NULL);
        fprintf(weird_pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y1_rk4[0], y1_rk4[1], y2_rk4[0], y2_rk4[1], y3_rk4[0], y3_rk4[1]);
    }
    fclose(weird_pos_file);
    return 0;
}

