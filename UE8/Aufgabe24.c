#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"


// Physikalische Parameter des Systems
const double BOX_X = 1;
const double BOX_Y = 1;
const double V_RAND = 0;


// Elektroden
const double R_ELECTRODE = 0.1;
const double x_LEFT_ELECTRODE = 0.25;          
const double x_RIGHT_ELECTRODE = 0.75;
const double y_LEFT_ELECTRODE = 0.5;
const double y_RIGHT_ELECTRODE = 0.5;
const double V_LEFT_ELECTRODE = -1;
const double V_RIGHT_ELECTRODE = 1;


// Hilfsparameter
const double SMOOTHING_FACTOR = 1e-10;


// Typ eines Gitterpunktes
enum Node_Type {
    INSIDE,
    DIRICHLET,
    NEUMANN,
    DISABLED
};


// Struktur eines Gitterpunktes: Typ und Wert
struct node {
    enum Node_Type type;
    double value;
};


// Initialisierung eines quadratischen Arrays F[N][N]
int InitDomain(int N, struct node F[][N], double N_val, double E_val, double S_val, double W_val) {
    for (int n_x = 0; n_x < N; n_x++) {
        for (int n_y = 0; n_y < N; n_y++) {
            if (n_x == 0 || n_x == N-1 || n_y == 0 || n_y == N-1) {
                F[n_x][n_y].type = DIRICHLET;
                F[n_x][n_y].value = V_RAND;
            } else {
                F[n_x][n_y].type = INSIDE;
                F[n_x][n_y].value = 0;
            }
        }
    }
    return 0;
}


// Einsetzten eines Potentialkreises mit Radius < R um c_x, c_y
int InsertCircle(int N, double c_x, double c_y, double R, struct node F[][N], double val) {
    double delta_x = BOX_X / N;
    double delta_y = BOX_Y / N;
    for (int n_x = 0; n_x < N; n_x++) {
        for (int n_y = 0; n_y < N; n_y++) {
            if (sqrt(pow((n_x*delta_x-c_x),2)+pow((n_y*delta_y-c_y),2)) < R) {
                F[n_x][n_y].type = DIRICHLET;
                F[n_x][n_y].value = val;
            }
        }
    }
    return 0;
}






// Gauß-Seidel-Verfahren für die Poisson-Gleichung
int PoissonGaussSeidel(int N, struct node F[][N], int iter_max, double tolerance) {
    int count = 0;                                                                                      // Anzahl der Schleifendurchläufe
    double R_squared;                                                                                   // quadratische relative Anderung aller Gitterpunkte
    do {
        R_squared = 0;
        for (int n_x = 0; n_x < N; n_x++) {
            for (int n_y = 0; n_y < N; n_y++) {   
                switch (F[n_x][n_y].type) {
                    case INSIDE:
                        double F_old = F[n_x][n_y].value;
                        double F_new = 1.0 / 4 * (F[n_x][n_y-1].value + F[n_x][n_y+1].value + F[n_x-1][n_y].value + F[n_x+1][n_y].value);
                        F[n_x][n_y].value = F_new;
                        R_squared += pow((F_new-F_old) / (F_old + SMOOTHING_FACTOR), 2);           // neue relative Änderung für innere Werte
                        break;

                    case DIRICHLET:
                        break;

                    default:
                        break;
                }
            }
        }
        count++;
    } while (count < iter_max && sqrt(R_squared) > tolerance);
    return count;
}


// SOR-Verfahren für die Poisson-Gleichung
int PoissonSOR(int N, struct node F[][N], int iter_max, double tolerance, double alpha) {
    int count = 0;                                                                                      // Anzahl der Schleifendurchläufe
    double R_squared;                                                                                   // quadratische relative Anderung aller Gitterpunkte
    do {
        R_squared = 0; 
        for (int n_x = 0; n_x < N; n_x++) {
            for (int n_y = 0; n_y < N; n_y++) {   
                switch (F[n_x][n_y].type) {
                    case INSIDE:
                        double F_old = F[n_x][n_y].value;
                        double F_new = 1.0 / 4 * (F[n_x][n_y-1].value + F[n_x][n_y+1].value + F[n_x-1][n_y].value + F[n_x+1][n_y].value);
                        F[n_x][n_y].value = F_old + alpha * (F_new - F_old);
                        R_squared += pow((F_new-F_old) / (F_old + SMOOTHING_FACTOR), 2);                                // neue relative Änderung für innere Werte
                        break;

                    case DIRICHLET:
                        break;

                    default:
                        break;
                }
            }
        }
        count++;
    } while (count < iter_max && sqrt(R_squared) > tolerance);
    return count;
}


int main(void) {
    // Parameter des SOR-Verfahrens
    double tolerance = 1e-4;
    int iter_max = 1e6;

    // Größen der Felder
    int N16 = 16;
    int N32 = 32;
    int N64 = 64;
    int N128 = 128;

    // Datei der Konvergenzdaten
    FILE* SOR_convergence_file = fopen("data/A24_SOR_convergence_data.csv", "w");
    fprintf(SOR_convergence_file, "alpha, convergence N16, N32, N64, N128\n");

    // Erstellen der Felder
    struct node (*F16)[N16] = malloc(sizeof(struct node[N16][N16]));
    struct node (*F32)[N32] = malloc(sizeof(struct node[N32][N32]));
    struct node (*F64)[N64] = malloc(sizeof(struct node[N64][N64]));
    struct node (*F128)[N128] = malloc(sizeof(struct node[N128][N128]));


    // Feld für SOR und Anzahl Iterationen bis Konvergenz
    for (double alpha = 1; alpha < 2 - 0.01; alpha += 0.0125) {

        // Initialisierung der Ränder
        InitDomain(N16, F16, 0, 0, 0, 0);
        InitDomain(N32, F32, 0, 0, 0, 0);
        InitDomain(N64, F64, 0, 0, 0, 0);
        InitDomain(N128, F128, 0, 0, 0, 0);

        // Hinzufügen der Elektroden
        InsertCircle(N16, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F16, V_LEFT_ELECTRODE);                      // Einfügen linke Elektrode
        InsertCircle(N16, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F16, V_RIGHT_ELECTRODE);                   // Einfügen rechte Elektrode

        InsertCircle(N32, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F32, V_LEFT_ELECTRODE);                      // Einfügen linke Elektrode
        InsertCircle(N32, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F32, V_RIGHT_ELECTRODE);                   // Einfügen rechte Elektrode

        InsertCircle(N64, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F64, V_LEFT_ELECTRODE);                      // Einfügen linke Elektrode
        InsertCircle(N64, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F64, V_RIGHT_ELECTRODE);                   // Einfügen rechte Elektrode

        InsertCircle(N128, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F128, V_LEFT_ELECTRODE);                    // Einfügen linke Elektrode
        InsertCircle(N128, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F128, V_RIGHT_ELECTRODE);                 // Einfügen rechte Elektrode

        // Konvergenz des SOR-Verfahrens
        printf("calculating alpha %g\t", alpha);
        printf("...F16");
        int n_iterations_SOR_16 = PoissonSOR(N16, F16, iter_max, tolerance, alpha);
        printf("...F32");
        int n_iterations_SOR_32 = PoissonSOR(N32, F32, iter_max, tolerance, alpha);
        printf("...F64");
        int n_iterations_SOR_64 = PoissonSOR(N64, F64, iter_max, tolerance, alpha);
        printf("...F128");
        int n_iterations_SOR_128 = PoissonSOR(N128, F128, iter_max, tolerance, alpha);
        printf("...finished\n");

        fprintf(SOR_convergence_file, "%g, %d, %d, %d, %d\n", alpha, n_iterations_SOR_16, n_iterations_SOR_32, n_iterations_SOR_64, n_iterations_SOR_128);               
    }

    fclose(SOR_convergence_file);
    free(F16);
    free(F32);
    free(F64);
    free(F128);
    return 0;
}
