

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"


/*------------------------  PHYSIKALISCHE KONSTANTEN  -------------------------*/

const int N = 1 ;              // Anzahl der Pendel
const double k = 100. ;           // Federhaerte [N/m]
const double base_length = 1. ; // Basislaenge der Federn [m]
const double mass = 1. ;        // Masse der Pendel [kg]
const double v_0 = 1 ;
/*-------------------------  SIMULATIONS-PARAMETER  ---------------------------*/

const double T_max = 20.0;
const double delta_t = 1e-2; // Zeitliche Schrittweite

static const char pos_file_name[] = "UE5/position15.csv";
static const char energy_file_name[] = "UE5/energy15.csv";

/*-------------------------  PHYSIKALISCHES SYSTEM  ---------------------------*/



/*int pendulumsODE(double t, const double y[], double f[], void *params)
{   
    for(int i = 0; i < 2*N; ++i)
        f[i] = 0;

    for (int i = 0; i < N; i++) {                   
        f[i] = y[N+i];
    }
    
    f[N] = -k*y[0] / mass;                                  // Beschleunigungen aus Positionen berechnen
    
    
    if (N != 1) {                                   // Position 0
        f[N] += k*(y[1] - y[0]) / mass;                            
    }
    for (int i = 1; i < N - 1; i++) {                       // Positionen 1 bis N - 2
        f[N+i] = -k*(y[i] - y[i-1]) / mass + k*(y[i+1] - y[i]) / mass; 
    }
    
    if (N-1 != 1) {                                         
        f[2*N-1] = -k*(y[2*N-1] - y[2*N-2]) / mass;         // Position N 
    }
    return 0;
}*/

int pendulumsODE(double t, const double y[], double f[], void *params) {
    // Geschwindigkeiten aus dem Zustandsarray übertragen
    for (int i = 0; i < N; i++) {                   
        f[i] = y[N+i];
    }

    // Beschleunigungen aus Positionen berechnen
    f[N] = -k*y[0] / mass;                                              // Position 0 generell
    if (N != 1) {
        f[N] += k*(y[1] - y[0] - base_length);                          // Position 1 auf Position 0 wenn mehr als 1 Pendel
        f[2*N-1] += -k*(y[N-1] - y[N-2] - base_length) / mass;          // Position N - 1 (wenn mehr als 1 Pendel)
    }
    for (int i = 1; i < N - 1; i++) {                                   // Positionen 1 bis N - 2
        f[N+i] += -k*(y[i] - y[i-1] - base_length) / mass + k*(y[i+1] - y[i] - base_length) / mass; 
    }
    return 0;
}


double pendulums_energy(const double y[]) {
    double E_pot = 0;
    double E_kin = 0;

    for (int i = 0; i < N; i++) {                           // Kinetische Energie berechnen
        E_kin += mass / 2 * y[N+i]*y[N+i];
    }

    E_pot += k / 2 * y[0]*y[0];                 // Potentielle Energie
    if (N != 1) {
        E_pot += k / 2 * (y[1] - y[0] - 1)*(y[1] - y[0] - 1);
        E_pot += k / 2 * (y[N-1] - y[N-2] - 1)*(y[N-1] - y[N-2] - 1);                                     
    }
    for (int i = 1; i < N - 1; i++) {                          
        E_pot += k / 2 * (y[i] - y[i-1] - 1)+(y[i] - y[i-1] - 1) + k / 2 * (y[i+1] - y[i] - 1)*(y[i+1] - y[i] - 1); 
    }
    return E_kin + E_pot;
}


int main(void)
{
    // analytische Parameter
    double T = 2 * M_PI * sqrt(mass / k);
    double x_T4 = v_0 * sqrt(mass / k);
    double T4 = T/4;

    // Dimensionalitaet des DGL-Systems 
    int dimension = 2*N ;

    // Initialisierung des Arrays
    double y[dimension];
    
    // fuellen der Arrays
    for (int i = 0; i < N; i++) {
        y[i] = i;
    }
    // Anfangsbedingung
    y[N] = v_0;
    for (int i = N+1; i<2*N; i++){
        y[i] = 0;
    }
    double t = 0;

    // Ausgabe Dateien 
    FILE *pos_file = fopen(pos_file_name, "w");
    FILE *energy_file = fopen(energy_file_name, "w");
    FILE* res_file = fopen("UE5/residuen15.csv", "w");
    fprintf(res_file, "Zeitschritt delta_t, Residuum RK2, , Residuum RK4, Residuum Euler, Residuum Verlet\n");

    // Parameter zur bestimmung der Nullstellen
    int count_root;
    double *rk4_root = (double*) malloc(sizeof(double));
    double even_root;
    double uneven_root; 

    // Parameter zu Residuenbestimmung 
    double d_t, exp = 0, residue_RK2, residue_RK4, residue_Euler, residue_verlet;

    // Pendel mit rk2
    while (t < T_max)
    {   
        double x_pre = y[N];

        // rk2 schritt
        rk2_step(t, delta_t, y, pendulumsODE, dimension, NULL);

        double x_post = y[N];

        if(x_pre*x_post<0){
            if(count_root % 2 == 0){
                even_root = t;
            }
            else{
                uneven_root = t;
            }
            count_root += 1;
        }

        //Energie 
        double energy = pendulums_energy( (const double*) y);

        // Ausgabe
        fprintf(pos_file, "%g", t);
        fprintf(energy_file, "%g, %g\n", t, energy);

        for (int i = 0; i < N; i++) {
            fprintf(pos_file, ", %g", y[i]);
        }
        fprintf(pos_file, "\n");

        t += delta_t;
    }
    
    double T_rk2;

    if (count_root < 2) {                                      // Prüfen ob ausreichend Nulldurchgänge (Euler)
        printf("Fehler! Nicht genügent Nullstellen");
    } 
    else {                                                                    // Berechnung der Periodendauer
            T_rk2 = 2.0*even_root/(count_root-1);
            printf("Periodendauer: %g\n", T_rk2);
                         
    }

    fclose(pos_file);
    fclose(energy_file);
    for (int i = 0; i < 100; i++) {
        t = 0;                                                          // Starten der jeweiligen Integration bei t = 0
        exp = -1 - i * (7.0/99);                                        // Anpassen des Exponentens des Schrittweite
        d_t = pow(10, exp);                                             // neue Schrittweite aus Exponent

        // Reinitialisierung Zustandsarray
        double y_rk2[2] = {0, v_0};
        double y_rk4[2] = {0, v_0};
        double y_euler[2] = {0, v_0};
        double y_verlet[2] = {0, v_0};

        // Numerische Integration bis vorletztem Integrationsschritt
        while (t + d_t < T4) {
            rk2_step(t, d_t, y_rk2, pendulumsODE, dimension, NULL);
            //printf("x (RK29: %g\n", y_rk2[0]);
            rk4_step(t, d_t, y_rk4, pendulumsODE, dimension, NULL);
            euler_step(t, d_t, y_euler, pendulumsODE, dimension, NULL);
            verlet_step(t, d_t, y_verlet, pendulumsODE, dimension, NULL);
            t += d_t;
        }

        // letzter Integrationsschritt mit verbleibender Zeit
        rk2_step(t, T4-t, y_rk2, pendulumsODE, dimension, NULL);
        rk4_step(t, T4-t, y_rk4, pendulumsODE, dimension, NULL);
        euler_step(t, T4-t, y_euler, pendulumsODE, dimension, NULL);
        verlet_step(t, T4-t, y_verlet, pendulumsODE, dimension, NULL);

        // Berechnung und Ausgabe des "Residuums"         
        residue_RK2 = fabs(y_rk2[0] - x_T4);
        residue_RK4 = fabs(y_rk4[0] - x_T4);
        residue_Euler = fabs(y_euler[0] - x_T4);
        residue_verlet = fabs(y_verlet[0] - x_T4);
        fprintf(res_file, "%g, %g, %g, %g, %g\n", d_t, residue_RK2, residue_RK4, residue_Euler, residue_verlet);
    }
    
    fclose(res_file);

    return EXIT_SUCCESS;
}

