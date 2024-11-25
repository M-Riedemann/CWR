#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../my_numerics.h"
#include <float.h>


typedef enum{
    INSIDE,
    DIRICHLET,
    NEUMANN, // We won’t be using these
    DISABLED
    } Node_Type;

typedef struct{
    Node_Type type;
    double val;
    } node;

int InitDomain(int N, node F[][N], double N_val){
    F[0][0].val = 100;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j == 0 || i==0 || j == N-1 || i == N-1){
                F[i][j].type = DIRICHLET; 
                F[i][j].val = N_val;
            } 
            else {
                F[i][j].type = INSIDE; 
                F[i][j].val = 0.; 
            }              
        }
    }
    return 0;
}

int InsertCircle(int N, double cx, double cy, double R, node F[][N], double val){
    double delta = 1.0/(N-1);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double r = sqrt(pow((i*delta-cx),2)+pow((j*delta-cy),2));
            if(r<R){
                F[i][j].type = DIRICHLET; 
                F[i][j].val = val;
            }
        }
    }
    return 0;
}


void PoissonGaussSeidel(int N, node F[][N], int iter_max, double tolerance){
    double L=0;
    double R_phi;
    double R_quad = pow((tolerance + 1.0),2);
    double phi;
    while(L < iter_max && sqrt(R_quad) > tolerance){
        R_quad = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                switch (F[i][j].type){
                    case INSIDE:
                    // Hier einen inneren Gitterpunkt aktualisieren 
                    phi = F[i][j].val; 
                    F[i][j].val = (1/4.0)*(F[i][j-1].val+F[i][j+1].val+F[i-1][j].val+F[i+1][j].val);
                    R_phi = ((F[i][j].val-phi)/(phi+DBL_MIN));
                    R_quad += pow(R_phi,2);
                    break;
                    case DIRICHLET:
                    // Hier einen Randwert aktualisieren 
                    
                    break;
                default:
                break;
                }
            }
        }
        L+=1;
    }
}



int main(void){
    int N = 128;
    double cx_l = 0.25;
    double cy_l = 0.5;
    double cx_r = 0.75;
    double cy_r = 0.5;
    double R = 0.1;
    double val_l = -1;
    double val_r = 1;
    double N_val = 0;

    int iter_max = 1000;
    double tolerance = 1e-3;

    node (*F)[N] = malloc(sizeof(node[N][N]));
    
    InitDomain(N, F, N_val);

    InsertCircle(N, cx_l, cy_l, R, F, val_l);

    InsertCircle(N, cx_r, cy_r, R, F, val_r);

    FILE* GS_file = fopen("A21_GS.csv", "w");
    
    for (int j = 0; j < N; j++){
        for (int i = 0; i < N-1; i++){
            fprintf(GS_file, "%g, ", F[i][j].val);
        }
        fprintf(GS_file, "%g \n", F[N-1][j]);
    }
    
    PoissonGaussSeidel(N, F, iter_max, tolerance);

    FILE* GS2_file = fopen("A21_GS2.csv", "w");
    
    for (int j = 0; j < N; j++){
        for (int i = 0; i < N-1; i++){
            fprintf(GS2_file, "%g, ", F[i][j].val);
        }
        fprintf(GS2_file, "%g \n", F[N-1][j]);
    }
    

    return 0;
}
