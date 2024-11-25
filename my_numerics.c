#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>
/*#include <gsl/gsl_rng.h>*/


// Strukturen
struct mlr_tuple_2 {
    double x1;
    double x2;
};

struct mlr_tuple_3 {
    double x1;
    double x2;
    double x3;
};




//erf
static double myerf(double y) {
    return 2 * exp(-pow(y, 2)) / sqrt(M_PI);
}



// Typ gewöhnliche DGL
typedef int ode_func(double, const double[], double[], void*);
typedef int sde_func(double, const double[], double[], double, void*);




//gauss 
double mygauss(double y) {
    return exp(-y*y);
}


// num. Integration (Mittelpunktsregel)
double erf_midpoint(double x, double delta_x) {
    int N = fabs(x)/delta_x;
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += myerf( (i-0.5)*delta_x) * delta_x;
    }
    return sum;
}


// num. Integration (Mittelpunktsregel)
double midpoint_integrator(double x, double delta_x, double func(double)) {
    int N = fabs(x)/delta_x;
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += func( (i-0.5)*delta_x) * delta_x;
    }
    return sum;
}


// 2-Dimensionale Integration (Mittelpunktsregel)
double integrate_midpoint_2D(int A(double, double), double a_x, double b_x, double a_y, double b_y, double delta_x, double f(double, double)) {
    int N_x = (b_x - a_x) / delta_x;
    int N_y = (b_y - a_y) / delta_x;
    double x, y, sum = 0;
    for (int i_x = 0; i_x < N_x; i_x++) {
        for (int i_y = 0; i_y < N_y; i_y++) {
            x = a_x + (i_x + 0.5) * delta_x;
            y = a_y + (i_y + 0.5) * delta_x;
            sum += f(x, y) * A(x, y);
        }
    }
    return pow(delta_x, 2) * sum;
}



// num. Integration (Simpsonregel)
double simpson_integrator(double x, double delta_x, double func(double)) {
    int N = fabs(x) / delta_x;
    double sum = 0, x_i, m_i, x_ii;
    for (int i = 0; i < N; i++) {
        x_i = delta_x * i;
        m_i = delta_x * (i+0.5);
        x_ii = delta_x * (i+1);
        sum += ( func(x_i) + 4*func(m_i) + func(x_ii) ) * delta_x / 6;
    }
    return sum;
}


// num. Integration (Simpsonregel mit Funktion func für Integralgrenzen (left, right), Schrittweite dx und Parameter *params
double integrate_simpson_2_params(double left, double right, double delta_x, double func(double,  double), double k) {
    double sum = 0, m_i, x_ii, ind = 1.0;
    if(right<left){
        double temp = right;
        right = left;
        left = temp;
        ind = -1.0;
    }
    double x_i = left;
    while (x_i < right) {
        m_i = x_i + 0.5*delta_x;
        x_ii = x_i + delta_x;
        sum += (func(x_i, k) + 4*func(m_i, k) + func(x_ii, k)) * delta_x / 6;
        x_i += delta_x;
    }
    return ind*sum;
}


// num. Integration (Simpsonregel)
double erf_simpson(double x, double delta_x) {
    int N = fabs(x) / delta_x;
    double sum = 0, x_i, m_i, x_ii;
    for (int i = 0; i < N; i++) {
        x_i = delta_x * i;
        m_i = delta_x * (i+0.5);
        x_ii = delta_x * (i+1);
        sum += ( myerf(x_i) + 4*myerf(m_i) + myerf(x_ii) ) * delta_x / 6;
    }
    return sum;
}


//einfache Ableitung
double diff(double x, double delta, double func(double)){
    return (func(x+delta)-func(x))/delta;
}

//zentrale differenz
double centered_diff(double x, double delta, double func(double)) {
    return ( func(x + delta) - func(x - delta) ) / ( 2 * delta);
}

//Newton_Raphson Methode der Nullstellensuche
double find_root_newton_raphson(double func(double), double x0, double delta, double rel_tol, int max_iter) {
    int i = 0;
    double x_old;
    while (i++ < max_iter) {
        x_old = x0;
        x0 = x0 - func(x0) / diff(x0, delta, func);
        if (fabs(x_old - x0) / fabs(x0) < rel_tol) {
            break;
        }
    }
    return x0;
}

// numerische Euler Integration des Zustandsarrays y mit gegebenen Parametern
void euler_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params) {
    double *f;
    f = (double*) calloc(dimension, sizeof(double));           // Reservierung des Ableitungsarrays
    func(t, y, f, params);                                      // Füllen des Ableitungsarrays über Aufruf der entprechenden ODE
    for (int i = 0; i < dimension; i++) {
        y[i] += f[i] * delta_t;
    }
    free(f);                                                    // Freigeben des Ableitungsarrays
    return;
}



// numerische Integration mittels Runge-Kutta 2. Ordnung
void rk2_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params) {
    double *support, *k1, *k2;
    support = (double*) calloc(dimension, sizeof(double));
    k1 = (double*) calloc(dimension, sizeof(double));
    k2 = (double*) calloc(dimension, sizeof(double));
    func(t, y, k1, params);                                     // Berechnung k1 = f(t, y)
    for (int i = 0; i < dimension; i++) {
        k1[i] *= delta_t;                                       // Berücksichtigung des Zeitschritts: k1 = f(t, y) * dt
        support[i] = y[i] + k1[i] / 2;                          // support = y + k1/2 (für nächsten Schritt)
    }
    func(t+delta_t/2, support, k2, params);                     // Berechnung k2 = f(t+dt/2, y+k1/2) und y_(i+1)
    for (int i = 0; i < dimension; i++) {
        k2[i] *= delta_t;
        y[i] += k2[i];
    }
    free(support), free(k1), free(k2);
    return;
}


// numerische Integration mittels Runge-Kutta 4. Ordnung 
void rk4_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params) {
    double *help, *k1, *k2, *k3, *k4;
    help = (double*) malloc(sizeof(double) * dimension);
    k1 = (double*) calloc(dimension, sizeof(double));
    k2 = (double*) calloc(dimension, sizeof(double));
    k3 = (double*) calloc(dimension, sizeof(double));
    k4 = (double*) calloc(dimension, sizeof(double));
    func(t, y, k1, params);                                       // Berechnung k1 = f(t, y) * dt und support = y + k1/2
    for (int i = 0; i < dimension; i++) {
        k1[i] *= delta_t;
        help[i] = y[i] + k1[i] / 2;
    }
    func(t+delta_t/2, help, k2, params);                     // Berechnung k2 = f(t+dt/2, y+k1/2) * dt und support = y + k2/2
    for (int i = 0; i < dimension; i++) {
        k2[i] *= delta_t;
        help[i] = y[i] + k2[i] / 2;
    }
    func(t+delta_t/2, help, k3, params);                     // Berechnung k3 = f(t+dt/2, y+k2/2) * dt und support = y + k3
    for (int i = 0; i < dimension; i++) {
        k3[i] *= delta_t;
        help[i] = y[i] + k3[i];
    }
    func(t+delta_t, help, k4, params);                       // Berechnung k4 = f(t+dt, y+k2) * dt und y_(i+1)
    for (int i = 0; i < dimension; i++) {
        k4[i] *= delta_t;
        y[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6; 
    }
    free(help), free(k1), free(k2), free(k3), free(k4);
    return;   
}


// numerische Integration
void verlet_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params) {
    int N = dimension / 2;
    double *a1, *a2;
    a1 = (double*) calloc(dimension, sizeof(double));
    a2 = (double*) calloc(dimension, sizeof(double));
    func(t, y, a1, params);                                     // Berechnung von a1 = f(t, y) * dt
    for (int i = 0; i < N; i++) {                               // Berechnung (erster Hälfte, Positionen) von y_(i+1) aus a1
        y[i] += a1[i] * delta_t + a1[i+N] * (delta_t * delta_t) /2;
    }
    func(t+delta_t, y, a2, params);                               // Berechnung von a2 = f(t+delta_t, y_(i+1)) aus Positionen von y_(i+1)
    for (int i = 0; i < N; i++) {                               // Berechnung (zweite Hälfte, Geschwindigkeiten) von y_(i+1) aus a1 und a2
        y[i+N] += (a1[i+N] + a2[i+N]) * delta_t / 2;
    }                                                                                                                 
    free(a1), free(a2);
    return;
}

