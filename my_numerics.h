#ifndef MYNUMERICS_H
#define MYNUMERICS_H

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


/*Gauss*/
double mygauss(double y);

/* Numerical Integration of scalar 1D function */
double integrate(double left, double right,
int N, double integrand(double));

/*Typ gewöhnliche DGL*/
typedef int cvc_ode_func(double, const double[], double[], void*);
typedef int cvc_sde_func(double, const double[], double[], double, void*);

/* Integratoren */
double midpoint_integrator(double x, double delta_x, double func(double));
double simpson_integrator(double x, double delta_x, double func(double));
double integrate_simpson_2_params(double left, double right, double delta_x, double func(double,  double), double k);
double integrate_midpoint_2D(int A(double, double), double a_x, double b_x, double a_y, double b_y, double delta_x, double f(double, double));



/* Implementations of the error function */
double erf_midpoint(double x, double delta_x);
double erf_simpson(double x, double delta_x);

/* Differenzieren */
double diff(double x, double delta, double func(double));
double centered_diff(double x, double delta, double func(double));

/*Nullstellensuche*/
double find_root_newton_raphson(double func(double), double x0, double delta, double rel_tol, int max_iter);

/* Euler Verfahren */
void euler_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);

/*Runge-Kutta 2. Ordnung*/
void rk2_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);

/*Runge Kutta 4. Ordnung*/
void rk4_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);

/*verlet*/
void verlet_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);




#endif

