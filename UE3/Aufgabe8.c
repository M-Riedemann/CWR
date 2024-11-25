#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Analytische Lösungen der gegebenen Gleichung
double erg1_analytical = pow(10,8);
double erg2_analytical = pow(10,-8);

// Struktur eines 2er Tupels
struct tu {
    double erg1;
    double erg2;
};

// Lösung der Mitternachts-Formel
struct tu Mi_Fo(double a, double b, double c) {
    double erg1, erg2;

    erg1 = (-b + sqrt(pow(b, 2) - 4*a*c)) / (2*a);
    erg2 = (-b - sqrt(pow(b, 2) - 4*a*c)) / (2*a);
    struct tu MiFo;
    MiFo.erg1 = erg1;
    MiFo.erg2 = erg2;
    return MiFo;
}

// Lösung der erweiterten Mitternachts-Formel
struct tu Mi_Fo_extended(double a, double b, double c) {
    double erg1, erg2;
    erg1 = (2*c) / (-b - sqrt(pow(b, 2) - 4*a*c));
    erg2 = (2*c) / (-b + sqrt(pow(b, 2) - 4*a*c));
    struct tu MiFo_extended;
    MiFo_extended.erg1 = erg1;
    MiFo_extended.erg2 = erg2;
    return MiFo_extended;
}


// kombinierte Lösungsmethode quadratischer Gleichungen
struct tu solve_quadratic(double a, double b, double c) {
    double erg1, erg2;

    if (b > 0) {
        erg1 = (2*c) / (-b - sqrt(pow(b, 2) - 4*a*c));
        erg2 = (-b - sqrt(pow(b, 2) - 4*a*c)) / (2*a);
    } else {
        erg1 = (-b + sqrt(pow(b, 2) - 4*a*c)) / (2*a);
        erg2 = (2*c) / (-b + sqrt(pow(b, 2) - 4*a*c));
    }

    struct tu quadratic_solution;
    quadratic_solution.erg1 = erg1;
    quadratic_solution.erg2 = erg2;
    return quadratic_solution;
}



int main(void) {
    double a, b, c;
    a = 1;
    b = - (1e16 + 1) / 1e8;
    c =  1;
    struct tu MiFo, MiFo_extended, quadratic_solution;
    MiFo = Mi_Fo(a,b,c);
    MiFo_extended = Mi_Fo_extended(a,b,c);
    quadratic_solution = solve_quadratic(a,b,c);
    double re_error_1 = fabs(erg1_analytical - MiFo.erg1);
    double re_error_2 = fabs(erg2_analytical - MiFo.erg2);
    double error_1 = fabs(erg1_analytical - MiFo_extended.erg1);
    double error_2 = fabs(erg2_analytical - MiFo_extended.erg2);
    double p_error_1 = fabs(erg1_analytical - quadratic_solution.erg1);
    double p_error_2 = fabs(erg2_analytical - quadratic_solution.erg2);
    printf("++++ Ergebnis_1: %g Fehler_1: %g ++++ \n++++ Ergebnis_2: %g Fehler_2: %g ++++ \n", MiFo.erg1, re_error_1, MiFo.erg2, re_error_2);
    printf("Extended: \n++++ Ergebnis_1: %g Fehler_1: %g ++++ \n++++ Ergebnis_2: %g Fehler_2: %g ++++", MiFo_extended.erg1, error_1, MiFo_extended.erg2, error_2);
    printf("Perfected \n++++ Ergebnis_1: %g Fehler_1: %g ++++ \n++++ Ergebnis_2: %g Fehler_2: %g ++++ \n", quadratic_solution.erg1, p_error_1, quadratic_solution.erg2, p_error_2);
    return 0;
}

