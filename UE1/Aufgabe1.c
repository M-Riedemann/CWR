#include <stdio.h>
#include <stdlib.h>

int main(void) {
    double left = 0.0;
    double right = 2.0;
    int N = 100;
    double sum = 0.0;
    double delta_y = (right - left)/ N;
    for(int k = 0; k<N; ++k) {
        double y = left + k * delta_y;
        double f = y*y*y - y/2;
        double A = f * delta_y;
        sum += A;
    }
    printf("Integralwert: %f\n", sum);

    return EXIT_SUCCESS;
}

