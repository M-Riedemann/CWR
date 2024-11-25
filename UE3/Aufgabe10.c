#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void double_100(int *p) {
    for (int i = 0; i < 101; i++) {
        p[i] *= 2;
        }
}


int main(void) {
    
    int a[100];
    for(int i = 0; i < 101; ++i){
        a[i] = 10+i;
        printf("a[%d] = %d\n", i, a[i]);}

    int *p = malloc(100*sizeof(int));
    for(int i = 0; i < 101; ++i){
        p[i] = 10+i;
        printf("p[%d] = %d\n", i, p[i]);}

    // double_100(p);
    // for(int i = 0; i < 101; ++i){
        // printf("p[%d] = %d\n", i, p[i]);}
    
    printf("pointer = %lu\n array = %lu \n", sizeof(p), sizeof(a));

    p += 1;

    printf("p[%d] = %d\n", 0, p[0]);

     // MEME 
    printf("a[10] = %d\n", a[10]);
    printf("10[a] = %d\n", 10[a]);

    return 0;
}

