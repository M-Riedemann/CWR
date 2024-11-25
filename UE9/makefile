#Flags
FLAGS := -lm -Wall -g -lgsl -lgslcblas #-fsanitize=address -v

#Aufgabe30
Aufgabe30: Aufgabe30.o my_numerics.o mlr_rand.o
	gcc Aufgabe30.o my_numerics.o mlr_rand.o -o Aufgabe30 $(FLAGS)
	rm *.o

Aufgabe30.o: Aufgabe30.c
	gcc -c Aufgabe30.c $(Flags)

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c $(Flags)

mlr_rand.o: ../mlr_rand.h ../mlr_rand.c
	gcc -c ../mlr_rand.c $(Flags)


run30: Aufgabe30
	./$< Aufgabe30
