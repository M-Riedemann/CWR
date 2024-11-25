#Flags
FLAGS := -lm -Wall #-v -lgsl -lgslcblas


#Aufgabe4
Aufgabe4: Aufgabe4.o my_numerics.o
	gcc Aufgabe4.o my_numerics.o -o Aufgabe4 $(FLAGS)
	rm *.o

Aufgabe4.o: Aufgabe4.c
	gcc -c Aufgabe4.c $(Flags)

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c $(Flags)

run4: Aufgabe4
	./$< Aufgabe4


#Aufgabe6
Aufgabe6: Aufgabe6.o my_numerics.o
	gcc Aufgabe6.o my_numerics.o -o Aufgabe6 $(FLAGS)
	rm *.o

Aufgabe6.o: Aufgabe6.c
	gcc -c Aufgabe6.c

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c

run6: Aufgabe6
	./$< Aufgabe6

