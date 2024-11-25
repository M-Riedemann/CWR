#Flags
FLAGS := -lm -Wall -g #-v -lgsl -lgslcblas -fsanitize=address

#Aufgabe21
Aufgabe21: Aufgabe21.o my_numerics.o
	gcc Aufgabe21.o my_numerics.o -o Aufgabe21 $(FLAGS)
	rm *.o

Aufgabe21.o: Aufgabe21.c
	gcc -c Aufgabe21.c $(Flags)

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c $(Flags)

run21: Aufgabe21
	./$< Aufgabe21


#Aufgabe22
Aufgabe22: Aufgabe22.o my_numerics.o
	gcc Aufgabe22.o my_numerics.o -o Aufgabe22 $(FLAGS)
	rm *.o

Aufgabe22.o: Aufgabe22.c
	gcc -c Aufgabe22.c $(Flags)

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c $(Flags)

run22: Aufgabe22
	./$< Aufgabe22


