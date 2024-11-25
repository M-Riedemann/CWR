#Flags
FLAGS := -lm -Wall -g #-v -lgsl -lgslcblas -fsanitize=address

#Aufgabe20
Aufgabe20: Aufgabe20.o my_numerics.o
	gcc Aufgabe20.o my_numerics.o -o Aufgabe20 $(FLAGS)
	rm *.o

Aufgabe20.o: Aufgabe20.c
	gcc -c Aufgabe20.c $(Flags)

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c $(Flags)

run20: Aufgabe20
	./$< Aufgabe20


