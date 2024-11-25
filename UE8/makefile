#Flags
FLAGS := -lm -Wall -g -lgsl -lgslcblas #-fsanitize=address -v

#Aufgabe24
Aufgabe24: Aufgabe24.o my_numerics.o
	gcc Aufgabe24.o my_numerics.o -o Aufgabe24 $(FLAGS)
	rm *.o

Aufgabe24.o: Aufgabe24.c
	gcc -c Aufgabe24.c $(Flags)

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c $(Flags)

run24: Aufgabe24
	./$< Aufgabe24


#Aufgabe25
Aufgabe25: Aufgabe25.o my_numerics.o
	gcc Aufgabe25.o my_numerics.o -o Aufgabe25 $(FLAGS)
	rm *.o

Aufgabe25.o: Aufgabe25.c
	gcc -c Aufgabe25.c $(Flags)

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -c ../my_numerics.c $(Flags)

run25: Aufgabe25
	./$< Aufgabe25

