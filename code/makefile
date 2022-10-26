##
# Project 3
# To run simply write
#
# make
#
# if you want to just do one of them
#
# make foo
#
# where foo is compile, link, run, etc.
# @file
# @version 0.1

all: compile link run clean

Particle:
	g++ -c -I include src/Particle.cpp

PenningTrap:
	g++ -c -I include src/PenningTrap.cpp

full:
	g++ -c -I include ex9.cpp

ex8:
	g++ -c -I include main.cpp

compile: Particle PenningTrap ex8

link:
	g++ -o main Particle.o PenningTrap.o main.o -larmadillo -O2

run:
	./main

clean:
	@rm -rf *.o
# end
