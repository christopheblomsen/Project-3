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

all: compile link clean

Particle:
	g++ -c -I include src/Particle.cpp

PenningTrap:
	g++ -c -I include src/PenningTrap.cpp

main:
	g++ -c -I include main.cpp

ex9:
	g++ -c -I include ex9.cpp

give_me_omega:
	g++ -c -I include give_me_omega.cpp

compile: Particle PenningTrap main ex9 

link:
	g++ -o ex9 Particle.o PenningTrap.o ex9.o -larmadillo -O2
	g++ -o main Particle.o PenningTrap.o main.o -larmadillo -O2

run:
	./ex9

clean:
	@rm -rf *.o
# end
