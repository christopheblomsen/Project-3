#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <armadillo>
#include <string>
#include <vector>

// Structure
struct Particle{

    // Definition of public so that we can access the variables from outside the class (?)
    double charge;
    double mass;
    arma::vec position;
    arma::vec velocity;

    Particle(double q, double m, arma::vec r, arma::vec v);

    // Returns the charge of the particle
    void get_charge();

    // Returns the mass of the particle
    void get_mass();

    // Returns the position of the particle
    void get_position();

    // Returns the velocity of the particle
    void get_velocity();

    // gets all info
    void get_info();
};

#endif // PARTICLE_H_
