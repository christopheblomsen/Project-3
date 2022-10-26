#include "Particle.hpp"

Particle::Particle(double q, double m, arma::vec r, arma::vec v){
    charge = q;
    mass = m;
    position = r;
    velocity = v;
}

void Particle::get_charge(){
    std::cout << charge << std::endl;
}

void Particle::get_mass(){
    std::cout << mass << std::endl;
}

void Particle::get_position(){
    std::cout << position << std::endl;
}

void Particle::get_velocity(){
    std::cout << velocity << std::endl;
}

void Particle::get_info(){
    get_charge();
    get_mass();
    get_position();
    get_velocity();
}
