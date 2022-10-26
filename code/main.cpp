#include <iostream>
#include <armadillo>
#include "PenningTrap.hpp"
#include "Particle.hpp"

int main(){

    // Initialize the particle
    double q = 1;
    double m = 40.08;

    arma::vec r_initial = {20., 0., 20.};
    arma::vec v_initial = {0., 25., 0.};

    arma::vec r2_initial = {25., 25., 0.};
    arma::vec v2_initial = {0., 40., 5.};

    Particle p1(q, m, r_initial, v_initial);
    Particle p2(q, m, r2_initial, v2_initial);

    // Initialize the trap
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500.;

    std::vector<Particle> particles_in{p1};

    PenningTrap trap(particles_in, B0, V0, d);
    // trap.interaction = true;

    // Initialize the time and position 
    std::vector<int> n{4000, 8000, 16000, 32000};
    for (int nk : n){
        double dt = 50. / nk;
        for (int i=0; i <= 1; i++){
            trap.solve(dt, 50, true, i);
        }
    }

    // ex8, uncomment bellow
    // particles_in.push_back(p2);
    // PenningTrap ex8(particles_in, B0, V0, d);

    // double dt = 50. / 32000;
    // ex8.solve(dt, 50, true, 0);
    // ex8.solve(dt, 50, true, 1);
    return 0;
}
