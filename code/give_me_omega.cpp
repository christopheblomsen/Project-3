#include "PenningTrap.hpp"
#include <armadillo>

int main(){

// Initialize the particle
    double q = 1;
    double m = 40.08;

    arma::vec r_initial = {20., 0., 20.};
    arma::vec v_initial = {0., 25., 0.};

    // arma::vec r2_initial = {25., 25., 0.};
    // arma::vec v2_initial = {0., 40., 5.};

    Particle p1(q, m, r_initial, v_initial);
    // Particle p2(q, m, r2_initial, v2_initial);

    // Initialize the trap
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500.;

    std::vector<Particle> particles_in{p1};

    PenningTrap trap(particles_in, B0, V0, d);
    arma::vec omega_pm = trap.omega_pm(q, m);
    std::cout << omega_pm << std::endl;

}
