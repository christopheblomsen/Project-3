#ifndef ODESOLVER_H_
#define ODESOLVER_H_
#include "PenningTrap.hpp"
#include "Particle.hpp"

#include <armadillo>
#include <vector>

struct ODESolver{
    arma::vec f_;
    int N;
    arma::cube r, v;

    // Constructor
    ODESolver(PenningTrap penningtrap);
    std::vector<Particle> particles_;
    PenningTrap penning_trap_;
    bool interaction_;

    void step();
    // Solves the system with given model time and timestep
    void solve(double t0, double t_end, double dt, int model, bool interaction);

    // i is the particle
    // j is the "time step"
    // Runge-Kutta4 step for solve
    void RK4_step(double dt, int i, int j, double t);

    // Forward-Euler step for solve
    void FE_step(double dt, int i, int j, double t);
};

// class ForwardEuler: public ODESolver{
//         public:
//             arma::vec advance();
// };

// class RungeKutta4: public ODESolver{
//         public:
//             arma::vec advance();
// };


#endif // ODESOLVER_H_
