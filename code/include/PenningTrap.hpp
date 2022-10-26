#ifndef PENNINGTRAP_H_
#define PENNINGTRAP_H_

#include "Particle.hpp"
#include <armadillo>
#include <string>
#include <vector>

struct PenningTrap{

        // Public variables
        double B0;  // Magnetic field strength
        double V0;  // Voltage
        double V02;
        double mass_ion;
        double charge_ion;
        double d;   // Distance between plates
        double ke = 1.3893533e5; // Coulomb constant
        std::vector<Particle> particles;   // Vector of particles
        bool interaction;   // Interaction between particles
        int par_left;

        arma::vec N_array;
        arma::vec t;
        arma::cube v;
        arma::cube r;

        // For EX9 not useful for now
        double frequency;   // Frequency of the oscillation
        double amplitude;   // Amplitude of the oscillation

        // Constructor
        PenningTrap(std::vector<Particle> particles, double B0_in, double V0_in, double d_in);

        // Constructor for frequency and amplitude
        PenningTrap(std::vector<Particle> particles, double B0_in, double V0_in, double d_in,
                        double omega, double f);

        //return the analytical freqencies for a single particle
        arma::vec omega_pm(double q, double m);

        // Add a particle to the trap
        void add_particle(Particle p_in);

        // Add several particles to the trap
        void add_several_particles(std::vector<Particle> p_ins);

        // Checks if particle norm is outside the dim of trap
        void particles_left();

        // Number of particle trapped
        int number_of_particles();

        // to add N particles at once to the trap
        void add_particles(int n_par);

        // calculate the external E_field without time dependence
        arma::vec external_E_field(arma::vec r);

        // caculate the external E field with time dependence
        arma::vec external_E_field(arma::vec r, double t);

        //calculate the external B_field without time dependence
        arma::vec external_B_field(arma::vec r);

        // Force on particle_i from particle_j
        arma::vec force_particle(int i, int j);

        // The total force on particle_i from the external fields Lorentz force
        arma::vec total_force_external(int i);

        // The total force on particle_i from the external fields Lorentz force
        // time dependent
        arma::vec total_force_external(int i, double t);

        // The total force on particle_i from the other particles
        arma::vec total_force_particles(int i);

        // The total force on particle_i from both external fields and other particles
        // time independent
        arma::vec total_force(int i);

        // The total force on particle_i from both external fields and other particles
        // time dependent
        arma::vec total_force(int i, double t);

        // Solves the system for all time steps
        void solve(double dt, double total_time, bool interaction_in, int model);

        // solves the system for all time step, for ex9
        void solve_with(double dt, double total_time, bool interaction_in, int model);

        // evolves one step with forward Euler
        void evolve_forward_Euler(double dt, int j, double t);

        // evolves one step with Runge-Kutta 4
        void evolve_RK4(double dt, int j, double t);
};

#endif // PENNINGTRAP_H_
