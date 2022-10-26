#include <iostream>
#include <armadillo>
#include "PenningTrap.hpp"
#include "Particle.hpp"

void fine_grain(int tot_particle, double tot_time, arma::vec omegas, arma::vec fraction_left, 
				double amps, int N, bool interaction){
	// Tar inn antall partikler, total tid, omega array og frequency
    // Initialize the particle
    double q = 1;
    double m = 40.08;

    arma::vec r_initial = {20., 0., 20.};
    arma::vec v_initial = {0., 25., 0.};

    // arma::vec r2_initial = {25., 25., 0.};
    // arma::vec v2_initial = {0., 40., 5.};

    Particle p1(q, m, r_initial, v_initial);
	int add_par = tot_particle - 1;
    // Particle p2(q, m, r2_initial, v2_initial);

    // Initialize the trap
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500.;

    std::vector<Particle> particles_in{p1};
    double dt = tot_time/N;
	int i = 0;
	for ( double omega : omegas ){
		PenningTrap trap(particles_in, B0, V0, d, omega, amps);
		trap.add_particles(add_par);
		trap.solve_with(dt, tot_time, interaction, 0);
		
		int particles_left = trap.par_left;
		double fraction_left_ = particles_left / tot_particle;
		fraction_left(i) = fraction_left_;
		i++;
	}
	
	std::string A_s = std::to_string(amps);
	
	if(interaction){
		fraction_left.save("fraction_left_fine_grain_f_" + A_s + "_interaction.bin");
		omegas.save("omega_fine_grain_f_" + A_s + "_interaction.bin");
	}
	else{
		fraction_left.save("fraction_left_fine_grain_f_" + A_s + "_no_interaction_.bin");
		omegas.save("omega_fine_grain_f_" + A_s + "_no_interaction_.bin");
	}
}

void part_one(int tot_particle, double tot_time, arma::vec omegas, arma::vec fraction_left, 
				std::vector<double> amps, int N){
	// Initialize the particle
    double q = 1;
    double m = 40.08;

    arma::vec r_initial = {20., 0., 20.};
    arma::vec v_initial = {0., 25., 0.};

    // arma::vec r2_initial = {25., 25., 0.};
    // arma::vec v2_initial = {0., 40., 5.};

    Particle p1(q, m, r_initial, v_initial);
	int add_par = tot_particle - 1;
    // Particle p2(q, m, r2_initial, v2_initial);

    // Initialize the trap
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500.;

    std::vector<Particle> particles_in{p1};
    double dt = tot_time/N; // should be 500 ms


	int i;
    for ( double A : amps ){
		i = 0;
        for ( double omega : omegas ){
            PenningTrap trap(particles_in, B0, V0, d, omega, A);
            trap.add_particles(add_par);
            trap.solve_with(dt, tot_time, false, 0); // no interaction, Runge Kutta 4
			
			int particles_left = trap.par_left;
			double fraction_left_ = particles_left / tot_particle;
			fraction_left(i) = fraction_left_;
			i++;
        }
		std::string A_s = std::to_string(A);
		std::string N_s = std::to_string(N);
		fraction_left.save("fraction_left_f_" + A_s + "_N_" + N_s + ".bin");
		omegas.save("omega_f_" + A_s + "_N_" + N_s + ".bin");
    }
}

int main(){
	std::cout << "Starting:" << std::endl;
	
	// Variables for both simulations
	int tot_particle = 100; // total number of particles in the simulation
	double tot_time = 500.; // the total time in ms
	
	// // Part one of the simulation: Broad exploration of frequencies
	// arma::vec fraction_left = arma::linspace(0.2, 2.5, 115); // very hacky
	// arma::vec omegas = arma::linspace(0.2, 2.5, 115);
	// std::vector<double> amps{0.1, 0.4, 0.7};
	
	// // Calling for different timesteps N: 4000, 8000, 16000, 32000
	// int N_temp = 4000;
	// part_one(tot_particle, tot_time, omegas, fraction_left, amps, N_temp);
	// int N_temp_ = 8000;
	// part_one(tot_particle, tot_time, omegas, fraction_left, amps, N_temp_);
	// int N = 16000; // middle from task 8
	// part_one(tot_particle, tot_time, omegas, fraction_left, amps, N);
	// int N_ = 32000;
	// part_one(tot_particle, tot_time, omegas, fraction_left, amps, N_);
	
	// Part two of the simulation: Zoomed in on the intervall [2.0, 2.5] for f = 0.1
	int N = 32000;
	double amp{0.1}; // The amplitude for the potential
	arma::vec omega_arr = arma::linspace(2.0, 2.5, 230); // dt of 0.002, is enough?
	arma::vec fraction_left_arr = arma::linspace(2.0, 2.5, 230);; // very hacky
	
	// Calling with and without interactions
	fine_grain(tot_particle, tot_time, omega_arr, fraction_left_arr, amp, N, false); // interaction
	fine_grain(tot_particle, tot_time, omega_arr, fraction_left_arr, amp, N, true); // no interaction
    return 0;
}
