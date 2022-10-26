#include "PenningTrap.hpp"
#include "Particle.hpp" 
  
// Constructor
PenningTrap::PenningTrap(std::vector<Particle> particles_in, double B0_in, double V0_in, double d_in){
  particles = particles_in;
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
  V02 = V0 / (d * d);

  // The mass is the same for all ions
  mass_ion = particles[0].mass;

  // The charge is the same for all ions
  charge_ion = particles[0].charge;
}

// Constructor for frequency and amplitude
PenningTrap::PenningTrap(std::vector<Particle> particles_in, double B0_in, double V0_in, double d_in,
                         double omega, double f){
  particles = particles_in;
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
  V02 = V0 / (d * d);

  // The mass is the same for all ions
  mass_ion = particles[0].mass;

  // The charge is the same for all ions
  charge_ion = particles[0].charge;

  frequency = omega;
  amplitude = f;

}

//return the analytical freqencies for a single particle
arma::vec PenningTrap::omega_pm(double q, double m){
  arma::vec omega_pm = {0, 0};
  double omega_0 = q*B0/m;
  double omega_z2 = 2*q*V0/(m*d*d);
  double omega_p = omega_0/2 + std::sqrt(std::pow(omega_0, 2) - 2*omega_z2)/2;
  double omega_n =  omega_0/2 - std::sqrt(std::pow(omega_0, 2) - 2*omega_z2)/2;
  omega_pm = {omega_p, omega_n};
return omega_pm;
}


// Checks if particle norm is outside the dim of trap
void PenningTrap::particles_left(){
  par_left = 0;
  for (int i=0; i < particles.size(); i++){
    if (arma::norm(particles[i].position) <= d){
      par_left += 1;
    }
  }
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in){
  PenningTrap::particles.push_back(p_in);
}

// Adds several particles to the trap
void PenningTrap::add_several_particles(std::vector<Particle> p_ins){
        particles.insert(particles.end(), p_ins.begin(), p_ins.end());
        particles_left();
}

// to add N particles at once to the trap
void PenningTrap::add_particles(int n_par){

  for (int i = 0; i < n_par; i++){
    arma::vec r = arma::vec(3).randn() * 0.1 * d;
    arma::vec v = arma::vec(3).randn() * 0.1 * d;

    Particle particle = Particle(1, 40.08, r, v);
    add_particle(particle);
  }
}

// number of particles in trap
int PenningTrap::number_of_particles(){
  int N = 0;
  for (int i = 0; i < particles.size(); i++){
    if (arma::norm(particles[i].position) < d){
      N += 1;
    }
  }
  return N;
}

// External electric field at point r=(x,y,z)
// time independent
arma::vec PenningTrap::external_E_field(arma::vec r){
  arma::vec E = {0, 0, 0};

  if (arma::norm(r) < d){  
    E(0) = r(0) * V02;
    E(1) = r(1) * V02;
    E(2) = -2. * r(2) * V02;
  }
  return E;
}

// E field with time dependence
arma::vec PenningTrap::external_E_field(arma::vec r, double t){
  arma::vec E = {0, 0, 0};

  if (arma::norm(r) < d){
    E(0) = r(0);
    E(1) = r(1);
    E(2) = - r(2);

    E = E*V02*(1 + amplitude*cos(frequency*t));
  }

  return E;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r){
  arma::vec B = {0, 0, 0};

  if (arma::norm(r) < d){
    return {0, 0, B0};
  }
  return B;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j){
  // Calculate the distance between the particles and its norm
  arma::vec r_ij = particles[i].position - particles[j].position;
  double r = arma::norm(r_ij);

  // The charge of each particle
  double q_i = particles[i].charge;
  double q_j = particles[j].charge;
  // In this exact solution they could be
  // double q_i, q_j = charge_ion;

  // Initialize the force 
  arma::vec F = {0, 0, 0};
  
  // Check if the particles are not the same maybe it is the problem
  if (i != j){
    F = ke * (q_i * q_j * r_ij) / (r * r * r);
  }
  
  return F;

}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
  
  arma::vec F = {0, 0, 0};
  // Calculate the force on particle i from the external fields
  F = particles[i].charge *
      (external_E_field(particles[i].position)
      + arma::cross(particles[i].velocity, external_B_field(particles[i].position)));
  
  return F;

}

// The total force on particle_i from the external fields with time dependence
arma::vec PenningTrap::total_force_external(int i, double t){
  arma::vec F = {0, 0, 0};
  arma::vec ri = particles[i].position;
  arma::vec vi = particles[i].velocity;
  double qi = particles[i].charge;
  // For this we could've also used
  // double qi = charge_ion;

  arma::vec E = external_E_field(ri, t);
  arma::vec B = external_B_field(ri);

  // Calculate the force on particle i from the external fields
  F = qi * (E + arma::cross(vi, B));

  return F;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){
  arma::vec F = {0, 0, 0};
  // Calculate the force on particle i from the other particles
  for (int j = 0; j < particles.size(); j++){
    // We check if the particle j is far out of the trap
    if (i != j){
      F += force_particle(i, j);
    }
  }
  return F;
}

// The total force on particle_i from both external fields and other particles
// time independent
arma::vec PenningTrap::total_force(int i){

  arma::vec F = {0, 0, 0};
  // Calculate the total force on particle i
  F = total_force_external(i);

  if(interaction){
    F += total_force_particles(i);
  }
  return F;
}

// The total force on particle_i from both external fields and other particles
// time dependent
arma::vec PenningTrap::total_force(int i, double t){
  arma::vec F = {0, 0, 0};
  F = total_force_external(i, t);

  if(interaction){
    F += total_force_particles(i);
  }

  return F;
}

// solve method for ex9
void PenningTrap::solve_with(double dt, double total_time, bool interaction_in, int model){
  int n = (int) (total_time/dt);
  int n_par = particles.size();
  std::string method;

  interaction = interaction_in;

  v = arma::cube(3, n_par, n).fill(0.);
  r = arma::cube(3, n_par, n).fill(0.);
  t = arma::linspace(0, total_time, n);
  for (int j=0; j < n-1; j++){
      switch(model){
        case 0:
          evolve_RK4(dt, j, t(j));
          method = "RK4";
          break;
        case 1:
          evolve_forward_Euler(dt, j, t(j));
          method = "FE";
          break;
        default:
          std::cout << "No method " << std::endl;
    }
    particles_left();
  }
  particles_left();
}

// Solves the system for all time steps
void PenningTrap::solve(double dt, double total_time, bool interaction_in, int model){
  int n = (int) (total_time/dt);
  int n_par = particles.size();
  std::string method;
  std::string n_string = std::to_string(n);

  interaction = interaction_in;

  v = arma::cube(3, n_par, n).fill(0.);
  r = arma::cube(3, n_par, n).fill(0.);
  t = arma::linspace(0, total_time, n);

  for (int j=0; j < n-1; j++){
      switch(model){
        case 0:
          evolve_RK4(dt, j, t(j));
          method = "RK4";
          break;
        case 1:
          evolve_forward_Euler(dt, j, t(j));
          method = "FE";
          break;
        default:
          std::cout << "No method " << std::endl;
      }
    particles_left();
  }

  t.save("time.bin");
  r.save("position_"+method+"_"+n_string+".bin");
  v.save("velocity_"+method+"_"+n_string+".bin");
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, int j, double t){
  arma::vec a = arma::vec(3);
  arma::vec r_old, v_old, v_new, r_new, v_1, r_1, F;
  int N = particles.size();

  std::vector<arma::vec> r_vector(N);
  std::vector<arma::vec> v_vector(N);

  // Iterate over all particles
  for (int i = 0; i < N; i++){
    v_1 = particles[i].velocity;
    r_1 = particles[i].position;

    r.slice(j).col(i) = r_1;
    v.slice(j).col(i) = v_1;

    F = total_force(i, t);

    a = F/mass_ion;

    v_new = v_1 + a*dt;
    r_new = r_1 + v_1*dt;

    v.slice(j+1).col(i) = v_new;
    r.slice(j+1).col(i) = r_new;

    particles[i].position = r_new;
    particles[i].velocity = v_new;
  }
}

// RK4 with time dependence
void PenningTrap::evolve_RK4(double dt, int j, double t){
    // IMPORTANT
    // We have to loop separetly the coefficients 
    // because we need to use the old positions and velocities to calculate the new ones
    // when we have interaction since we need the old position at all time 
    // in the calculation of the force


    // Initialization
    int N = particles.size();

    // Create a initial 'pictures' of our particle 
    // Like this we can always access the intial position and velocity
    std::vector<Particle> particles_initial = particles;

    // "Push" all particles one time step
    evolve_forward_Euler(dt, j, t);
 
    // l is for the velocity and k for the position
    std::vector<arma::vec> k1(N), k2(N), k3(N), k4(N);
    std::vector<arma::vec> l1(N), l2(N), l3(N), l4(N);

    // Calculate the k1 and l1 for each particle
    for (int i = 0; i < N; i++){
      // Calculate k1 and l1
      l1[i] = dt * total_force(i, t) / mass_ion;
      k1[i] = dt * particles_initial[i].velocity;
    }

    // Calculates the corresponding (y_i + k1/2) in alg
    for (int i = 0; i < N; i++){
      //Calculate new r and v
      arma::vec r_new = particles_initial[i].position + k1[i] / 2.;
      arma::vec v_new = particles_initial[i].velocity + l1[i] / 2.;

      // Update the particle
      particles[i].position = r_new;
      particles[i].velocity = v_new;
    }

    // Calculate the k2 and l2 for each particle
    for (int i = 0; i < N; i++){
      // Calculate k2 and l2
      l2[i] = dt * total_force(i, t + 0.5 * dt) / mass_ion;
      k2[i] = dt * particles[i].velocity;
    }

    // Calculates the (y_i + k2/2) in alg
    for (int i = 0; i < N; i++){
      //Calculate new r and v
      arma::vec r_new = particles_initial[i].position + k2[i] / 2.;
      arma::vec v_new = particles_initial[i].velocity + l2[i] / 2.;

      // Update the particle
      particles[i].position = r_new;
      particles[i].velocity = v_new;
    }

    // Calculate the k3 and l3 for each particle
    for (int i = 0; i < N; i++){
      // Calculate k3 and l3
      l3[i] = dt * total_force(i, t + 0.5*dt) / mass_ion;
      k3[i] = dt * particles[i].velocity;
    }

    for (int i = 0; i < N; i++){
      //Calculate new r and v
      arma::vec r_new = particles_initial[i].position + k3[i];
      arma::vec v_new = particles_initial[i].velocity + l3[i];

      // Update the particle
      particles[i].position = r_new;
      particles[i].velocity = v_new;
    }

    // Calculate the k4 and l4 for each particle
    for (int i = 0; i < N; i++){
      // Calculate k4 and l4
      l4[i] = dt * total_force(i, t + dt) / mass_ion;
      k4[i] = dt * particles[i].velocity;
    }

    for (int i = 0; i < N; i++){
      //Calculate new r and v
      arma::vec r_new = particles_initial[i].position + (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]) / 6.;
      arma::vec v_new = particles_initial[i].velocity + (l1[i] + 2. * l2[i] + 2. * l3[i] + l4[i]) / 6.;

      // Update the particle
      particles[i].position = r_new;
      particles[i].velocity = v_new;

      r.slice(j+1).col(i) = r_new;
      v.slice(j+1).col(i) = v_new;
    }
}
