#ifndef SIMULATION_H
#define SIMULATION_H


#include "parameters.hpp"
#include "ray.hpp"

#include <vector>
#include <random>
#include <iostream>
#include <cmath>
#include <utility>
#include <tuple>


class SimulationSerial
{
    public: 
    // SimulationSerial() = delete;
    
    SimulationSerial(const int N, const double density, const int ray_count); /// will this work? for static arrays?


    /// Initalize densities
    void initialize_densities_random();
    void initialize_densities_constant(const double density);

    /// Run simulation for a given number of primary rays
    void run(int num_primary_rays);

    /// Get array of doses
    double* get_doses();

    // Print data
    void print_m_densities();
    void print_m_doses();

    std::default_random_engine random_engine {0};  // seeded random number generator 
    std::uniform_real_distribution<double> uniform_dist {0.0, 1.0}; // uniform distribution, pass {lowerbound, upperbound}
    std::normal_distribution<double> normal_dist {MEAN, SIGMA}; // normal distribribution, pass {mean, stddev}
    
    private:
    const int N;
    std::vector< std::vector<double>> m_densities;      /// pixel density values
    std::vector< std::vector<double>> m_doses;                   /// pixel dose values
    //Ray test;
    std::vector<Ray> m_rays;                /// active rays

    /// Randomly sample source angle for primary rays
    double _random_source_angle();

    /// Generate new primary ray
    void _spawn_primary_ray();

    /// Evolve all rays by one step
    void _evolve_rays();

    /// Evolve all rays until complete
    void _evolve_to_completion();

    

};


#endif