#ifndef SIMULATION_H
#define SIMULATION_H

#include "parameters.hpp"
#include "ray.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <thread>
#include <time.h>
#include <tuple>
#include <utility>
#include <vector>

class Simulation
{
public:
    Simulation() = delete;

    /// Constructor
    Simulation(const int N);

    /// Initialize densities randomly between 0 and 1
    void initialize_densities_random();

    /// Initialize densities to constant value
    void initialize_densities_constant(const double density);

    /// Initialize densities to 2D Gaussian centered on the grid,
    //  with spread (std dev) given as fraction of grid size
    void initialize_densities_centered_gaussian(const double max_density, const double spread);

    /// Initialize densities with multiple Gaussians at random positions,
    //  with highest resulting density normalized to given max_density
    void initialize_densities_random_gaussians(const int n_gaussians, const double max_density, const double spread);

    /// Run simulation for a given number of primary rays, in serial
    void run_serial(const int num_primary_rays);

    /// Write data
    void write_densities_to_file(const std::string& filename);
    void write_doses_to_file(const std::string& filename);

    /// Random distributions
    std::default_random_engine random_engine{(uint_fast32_t)time(0)}; // seeded random number generator
    std::uniform_real_distribution<double> uniform_dist{0.0,
                                                        1.0}; // uniform distribution, pass {lowerbound, upperbound}
    std::uniform_real_distribution<double> uniform_angle_dist{0.0,
                                                              2 * M_PI}; // uniform distributioin between 0 and 2 pi
    std::normal_distribution<double> normal_dist{PARAM_MEAN,
                                                 PARAM_SIGMA}; // normal distribribution, pass {mean, stddev}

private:
    const int m_N;                                // number of pixels per side
    std::vector<std::vector<double>> m_densities; // grid of density values
    std::vector<std::vector<double>> m_doses;     // grid of dose values
    std::vector<Ray> m_rays;                      // active rays

    /// Randomly sample source angle for primary rays
    double _random_source_angle(bool normal);

    /// Check if a given pixel lies outside the bounds of the grid
    bool _out_of_bounds(const PIXEL& current_pixel);

    /// Generate new primary ray from source
    void _spawn_primary_ray();

    /// Generate secondary rays from interaction point
    void _spawn_secondary_rays(PIXEL spawn_pixel, double total_energy);

    /// Determine whether primary ray interacts at visited pixel based on distance travelled
    bool _random_interact(PIXEL target_pixel, double distance);

    /// Transfer energy from ray to target pixel,
    //  where the actual energy transferred is scaled by the pixel density
    void _transfer_energy(Ray* ray, PIXEL target_pixel, double unscaled_energy);

    /// Evolve all active rays by one step, return the number of rays evolved
    int _evolve_rays();

    /// Evolve all rays until complete, i.e. none are active
    void _evolve_to_completion();
};

#endif
