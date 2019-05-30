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

class SimulationSerial
{
public:
    SimulationSerial() = delete;
    SimulationSerial(const int N);

    /// Initalize densities
    void initialize_densities_random();
    void initialize_densities_constant(const double density);

    /// Run simulation for a given number of primary rays
    void run(int num_primary_rays);

    /// Get array of doses
    /// std::vector< std::vector<double>> get_doses();

    /// Print data
    void print_densities();
    void print_doses();
    void write_to_file();

    /// Random distributions
    std::default_random_engine random_engine{(uint_fast32_t)time(0)}; // seeded random number generator
    std::uniform_real_distribution<double> uniform_dist{0.0,
                                                        1.0}; // uniform distribution, pass {lowerbound, upperbound}
    std::uniform_real_distribution<double> uniform_angle_dist{0.0,
                                                              2 * M_PI}; // uniform distributioin between 0 and 2 pi
    std::normal_distribution<double> normal_dist{MEAN, SIGMA};           // normal distribribution, pass {mean, stddev}

private:
    const int m_N; /// number of pixels per side

    // JK: is it ok to store these as std vectors rather than arrays? for GPU impl?
    std::vector<std::vector<double>> m_densities; /// grid of density values
    std::vector<std::vector<double>> m_doses;     /// grid of dose values
    std::vector<Ray> m_rays;                      /// active rays

    /// Randomly sample source angle for primary rays
    double _random_source_angle(bool normal);

    /// Generate new primary ray
    void _spawn_primary_ray();

    /// Evolve all rays by one step
    int _evolve_rays();

    /// Evolve all rays until complete
    void _evolve_to_completion();

    /// Determine whether primary ray interacts at current pixel
    bool _random_interact(Ray* r, PIXEL visited, double distance);

    /// Generate secondary rays from primary ray
    std::vector<Ray> _spawn_secondary_rays(Ray* primary); // TODO no ray argument

    /// Enforce limits: 0 <= angle < 2*pi
    double _normalize_angle(double angle);

    /// Deposit energy from ray to pixel visited
    void _deposit_energy(Ray* r, PIXEL visited, double distance);

    /// Fixes position discrepency when spawning secondary rays that are going in the opposite direction of the primary
    /// ray Returns corrections to current pixel
    PIXEL _fix_position(PIXEL_EDGE edge, double current_angle, double new_angle);
};

#endif
