#ifndef SIMULATION_H
#define SIMULATION_H

#include "parameters.hpp"
#include "ray.cuh"

// TODO do we need to include all these?
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

#include <curand.h>
#include <curand_kernel.h>

// TODO I made N a function argument for now, but we could just define it globally
// Passing it along is a little annoying I suppose

/// TODO I just made this a struct for now
/// Struct to hold and manage an array of rays
/// Assuming my_rays is initialized to some known length,
/// but you might only have my_size rays actually stored
struct RayGroup
{
    Ray* my_rays; // pointer to beginning of array
    int my_size;  // number of rays in array
};

/// Initialize doses (to zero)
__host__ void initialize_doses(double* doses, int N);

/// Initialize densities randomly between 0 and 1
__host__ void initialize_densities_random(double* densities, int N);

/// Initialize densities to constant value
__host__ void initialize_densities_constant(double* densities, int N, double density);

/// Initialize densities to 2D Gaussian centered on the grid,
//  with spread (std dev) given as fraction of grid size
__host__ void initialize_densities_centered_gaussian(double* densities, int N, double max_density, double spread);

/// Initialize densities with multiple Gaussians at random positions,
//  with highest resulting density normalized to given max_density
__host__ void
initialize_densities_random_gaussians(double* densities, int N, int n_gaussians, double max_density, double spread);

/// Write data
__host__ void write_to_csv_file(double* grid_data, int N, const std::string& filename);

/// Random distributions
std::default_random_engine random_engine{(uint_fast32_t)time(0)}; // seeded random number generator
std::uniform_real_distribution<double> uniform_dist{0.0, 1.0};    // uniform distribution, pass {lowerbound, upperbound}
std::uniform_real_distribution<double> uniform_angle_gen{0.0,
    2 * M_PI}; // uniform distributioin between 0 and 2 pi

__device__ void init_curand_state(curandState_t* state); // initialize curand state based on clock, thread id

__device__ double uniform_angle_dist(curandState_t* state); // uniform between 0 and 2 pi

__device__ double normal_dist(curandState_t* state, double mean, double std_dev); // normal with given mean and std dev

/// Randomly sample source angle for primary rays
__device__ double random_source_angle_normal();
__host__ double random_source_angle();

/// Randomly sample a coordinate for primary rays
__host__ Pixel random_pixel();

/// Check if a given pixel lies outside the bounds of the grid
__device__ bool out_of_bounds(Pixel current_pixel, int N);

/// Generate new primary ray from source
__host__ void spawn_primary_rays(Ray* rays, int N, int num_primary_rays);

/// Generate secondary rays from interaction point
__device__ void spawn_secondary_rays(RayGroup* group, Pixel spawn_pixel, double total_energy, int N);

/// Determine whether primary ray interacts at visited pixel based on distance travelled
__device__ bool random_interact(Pixel target_pixel, double distance, double* densities, int N);

/// Transfer energy from ray to target pixel,
//  where the actual energy transferred is scaled by the pixel density
__device__ void
transfer_energy(Ray* ray, Pixel target_pixel, double unscaled_energy, double* densities, double* doses, int N);

/// Evolve all active rays by one step, return the number of rays evolved
__device__ int evolve_rays(RayGroup* group, double* densities, double* doses, int N);

/// Evolve all rays until complete, i.e. none are active
__device__ void evolve_to_completion(RayGroup* group, double* densities, double* doses, int N);

/// Run simulation for a given number  primary rays, in serial
__device__ void run_serial(int num_primary_rays, double* densities, double* doses, int N);

#endif
