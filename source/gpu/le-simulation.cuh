#ifndef SIMULATION_H
#define SIMULATION_H

#include "../parameters.hpp"
#include "ray.cuh"

// TODO do we need to include all these?
#include <chrono>
#include <cmath>
#include <curand.h>
#include <curand_kernel.h>
#include <fstream>
#include <iostream>
#include <random>
#include <thread>
#include <time.h>
//#include <tuple>
//#include <utility>
#include <vector>

#define error_check(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


////////// STRUCTS AND TYPEDEFS //////////

/// Struct to hold and manage an array of rays
struct RayGroup
{
    Ray* my_rays; // pointer to beginning of array
    int my_size;  // number of rays in array
    int max_size; // max number of rays allocated
};

/// Shorthand for a vector of ray groups
typedef std::vector<RayGroup> RegionGroup;

/// Struct to hold coordinates of a region within the grid
struct Region
{
    int first;
    int second;
};

/// Struct to handle buffer for regrouping rays
//  The arrays here consist of p sections of size q
//  In use, each thread writes to its own section
struct RegroupBuffer
{
    int section_size;    // this is q, how much space is allocated to each section
    Ray* rays;           // size p*q, rays to be regrouped
    int* region_indices; // size p*q, linear indices of regions to which to move rays
    int* ray_counts;     // size p, counter of how many rays are in each section of buffer
};

/// Struct to handle buffer for depositing energy into pixels
// I think a key difference with this buffer is that there's no data dependencies for the next regions in the schedule
// so when GPU finishes a region it will fill this buffer and CPU can write it's values while GPU goes to next region
// essentially CPU is stagger writing
struct DoseBuffer
{
    int section_size; // how much space is allocated to each section. max number of pixels a ray can travel in a region
                      // (this is the number of pixels in the region diagonal * number of rays right?)
    double* energy_deposited; // energy deposited to a pixel
    Pixel* pixel_traversed;   // pixel traversed by single ray where energy was deposited
    int* pixel_counts;        // total count of pixels traversed by all rays in a ray group
};

////////// GRID INITIALIZATION //////////

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

////////// OUTPUT //////////

/// Write data
__host__ void write_to_csv_file(double* grid_data, int N, const std::string& filename);

////////// RANDOMIZATION //////////

/// Random distributions
std::default_random_engine random_engine{(uint_fast32_t)time(0)}; // seeded random number generator
std::uniform_real_distribution<double> uniform_dist{0.0, 1.0};    // uniform distribution, pass {lowerbound, upperbound}
std::normal_distribution<double> normal_dist{PARAM_MEAN, PARAM_SIGMA}; // normal distribribution, pass {mean, stddev}

__device__ void init_curand_state(curandState_t* state);    // initialize curand state based on clock, thread id
__device__ double uniform_angle_dist(curandState_t* state); // uniform between 0 and 2 pi
__device__ double normal_angle_dist(curandState_t* state, double mean, double std_dev); // normal distribution

/// Randomly sample source angle for primary and secondary rays
__device__ double random_source_angle_uniform();
__host__ double random_source_angle_normal();

////////// RAY LOCATION CHECKING //////////

/// Return the region containing a give pixel
__host__ __device__ Region get_region(Pixel position, int N, int M);

/// Return the linear index of the region containing a given pixel
__host__ __device__ int get_region_index(Pixel position, int N, int M);

/// Check if a given pixel lies outside the bounds of the grid
__device__ bool out_of_bounds(Pixel current_pixel, int N);

////////// RAY CREATION //////////

/// Generate primary rays from source and add as ray groups within appropriate region groups
__host__ void spawn_primary_rays(
    std::vector<RegionGroup>& region_groups, int num_primary_rays, int max_rays_per_ray_group, int N, int M);

/// Generate secondary rays from given interact point and add them to ray group
__device__ void spawn_secondary_rays(RayGroup* group, Pixel spawn_pixel, double total_energy, int N);

////////// THREAD GROUP EVOLUTION //////////

/// Determine whether primary ray interacts at visited pixel based on distance travelled
__device__ bool random_interact(Pixel target_pixel, double distance, double* densities, int N);

/// Transfer energy from ray to target pixel,
//  where the actual energy transferred is scaled by the pixel density
__device__ void
transfer_energy(Ray* ray, Pixel target_pixel, double unscaled_energy, double* densities, double* doses, int N);

/// Evolve all active rays in group by one step, return the number of rays evolved
__device__ int
evolve_rays(RayGroup* group, int region_index, double* densities, double* doses, int N, int M, RegroupBuffer* g_buffer_cuda, int num_ray_groups);

/// Evolve all rays in group until complete, i.e. none are active
__device__ void evolve_to_completion(
    RayGroup* group, int region_index, double* densities, double* doses, int N, int M, RegroupBuffer* g_buffer_cuda, int num_ray_groups);

/// Kernel function: Run each thread group in the given region group array in parallel,
//  rays within each thread group are run in serial
__global__ void run_rays(RayGroup* region_group_arr,
                         int region_group_arr_size,
                         int region_index,
                         double* densities,
                         double* doses,
                         int N,
                         int M,
                         RegroupBuffer* g_buffer_cuda);

////////// REGION GROUP RUNNING AND PROCESSING //////////
// Use data in regroup buffer to add vectors to correct region group
__host__ void
regroup(std::vector<RegionGroup>& region_groups, RegroupBuffer* g_buffer, int num_ray_groups);

// allocate regroup buffer on device
__host__ void init_regroup_buffer_cuda(RegroupBuffer* &g_buffer_cuda, int max_num_rays, int num_ray_groups);

// allocate regroup buffer on host and copy device's regroup buffer to host's version
__host__ void
copy_regroup_buffer_host(RegroupBuffer* &g_buffer, RegroupBuffer* &g_buffer_cuda, int max_num_rays, int num_ray_groups);

/// Get list of linear indices for traversing L x L array of tasks in diagonal order
//  e.g. [0,1,3,2,4,6,5,7,8] for 3 x 3
__host__ std::vector<int> get_forward_schedule(int L);

/// Run region groups, one at a time, in back-and-forth diagonal fashion until all rays are done
__host__ void
run_region_groups(std::vector<RegionGroup>& region_groups, double* densities, double* doses, int N, int M);

/// Run each ray group within a single region group on its own thread
__host__ void run_region_group(RegionGroup& region_group,
                               int region_index,
                               double* densities,
                               double* doses,
                               int N,
                               int M,
                               RegroupBuffer* &g_buffer);

#endif
