#include "base-simulation.cuh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

// Needed to move this to .cu
std::default_random_engine random_engine{(uint_fast32_t)time(0)}; // seeded random number generator
std::uniform_real_distribution<double> uniform_dist{0.0, 1.0};    // uniform distribution, pass {lowerbound, upperbound}

__host__ void initialize_doses(double* doses, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            doses[i * N + j] = 0;
        }
    }
    return;
}

__host__ void initialize_densities_random(double* densities, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            densities[i * N + j] = uniform_dist(random_engine);
        }
    }
    return;
}

__host__ void initialize_densities_constant(double* densities, int N, double density)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            densities[i * N + j] = density;
        }
    }
    return;
}

__host__ void initialize_densities_centered_gaussian(double* densities, int N, double max_density, double spread)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int mid = N / 2; // middle pixel
            int x = i - mid;
            int y = j - mid;
            double std_dev = spread * N;
            densities[i * N + j] = max_density * exp(-(x * x + y * y) / (2 * std_dev * std_dev));
        }
    }
    return;
}

__host__ void
initialize_densities_random_gaussians(double* densities, int N, int n_gaussians, double max_density, double spread)
{
    double std_dev = spread * N;
    double highest = 0;

    // Add the Gaussians
    for (int k = 0; k < n_gaussians; k++)
    {
        int mid_x = floor(uniform_dist(random_engine) * N);
        int mid_y = floor(uniform_dist(random_engine) * N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                int x = i - mid_x;
                int y = j - mid_y;
                densities[i * N + j] += max_density * exp(-(x * x + y * y) / (2 * std_dev * std_dev));
                highest = fmax(highest, densities[i * N + j]);
            }
        }
    }

    // Normalize the resulting density distribution
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            densities[i * N + j] = max_density * densities[i * N + j] / highest;
        }
    }
    return;
}

__host__ void write_to_csv_file(double* grid_data, int N, const std::string& filename)
{
    std::ofstream output;
    output.open(filename);
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N - 1; i++)
        {
            output << grid_data[i * N + j] << ",";
        }
        output << grid_data[N - 1 + j] << "\n";
    }
    output.close();
    return;
}

__device__ void init_curand_state(curandState_t* state)
{
    // Initialize random kernel
    int tId = threadIdx.x + (blockIdx.x * blockDim.x);
    curand_init((unsigned long long)clock(), tId, 0, state);

    return;
}

__device__ double uniform_angle_dist(curandState_t* state) { return 2 * M_PI * curand_uniform_double(state); }

__device__ double normal_dist(curandState_t* state, double mean, double std_dev)
{
    return mean + std_dev * curand_normal_double(state);
}
__device__ double random_source_angle(bool normal)
{
    curandState state;
    init_curand_state(&state);

    double angle;
    if (normal)
    { // Normal distribution
        angle = normal_dist(&state, PARAM_MEAN, PARAM_SIGMA);
    }
    else
    { // Uniform distribution between 0 and 2 pi
        angle = uniform_angle_dist(&state);
    }

    // Normalize angle
    while (angle < 0.0)
    {
        angle += 2 * M_PI;
    }
    while (angle >= 2 * M_PI)
    {
        angle -= 2 * M_PI;
    }
    return angle;
}

__device__ bool out_of_bounds(Pixel current_pixel, int N)
{
    return (current_pixel.first < 0 || current_pixel.first >= N || current_pixel.second < 0 ||
            current_pixel.second >= N);
}

__device__ void spawn_primary_ray(RayGroup* group, int N)
{
    // Randomly select source angle from normal distribution
    double source_angle = random_source_angle(true);

    //printf("block %d, thread %d spawning primary ray with angle %0.5f\n", blockIdx.x, threadIdx.x, source_angle);

    // Calculate initial ray position
    double horiz_dist_from_center = PARAM_D * N * tan(source_angle); // horizontal distance from center of top edge
    int middle_pixel = N / 2;
    double horiz_dist_from_left = middle_pixel + horiz_dist_from_center;

    // If ray does not miss grid entirely, spawn it
    if (horiz_dist_from_left < 0 || horiz_dist_from_left >= N ||
        (source_angle >= M_PI / 2 && source_angle <= 3 * M_PI / 2))
    {
        return;
    }
    else
    {
        double horiz_dist_from_left_rounded = floor(horiz_dist_from_left);
        Pixel spawn_pixel;
        spawn_pixel.first = horiz_dist_from_left_rounded;
        spawn_pixel.second = 0; // always starts from top of grid
        double edge_dist = horiz_dist_from_left - horiz_dist_from_left_rounded;
        group->my_rays[group->my_size] = Ray::primary(source_angle, spawn_pixel, PIXEL_EDGE::TOP, edge_dist);
        group->my_size++;
    }
    return;
}

__device__ void spawn_secondary_rays(RayGroup* group, Pixel spawn_pixel, double total_energy, int N)
{
    for (int i = 0; i < PARAM_KS; i++)
    {
        double source_angle = random_source_angle(false); // uniform random source angle
        double partial_energy = total_energy / PARAM_KS;
        Ray new_ray = Ray::secondary_from_center(source_angle, spawn_pixel, partial_energy);
        Pixel current_pixel = new_ray.get_current_pixel();
        if (out_of_bounds(current_pixel, N))
        {
            // Ray is out of bounds, not adding
            continue;
        }
        else
        {
            group->my_rays[group->my_size] = new_ray;
            group->my_size++; 
        }
    }

    return;
}

__device__ bool random_interact(Pixel target_pixel, double distance, double* densities, int N)
{
    int i = target_pixel.first, j = target_pixel.second;
    double density = densities[i * N + j];
    double l_ep = density * distance; // effective path length travelled in pixel
    double probability = std::exp(-PARAM_A / l_ep);

    curandState state;
    init_curand_state(&state);
    double rand = curand_uniform_double(&state);
    return (rand < probability);
}

__device__ void
transfer_energy(Ray* ray, Pixel target_pixel, double unscaled_energy, double* densities, double* doses, int N)
{
    int i = target_pixel.first, j = target_pixel.second;
    double density = densities[i * N + j];
    double energy_to_transfer = unscaled_energy * density; // scale energy by pixel density
    double current_ray_energy = ray->get_current_energy();

    // Ray cannot transfer more energy that it has
    energy_to_transfer = fmin(energy_to_transfer, current_ray_energy);

    // Remove energy from ray and add it to pixel dose
    ray->set_current_energy(current_ray_energy - energy_to_transfer);

    doses[i * N + j] += energy_to_transfer;

    return;
}

__device__ int evolve_rays(RayGroup* group, double* densities, double* doses, int N)
{
    int rays_evolved = 0;

    for (int i = 0; i < group->my_size; i++)
    {
        Ray* r = &group->my_rays[i];
        // Only evolve active rays
        if (r->is_active())
        {
            // Trace ray
            TraceHistory rtrace = r->trace();
            Pixel visited_pixel = rtrace.visited;
            double travel_distance = rtrace.distance; // distance traveled in visited pixel
            rays_evolved++;

            if (r->is_primary()) // primary ray
            {
                if (random_interact(visited_pixel, travel_distance, densities, N))
                {
                    // Deposit energy to pixel
                    double energy_to_deposit = PARAM_F * travel_distance * r->get_current_energy();
                    transfer_energy(r, visited_pixel, energy_to_deposit, densities, doses, N);

                    // Spawn secondary rays, transferring remaining energy to them
                    spawn_secondary_rays(group, visited_pixel, r->get_current_energy(), N);
                    r->set_current_energy(0);
                }
            }
            else // secondary ray
            {

                double energy_to_deposit = PARAM_G * travel_distance;

                transfer_energy(r, visited_pixel, energy_to_deposit, densities, doses, N);
            }

            // Deactivate ray if out of energy or outside of the grid bounds
            if (r->get_current_energy() < PARAM_EPSILON || out_of_bounds(r->get_current_pixel(), N))
            {
                r->deactivate();
            }
        }
    }
    return rays_evolved;
}

__device__ void evolve_to_completion(RayGroup* group, double* densities, double* doses, int N)
{
    int rays_evolved = group->my_size;

    while (rays_evolved > 0)
    {
        rays_evolved = evolve_rays(group, densities, doses, N);
    }

    return;
}
__device__ void run_serial(int num_primary_rays, double* densities, double* doses, int N)
{
    // Each primary ray is done serially as its own individual ray group
    for (int i = 0; i < num_primary_rays; i++)
    {
        // printf("Hello from block %d, thread %d. I'm running primary ray %d\n", blockIdx.x, threadIdx.x, i);
        RayGroup primary_ray_group;
        int max_num_rays = 1 + PARAM_KS; // at most one primary ray plus all secondaries
        Ray* rays;
        rays = (Ray*)malloc(max_num_rays * sizeof(Ray));
        primary_ray_group.my_rays = rays;
        primary_ray_group.my_size = 0;

        spawn_primary_ray(&primary_ray_group, N);

        evolve_to_completion(&primary_ray_group, densities, doses, N);

        free(primary_ray_group.my_rays);
    }
    return;
}

/// Kernel function for base-GPU
__global__ void run_rays(int total_num_primary_rays, double* densities, double* doses, int N)
{
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_index < total_num_primary_rays)
    {
        run_serial(1, densities, doses, N);
    }
}

