#include "simulation.cuh"
#include <stdio.h>

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
    //printf("myID: %d\n", tId);
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

    printf("block %d, thread %d spawning primary ray with angle %0.5f\n", blockIdx.x, threadIdx.x, source_angle);

    // Calculate initial ray position
    double horiz_dist_from_center = PARAM_D * N * tan(source_angle); // horizontal distance from center of top edge
    int middle_pixel = N / 2;
    double horiz_dist_from_left = middle_pixel + horiz_dist_from_center;

    // If ray does not miss grid entirely, spawn it
    if (horiz_dist_from_left < 0 || horiz_dist_from_left >= N ||
        (source_angle >= M_PI / 2 && source_angle <= 3 * M_PI / 2))
    {
        // DEBUG(DB_INIT_PRI, std::cout << "New primary ray missed the grid, not adding" << std::endl);
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
        // DEBUG(DB_INIT_PRI, std::cout << "New primary ray added at pixel " << pixel.first << "," << pixel.second
        //                                     << " with angle " << source_angle << std::endl);
    }
    return;
}

__device__ void spawn_secondary_rays(RayGroup* group, Pixel spawn_pixel, double total_energy, int N)
{
    // DEBUG(DB_INIT_SEC, std::cout << "Spawning " << PARAM_KS << " secondary rays from pixel " << spawn_pixel.first <<
    // ","
    //                                 << spawn_pixel.second << "..." << std::endl);
    for (int i = 0; i < PARAM_KS; i++)
    {
        double source_angle = random_source_angle(false); // uniform random source angle
        double partial_energy = total_energy / PARAM_KS;
        Ray new_ray = Ray::secondary_from_center(source_angle, spawn_pixel, partial_energy);
        Pixel current_pixel = new_ray.get_current_pixel();
        if (out_of_bounds(current_pixel, N))
        {
            // DEBUG(DB_INIT_SEC, std::cout << "Ray is out of bounds, not adding" << std::endl);
        }
        else
        {
            group->my_rays[group->my_size] = new_ray;
            group->my_size++;

            // DEBUG(DB_INIT_SEC, std::cout << "Ray is in bounds, adding" << std::endl);
        }
    }
    // DEBUG(DB_INIT_SEC, std::cout << "Done" << std::endl << std::endl);
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
    //printf("Attempting to transfer %0.6f energy\n", energy_to_transfer);
    //printf("Pixel energy before transfer: %0.6f\n", doses[i * N + j]);
    doses[i * N + j] += energy_to_transfer;
    //printf("Pixel energy after transfer: %0.6f\n", doses[i * N + j]);
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
            // DEBUG(DB_TRACE, std::cout << "Tracing ray " << i << std::endl);
            TraceHistory rtrace = r->trace();
            Pixel visited_pixel = rtrace.visited;
            double travel_distance = rtrace.distance; // distance traveled in visited pixel
            rays_evolved++;

            if (r->is_primary()) // primary ray
            {
                // DEBUG(DBevolve_PRI, std::cout << "Primary ray " << i << "  visited pixel " << visited_pixel.first
                //                                               << "," << visited_pixel.second << " with travel dist "
                //                                               << travel_distance
                //                                               << std::endl);
                if (random_interact(visited_pixel, travel_distance, densities, N))
                {
                    printf("block %d, thread %d primary ray interacted at pixel %d,%d\n", blockIdx.x, threadIdx.x,
                           visited_pixel.first, visited_pixel.second);
                    // DEBUG(DBevolve_PRI, std::cout << "Primary ray " << i << " interacted" << std::endl);
                    // Deposit energy to pixel
                    // DEBUG(DBevolve_PRI, std::cout << "Depositing energy to pixel" << std::endl);
                    // DEBUG(DBevolve_PRI, std::cout << "Starting energy " << r->get_current_energy() << std::endl);
                    printf("distance traveled %0.6f\n", travel_distance);
                    double energy_to_deposit = PARAM_F * travel_distance * r->get_current_energy();
                    transfer_energy(r, visited_pixel, energy_to_deposit, densities, doses, N);
                    // DEBUG(DBevolve_PRI, std::cout << "Energy after deposit " << r->get_current_energy() <<
                    // std::endl); DEBUG(DBevolve_PRI, std::cout << "Will spawn secondary rays next" << std::endl <<
                    // std::endl); Spawn secondary rays, transferring remaining energy to them
                    spawn_secondary_rays(group, visited_pixel, r->get_current_energy(), N);
                    r->set_current_energy(0);
                }
                else
                {
                    // DEBUG(DBevolve_PRI, std::cout << "No interaction" << std::endl << std::endl);
                }
            }
            else // secondary ray
            {
                // DEBUG(DBevolve_SEC, std::cout << "Secondary ray " << i << " visited pixel " << visited_pixel.first
                //                                               << "," << visited_pixel.second << " with travel dist "
                //                                               << travel_distance
                //                                               << std::endl);
                double energy_to_deposit = PARAM_G * travel_distance;
                // DEBUG(DBevolve_SEC, std::cout << "Depositing energy to pixel" << std::endl);
                // DEBUG(DBevolve_SEC, std::cout << "Starting energy " << r->get_current_energy() << std::endl);
                // DEBUG(DBevolve_SEC, std::cout << "Unscaled energy to deposit " << energy_to_deposit << std::endl);
                transfer_energy(r, visited_pixel, energy_to_deposit, densities, doses, N);
                // DEBUG(DBevolve_SEC, std::cout << "Energy after deposit " << r->get_current_energy() << std::endl
                //                                               << std::endl);
            }

            // Deactivate ray if out of energy or outside of the grid bounds
            if (r->get_current_energy() < PARAM_MINERGY || out_of_bounds(r->get_current_pixel(), N))
            {
                // DEBUG(DBevolve_SEC, std::cout << "Ray " << i << " is out of energy or bounds, deactivating"
                //                                               << std::endl
                //                                               << std::endl);
                r->deactivate();
            }
        }
    }
    return rays_evolved;
}

__device__ void evolve_to_completion(RayGroup* group, double* densities, double* doses, int N)
{
    int rays_evolved = group->my_size;
    //  int raysevolved = ray_group->size();
    while (rays_evolved > 0)
    {
        rays_evolved = evolve_rays(group, densities, doses, N);
        // DEBUG(DB_GENERAL, std::cout << raysevolved << " rays evolved" << std::endl);
        // DEBUG(DB_GENERAL, std::cout << (ray_group->size() - prev_num_rays) << " rays added" << std::endl <<
        // std::endl)
    }

    return;
}
__device__ void run_serial(int num_primary_rays, double* densities, double* doses, int N)
{
    // Each primary ray is done serially as its own individual ray group
    for (int i = 0; i < num_primary_rays; i++)
    {
        printf("Hello from block %d, thread %d. I'm running primary ray %d\n", blockIdx.x, threadIdx.x, i);
        // Just running one primary ray (and its secondaries) at a time for now
        // TODO: might want to change how memory allocation/reallocation is done
        RayGroup primary_ray_group;
        int max_num_rays = 1 + PARAM_KS; // at most one primary ray plus all secondaries
        Ray* rays;
        rays = (Ray*)malloc(max_num_rays * sizeof(Ray));
        primary_ray_group.my_rays = rays;
        primary_ray_group.my_size = 0;

        // DEBUG(DB_INIT_PRI, std::cout << "Spawning " << num_primary_rays << " primary rays..." << std::endl);
        spawn_primary_ray(&primary_ray_group, N);
        // DEBUG(DB_INIT_PRI, std::cout << "Done, " << ray_group.size() << " rays added" << std::endl << std::endl);
        // DEBUG(DB_GENERAL, std::cout << "Evolving rays..." << std::endl)
        evolve_to_completion(&primary_ray_group, densities, doses, N);
        // DEBUG(DB_GENERAL, std::cout << "Done" << std::endl);
        free(primary_ray_group.my_rays);
    }
    return;
}

/// Kernel function for base-GPU
__global__ void run_rays(int num_primary_rays, double* densities, double* doses, int N)
{
    run_serial(num_primary_rays, densities, doses, N);
}

// TODO could actually make N a command line argument, right?
int main(void)
{
    DEBUG(DB_GPU, std::cout << "Starting simulation, allocating grids" << std::endl);
    int N = 10000; // grid size in pixels per side

    // Storing the N by N grid data as 1D arrays of length N*N
    // such that element i,j is at index i * N + j
    // Currently for the model I think i corresponds to x and j corresponds to y
    double *densities, *doses;
    cudaMallocManaged(&densities, N * N * sizeof(double));
    cudaMallocManaged(&doses, N * N * sizeof(double));

    DEBUG(DB_GPU, std::cout << "Initializing densities" << std::endl);
    // initialize_densities_random(densities, N);
    initialize_densities_constant(densities, N, 0.3);

    DEBUG(DB_GPU, std::cout << "Initializing doses to zero" << std::endl);
    initialize_doses(doses, N);
    size_t heap_limit = pow(2, 26); // default 8MB, this sets to 32 MB
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, heap_limit);
 
    DEBUG(DB_GPU, std::cout << "Writing densities" << std::endl);
    //write_to_csv_file(densities, N, "../../plot/densities.csv");
    write_to_csv_file(densities, N, "densities.csv");
    
    int grid_size = 256; // number of thread blocks
    int block_size = 128;   // TODO 1 thread per block, does this make sense?

    DEBUG(DB_GPU, std::cout << "Running rays on threads" << std::endl);
    int primary_rays_per_thread = 32; // each thread does this many rays in serial
    run_rays<<<grid_size, block_size>>>(primary_rays_per_thread, densities, doses, N);

    // Wait for GPU computation to finish
    cudaDeviceSynchronize();

    DEBUG(DB_GPU, std::cout << "Writing doses" << std::endl);
    //write_to_csv_file(doses, N, "../../plot/doses.csv");
    write_to_csv_file(doses, N, "doses.csv");
    DEBUG(DB_GPU, std::cout << "I'm the main function, look at me!" << std::endl);

    cudaFree(densities);
    cudaFree(doses);

    return 0;
}

