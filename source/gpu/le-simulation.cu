#include "le-simulation.cuh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

// Needed to move this to .cu
std::default_random_engine random_engine{(uint_fast32_t)time(0)}; // seeded random number generator

std::uniform_real_distribution<double> uniform_dist{0.0, 1.0}; // uniform distribution, pass {lowerbound, upperbound}

std::normal_distribution<double> normal_dist{PARAM_MEAN, PARAM_SIGMA}; // normal distribribution, pass {mean, stddev}

////////// GRID INITIALIZATION //////////

__host__ void initialize_doses(double* doses, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            doses[i + j * N] = 0;
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
            densities[i + j * N] = uniform_dist(random_engine);
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
            densities[i + j * N] = density;
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
            densities[i + j * N] = max_density * exp(-(x * x + y * y) / (2 * std_dev * std_dev));
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
                densities[i + j * N] += max_density * exp(-(x * x + y * y) / (2 * std_dev * std_dev));
                highest = fmax(highest, densities[i * N + j]);
            }
        }
    }

    // Normalize the resulting density distribution
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            densities[i + j * N] = max_density * densities[i * N + j] / highest;
        }
    }
    return;
}

////////// OUTPUT //////////

__host__ void write_to_csv_file(double* grid_data, int N, const std::string& filename)
{
    std::ofstream output;
    output.open(filename);
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N - 1; i++)
        {
            output << grid_data[i + j * N] << ",";
        }
        output << grid_data[N - 1 + j * N] << "\n";
    }
    output.close();
    return;
}

////////// RANDOMIZATION //////////
__device__ void init_curand_state(curandState_t* state)
{
    // Initialize random kernel
    int tId = threadIdx.x + (blockIdx.x * blockDim.x);
    curand_init((unsigned long long)clock(), tId, 0, state);
    return;
}

__device__ double uniform_angle_dist(curandState_t* state) { return 2 * M_PI * curand_uniform_double(state); }

__device__ double normal_angle_dist(curandState_t* state, double mean, double std_dev)
{
    return mean + std_dev * curand_normal_double(state);
}

__device__ double random_source_angle_uniform()
{
    curandState state;
    init_curand_state(&state);

    double angle = uniform_angle_dist(&state);

    // Angle of 2 pi goes to 0
    if (abs(angle - 2 * M_PI) < PARAM_EPSILON)
    {
        angle = 0.0;
    }
    return angle;
}

__host__ double random_source_angle_normal()
{
    double angle = normal_dist(random_engine);

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

////////// RAY LOCATION CHECKING //////////

__host__ __device__ Region get_region(Pixel position, int N, int M)
{
    int px = position.first;
    int py = position.second;
    Region region;
    region.first = px / M;
    region.second = py / M;
    return region;
}

__host__ __device__ int get_region_index(Pixel position, int N, int M)
{
    Region region = get_region(position, N, M);
    int L = N / M; // number of regions per side
    return region.first + region.second * L;
}

__device__ bool out_of_bounds(Pixel current_pixel, int N)
{
    return (current_pixel.first < 0 || current_pixel.first >= N || current_pixel.second < 0 ||
            current_pixel.second >= N);
}

////////// RAY CREATION //////////

__host__ void spawn_primary_rays(
    std::vector<RegionGroup>& region_groups, int num_primary_rays, int max_rays_per_ray_group, int N, int M)
{
    int L = N / M; // number of regions per side
    for (int i = 0; i < num_primary_rays; i++)
    {
        // Randomly select source angle from normal distribution
        double source_angle = random_source_angle_normal();

        // Calculate initial ray position
        double horiz_dist_from_center = PARAM_D * N * tan(source_angle); // horizontal distance from center of top edge
        int middle_pixel = N / 2;
        double horiz_dist_from_left = middle_pixel + horiz_dist_from_center;

        // Check if ray missed the grid entirely
        if (horiz_dist_from_left < 0 || horiz_dist_from_left >= N ||
            (source_angle >= M_PI / 2 && source_angle <= 3 * M_PI / 2))
        {
            continue;
        }

        // If not, spawn it
        double horiz_dist_from_left_rounded = floor(horiz_dist_from_left);
        Pixel spawn_pixel;
        spawn_pixel.first = horiz_dist_from_left_rounded;
        spawn_pixel.second = 0; // always starts from top of grid
        double edge_dist = horiz_dist_from_left - horiz_dist_from_left_rounded;
        Ray r = Ray::primary(source_angle, spawn_pixel, PIXEL_EDGE::TOP, edge_dist);

        // Create new ray group for primary ray and add it
        RayGroup primary_ray_group;
        primary_ray_group.my_rays = (Ray*)malloc(max_rays_per_ray_group * sizeof(Ray));
        primary_ray_group.max_size = max_rays_per_ray_group;
        primary_ray_group.my_rays[0] = r; // add the new ray
        primary_ray_group.my_size = 1;

        // Add the new ray group to the appropriate region
        Region region = get_region(spawn_pixel, N, M);
        int region_index = region.first + region.second * L; // index of region within vector of region groups
        region_groups[region_index].push_back(primary_ray_group);
    }
    return;
}

__device__ void spawn_secondary_rays(RayGroup* group, Pixel spawn_pixel, double total_energy, int N)
{
    // int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = 0; i < PARAM_KS; i++)
    {
        double source_angle = random_source_angle_uniform(); // uniform random source angle
        double partial_energy = total_energy / PARAM_KS;
        Ray new_ray = Ray::secondary_from_center(source_angle, spawn_pixel, partial_energy);
        Pixel current_pixel = new_ray.get_current_pixel();
        if (out_of_bounds(current_pixel, N))
        {
            continue;
        }
        group->my_rays[group->my_size] = new_ray;
        group->my_size++;
    }
    return;
}

////////// THREAD GROUP EVOLUTION //////////

__device__ bool random_interact(Pixel target_pixel, double distance, double* densities, int N)
{
    int i = target_pixel.first, j = target_pixel.second;
    double density = densities[i + j * N];
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
    double density = densities[i + j * N];
    double energy_to_transfer = unscaled_energy * density; // scale energy by pixel density
    double current_ray_energy = ray->get_current_energy();

    // Ray cannot transfer more energy that it has
    energy_to_transfer = fmin(energy_to_transfer, current_ray_energy);

    // Remove energy from ray and add it to pixel dose
    ray->set_current_energy(current_ray_energy - energy_to_transfer);

    doses[i + j * N] += energy_to_transfer;

    return;
}

/*
__device__ void print_buffer(RegroupBuffer*& g_buffer_cuda, int num_ray_groups)
{
    int section_size = g_buffer_cuda->section_size;
    ////printf("section size on buffer: %d\n", section_size);
    for (int i = 0; i < num_ray_groups; i++)
    {
        for (int j = 0; j < section_size; j++)
        {
            Pixel pixel = g_buffer_cuda->rays[i + j].get_current_pixel();

            // printf("%d, %d\t", pixel.first, pixel.second);
        }
        // printf("\n");
    }
    // printf("\n");
    for (int i = 0; i < num_ray_groups; i++)
    {
        for (int j = 0; j < section_size; j++)
        {
            int region_index = g_buffer_cuda->region_indices[i + j];

            // printf("%d\t", region_index);
        }
        // printf("\n");
    }
    // printf("\n");
    for (int i = 0; i < num_ray_groups; i++)
    {
        int ray_count = g_buffer_cuda->ray_counts[i];
        // printf("%d\t", ray_count);
    }
    // printf("\n***************************\n");
}
*/

__device__ int evolve_rays(RayGroup* group,
                           int region_index,
                           double* densities,
                           double* doses,
                           int N,
                           int M,
                           RegroupBuffer* g_buffer_cuda,
                           int num_ray_groups)
{
    int rays_evolved = 0;
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
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
                    double energy_to_deposit = PARAM_F * travel_distance * r->get_current_energy();
                    transfer_energy(r, visited_pixel, energy_to_deposit, densities, doses, N);
                    spawn_secondary_rays(group, visited_pixel, r->get_current_energy(), N);
                    r->set_current_energy(0);
                }
            }
            else // secondary ray
            {
                double energy_to_deposit = PARAM_G * travel_distance;
                transfer_energy(r, visited_pixel, energy_to_deposit, densities, doses, N);
            }
            // Check if the ray is still in the region
            int new_region_index = get_region_index(r->get_current_pixel(), N, M);
            // Deactivate ray if out of energy or outside of the grid bounds
            if (r->get_current_energy() < PARAM_EPSILON || out_of_bounds(r->get_current_pixel(), N))
            {
                r->deactivate();
            }
            else if (new_region_index != region_index)
            {
                // printf("thread index: %d, Ray is out of region %d, deactivating and adding to buffer\n",
                r->deactivate();

                int buffer_index = thread_index * g_buffer_cuda->section_size +
                                   g_buffer_cuda->ray_counts[thread_index]; // this thread's next index in buffer

                g_buffer_cuda->rays[buffer_index] = *r; // add ray to buffer
                g_buffer_cuda->region_indices[buffer_index] =
                    new_region_index;                      // add destination region index to buffer
                g_buffer_cuda->ray_counts[thread_index]++; // update buffer size
            }
        }
    }
    return rays_evolved;
}

__device__ void evolve_to_completion(RayGroup* group,
                                     int region_index,
                                     double* densities,
                                     double* doses,
                                     int N,
                                     int M,
                                     RegroupBuffer* g_buffer_cuda,
                                     int num_ray_groups)
{
    // int thread_index = blockIdx.x * blockDim.x + threadIdx.x;

    int rays_evolved = group->my_size;
    while (rays_evolved > 0)
    {
        rays_evolved = evolve_rays(group, region_index, densities, doses, N, M, g_buffer_cuda, num_ray_groups);
        // printf("thread %d evolved %d rays\n", thread_index, rays_evolved);
    }
    // printf("thread %d finished evolving\n", thread_index);
    return;
}

__global__ void run_rays(RayGroup* region_group_arr,
                         int region_group_arr_size,
                         int region_index,
                         double* densities,
                         double* doses,
                         int N,
                         int M,
                         RegroupBuffer* g_buffer_cuda)
{
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_index < region_group_arr_size)
    {
        RayGroup* cur_ray_group = &region_group_arr[thread_index];
        // printf("thread index %d: evolving %d rays\n", thread_index, cur_ray_group->my_size);
        evolve_to_completion(cur_ray_group, region_index, densities, doses, N, M, g_buffer_cuda, region_group_arr_size);
    }

    __syncthreads();
    return;
}
////////// REGION GROUP RUNNING AND PROCESSING //////////

// add rays from g_buffer to new regions
//
__host__ void regroup(std::vector<RegionGroup>& region_groups, RegroupBuffer* g_buffer, int num_ray_groups)
{
    int max_num_rays = g_buffer->section_size;
    DEBUG(DB_HOST, std::cout << "Starting regrouping. Num_ray_groups is " << num_ray_groups << std::endl);
    for (int i = 0; i < num_ray_groups; i++)
    {

        int buffer_index = i * g_buffer->section_size; // index into single ray groups data
        Ray* rays = &g_buffer->rays[buffer_index];     // array of rays to be regrouped
        int* region_indices =
            &g_buffer->region_indices[buffer_index]; // array of regions corresponding to array of rays
        int num_rays = g_buffer->ray_counts[i];      // number of rays to be regrouped
        DEBUG(DB_HOST, std::cout << "Regrouping " << num_rays << " rays." << std::endl);

        for (int j = 0; j < num_rays; j++)
        {

            Ray cur_ray = rays[j]; // current ray to be regrouped
            cur_ray.reactivate();
            int new_region = region_indices[j];                       // region ray is entering
            int region_group_size = region_groups[new_region].size(); // size of the region's last ray group
            if (cur_ray.is_primary() || region_group_size == 0)
            {
                RayGroup new_group;
                new_group.my_rays = (Ray*)malloc(max_num_rays * sizeof(Ray));
                new_group.my_rays[0] = cur_ray;
                new_group.my_size = 1;
                new_group.max_size = max_num_rays;
                region_groups[new_region].push_back(new_group);
            }
            else
            {

                RayGroup last_ray_group =
                    region_groups[new_region][region_group_size - 1]; // last ray group in the region group

                DEBUG(DB_HOST, std::cout << "regrouping a ray into region " << new_region << std::endl);

                int ray_group_size = last_ray_group.my_size;

                if ((ray_group_size == 1 && last_ray_group.my_rays[0].is_primary()) ||
                    ray_group_size == max_num_rays) // last ray group contains single primary ray
                {
                    RayGroup new_group;
                    new_group.my_rays = (Ray*)malloc(max_num_rays * sizeof(Ray));
                    new_group.my_rays[0] = cur_ray;
                    new_group.my_size = 1;
                    new_group.max_size = max_num_rays;
                    region_groups[new_region].push_back(new_group);
                }
                else if (ray_group_size < max_num_rays) // if last ray group is not full, add ray to ray group
                {
                    last_ray_group.my_rays[ray_group_size] = cur_ray;
                    last_ray_group.my_size++;
                    region_groups[new_region][region_group_size - 1] = last_ray_group;
                }
            }
        }
    }

    // Free memory of RegroupBuffer data members on host
    free(g_buffer->rays);
    free(g_buffer->region_indices);
    free(g_buffer->ray_counts);

    return;
}

// allocate a regroup buffer on device
__host__ void init_regroup_buffer_cuda(RegroupBuffer*& g_buffer_cuda, int max_num_rays, int num_ray_groups)
{
    DEBUG(DB_HOST, std::cout << "Initializing regroup buffer on device" << std::endl);
    // Need to allocate device memory for the members of struct BEFORE copying struct
    // Host can't modify device memory
    RegroupBuffer temp; // temporary buffer on host
    temp.section_size = max_num_rays;

    Ray* rays_cuda;
    cudaMalloc(&rays_cuda, num_ray_groups * max_num_rays * sizeof(Ray));
    cudaMemset(rays_cuda, 0, num_ray_groups * max_num_rays * sizeof(Ray));
    temp.rays = rays_cuda;

    int* region_indices_cuda;
    cudaMalloc(&region_indices_cuda, num_ray_groups * max_num_rays * sizeof(int));
    cudaMemset(region_indices_cuda, 0, num_ray_groups * max_num_rays * sizeof(int));
    temp.region_indices = region_indices_cuda;

    int* ray_counts_cuda;
    cudaMalloc(&ray_counts_cuda, num_ray_groups * sizeof(int));
    cudaMemset(ray_counts_cuda, 0, num_ray_groups * sizeof(int)); // initialize all ray counts to zero
    temp.ray_counts = ray_counts_cuda;

    // Now allocate device memory for buffer and copy from host to device
    cudaMalloc(&g_buffer_cuda, sizeof(RegroupBuffer));
    error_check(cudaMemcpy(g_buffer_cuda, &temp, sizeof(RegroupBuffer), cudaMemcpyHostToDevice));

    DEBUG(DB_HOST, std::cout << "Done" << std::endl);
    return;
}

// allocate a regroup buffer on host and copy the contents of device's regroup buffer to it
__host__ void
copy_regroup_buffer_host(RegroupBuffer*& g_buffer, RegroupBuffer*& g_buffer_cuda, int max_num_rays, int num_ray_groups)
{
    DEBUG(DB_HOST, std::cout << "Allocating memory for buffer on host" << std::endl);
    g_buffer = (RegroupBuffer*)malloc(sizeof(RegroupBuffer));
    g_buffer->rays = (Ray*)malloc(num_ray_groups * max_num_rays * sizeof(Ray));
    g_buffer->region_indices = (int*)malloc(num_ray_groups * max_num_rays * sizeof(int));
    g_buffer->ray_counts = (int*)malloc(num_ray_groups * sizeof(int));

    RegroupBuffer* temp_buffer = (RegroupBuffer*)malloc(sizeof(RegroupBuffer));
    error_check(cudaMemcpy(temp_buffer, g_buffer_cuda, sizeof(RegroupBuffer), cudaMemcpyDeviceToHost));

    g_buffer->section_size = temp_buffer->section_size;

    DEBUG(DB_HOST, std::cout << "Copying buffer from device to host" << std::endl);
    error_check(cudaMemcpy(g_buffer->rays, temp_buffer->rays, num_ray_groups * max_num_rays * sizeof(Ray),
                           cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(g_buffer->region_indices, temp_buffer->region_indices,
                           num_ray_groups * max_num_rays * sizeof(int), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(g_buffer->ray_counts, temp_buffer->ray_counts, num_ray_groups * sizeof(int),
                           cudaMemcpyDeviceToHost));

    // Free memory of RegroupBuffer data members on device, and temp buffer on host
    cudaFree(temp_buffer->rays);
    cudaFree(temp_buffer->region_indices);
    cudaFree(temp_buffer->ray_counts);
    free(temp_buffer);

    DEBUG(DB_HOST, std::cout << "Done" << std::endl);
    return;
}

__host__ std::vector<int> get_forward_schedule(int L)
{
    std::vector<int> schedule;
    for (int a = 0; a < (2 * L - 1); a++)
    {
        int i = min(a, L - 1);
        int j = max(0, a + 1 - L);
        while (i >= 0 && j < L)
        {
            int task_index = i + j * L;
            i--;
            j++;
            schedule.push_back(task_index);
        }
    }
    return schedule;
}

__host__ void run_region_groups(std::vector<RegionGroup>& region_groups, double* densities, double* doses, int N, int M)
{
    // Get the schedule
    int L = N / M;                                            // number of regions per side
    std::vector<int> forward_sched = get_forward_schedule(L); // linear indices of regions in diagonal order
    std::vector<int> full_sched = forward_sched;
    full_sched.insert(full_sched.end(), forward_sched.rbegin(), forward_sched.rend()); // concat forward and reverse

    // Go through schedule in passes until complete
    bool active_rays = true;
    while (active_rays) // keep going until no rays are left
    {
        DEBUG(DB_HOST, std::cout << "Pass" << std::endl);
        active_rays = false;
        for (std::vector<int>::iterator it = full_sched.begin(); it != full_sched.end(); it++)
        {

            int region_index = *it;
            RegionGroup& cur_region_group = region_groups[region_index]; // current region group REFERENCE
            DEBUG(DB_HOST, std::cout << "Running region group " << region_index << std::endl);
            DEBUG(DB_HOST, std::cout << "It has " << region_groups[region_index].size() << " ray groups" << std::endl);

            // Only do things if the region has ray groups
            int cur_region_group_size = cur_region_group.size();
            if (cur_region_group_size > 0)
            {
                active_rays = true;
                RegroupBuffer* g_buffer; // empty host regroup buffer. Will be filled in by run_region_group

                run_region_group(cur_region_group, region_index, densities, doses, N, M, g_buffer); // run group

                regroup(region_groups, g_buffer, cur_region_group_size); // regroup

                free(g_buffer); // free host buffer
            }
        }
    }
    return;
}

__host__ void run_region_group(RegionGroup& region_group,
                               int region_index,
                               double* densities,
                               double* doses,
                               int N,
                               int M,
                               RegroupBuffer*& g_buffer)
{
    // Set device memory limits
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, GPU_HEAP_LIMIT);

    // First copy rays to ray groups on device, done by just replacing host pointers with device pointers
    int num_ray_groups = region_group.size(); // number of ray groups in current region group
    DEBUG(DB_HOST, std::cout << "Copying rays from host to device" << std::endl);
    int max_num_rays = region_group[0].max_size; // all ray groups have same max size so just get any max size
    for (int g = 0; g < num_ray_groups; g++)
    {

        Ray* rays_cuda;
        cudaMalloc(&rays_cuda, max_num_rays * sizeof(Ray)); // allocated memory on device
        cudaMemcpy(rays_cuda, region_group[g].my_rays, max_num_rays * sizeof(Ray),
                   cudaMemcpyHostToDevice);               // copy from host to device
        Ray* old_host_rays_ptr = region_group[g].my_rays; // pointer to rays on host
        region_group[g].my_rays = rays_cuda;              // this is now a device pointer NOT a host pointer
        free(old_host_rays_ptr);                          // free host memory
    }

    DEBUG(DB_HOST, std::cout << "Copying ray groups from host to device" << std::endl);
    // Copy region group to GPU (std::vector on host to array on device)
    RayGroup* region_group_cuda_arr;
    cudaMalloc(&region_group_cuda_arr, num_ray_groups * sizeof(RayGroup)); // allocated memory on device
    cudaMemcpy(region_group_cuda_arr, region_group.data(), num_ray_groups * sizeof(RayGroup),
               cudaMemcpyHostToDevice); // copy from host to device

    // Allocate regroup buffer on DEVICE
    RegroupBuffer* g_buffer_cuda;                                          // empty device regroup buff0er
    init_regroup_buffer_cuda(g_buffer_cuda, max_num_rays, num_ray_groups); // allocate device regroup buffer

    // Run thread groups in parallel
    int grid_size = 1 + num_ray_groups / GPU_BLOCK_SIZE;
    int block_size = GPU_BLOCK_SIZE;
    DEBUG(DB_HOST,
          std::cout << "Calling run_rays with grid_size,block_size " << grid_size << "," << block_size << std::endl);
    run_rays<<<grid_size, block_size>>>(region_group_cuda_arr, num_ray_groups, region_index, densities, doses, N, M,
                                        g_buffer_cuda);

    // Wait for GPU computation to finish
    error_check(cudaDeviceSynchronize());
    // printf("Devices are synchronized\n");

    // Copy g_buffer back to host buffer
    copy_regroup_buffer_host(g_buffer, g_buffer_cuda, max_num_rays,
                             num_ray_groups); // copy g_buffer back to host buffer
    // Free device memory
    // First free ray group pointers on device, which are still stored on host
    for (int g = 0; g < num_ray_groups; g++)
    {
        cudaFree(region_group[g].my_rays);
    }
    cudaFree(g_buffer_cuda);
    cudaFree(region_group_cuda_arr);

    // Clear region group vector because we messed with its memory, and its rays are all going to be run
    region_group.clear();

    return;
}

