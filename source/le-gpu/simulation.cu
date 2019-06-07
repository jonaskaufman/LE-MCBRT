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

__host__ void extend_ray_groups(RegionGroup *region_groups, int rays_per_group, int region_index, int max_index){
    RayGroup *ray_groups = region_groups[region_index].my_ray_groups;
    ray_groups = (RayGroup *) realloc(ray_groups, (max_index + 0)* 2 * sizeof(RayGroup));
    region_groups[region_index].my_ray_groups = ray_groups;
    for (int i = max_index; i < max_index * 2; i++){
        region_groups[region_index].my_ray_groups[i].my_rays = (Ray *) malloc(rays_per_group * sizeof(Ray));
        region_groups[region_index].my_ray_groups[i].my_size = 0;
        region_groups[region_index].my_ray_groups[i].max_size = rays_per_group;
    }
    region_groups[region_index].max_size = max_index * 2;
    
}

__host__ void initialize_region_groups(RegionGroup *region_groups, int num_regions, int ray_groups_per_region, int rays_per_group){
    for (int i = 0; i < num_regions; i++){
        region_groups[i].my_ray_groups = (RayGroup *) malloc(ray_groups_per_region * sizeof(RayGroup));
        region_groups[i].my_size = -1;
        region_groups[i].max_size = ray_groups_per_region;
        for (int j = 0; j < ray_groups_per_region; j++){
            Ray *r = (Ray *) malloc(rays_per_group * sizeof(Ray));
            region_groups[i].my_ray_groups[j].my_rays = r;
            region_groups[i].my_ray_groups[j].my_size = 0;
            region_groups[i].my_ray_groups[j].max_size = rays_per_group;
        }
    }
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

__device__ double normal_angle_dist(curandState_t* state, double mean, double std_dev)
{
    return mean + std_dev * curand_normal_double(state);
}

__device__ double random_source_angle_uniform()
{
    curandState state;
    init_curand_state(&state);

    double angle = uniform_angle_dist(&state);
    
    // Normalize angle
    if (abs(angle - 2 * M_PI) < PARAM_EPSILON)// if angle == 2pi, angle = 0
    {
        angle = 0;
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


__device__ bool out_of_bounds(Pixel current_pixel, int N)
{
    return (current_pixel.first < 0 || current_pixel.first >= N || current_pixel.second < 0 ||
            current_pixel.second >= N);
}

__host__ void spawn_primary_rays(RegionGroup *region_groups, int num_primary_rays, int N, int M)
{
    for (int i = 0; i < num_primary_rays; i++){

        // Randomly select source angle from normal distribution
        double source_angle = random_source_angle_normal();
        // Calculate initial ray position
        double horiz_dist_from_center = PARAM_D * N * tan(source_angle); // horizontal distance from center of top edge
        int middle_pixel = N / 2;
        double horiz_dist_from_left = middle_pixel + horiz_dist_from_center;
        
        // If ray does not miss grid entirely, spawn it
        if (horiz_dist_from_left < 0 || horiz_dist_from_left >= N ||
            (source_angle >= M_PI / 2 && source_angle <= 3 * M_PI / 2))
        {
            continue;
            // DEBUG(DB_INIT_PRI, std::cout << "New primary ray missed the grid, not adding" << std::endl);
        }
        
        double horiz_dist_from_left_rounded = floor(horiz_dist_from_left);
        Pixel spawn_pixel;
        spawn_pixel.first = horiz_dist_from_left_rounded;
        spawn_pixel.second = 0; // always starts from top of grid
        double edge_dist = horiz_dist_from_left - horiz_dist_from_left_rounded;

        Region region = get_region(spawn_pixel, N, M);
        Ray r = Ray::primary(source_angle, spawn_pixel, PIXEL_EDGE::TOP, edge_dist, region);
        
        int Rx = floor((float) N / M);
        int region_index = region.second * Rx + region.first;
        int max_index = region_groups[region_index].max_size;
        int group_index = region_groups[region_index].my_size + 1;
        

        if (group_index > max_index - 1){
            //printf("resizing array of ray groups\n");
            //printf("max_index: %d\t group_index: %d\n", max_index, group_index);
            //printf("Group: %d\t Region: %d, %d\n\n", region_index, region.first, region.second);
            extend_ray_groups(region_groups, PARAM_KS + 1, region_index, max_index);
            //printf("finished resizing\n");
        }
        //printf("max_index: %d\t group_index: %d\n", region_groups[region_index].max_size, group_index);
        //printf("Group: %d\t Region: %d, %d\n\n", region_index, region.first, region.second);
        region_groups[region_index].my_ray_groups[group_index].my_rays[0] = r; // add new ray to the beginning of ray group
        region_groups[region_index].my_ray_groups[group_index].my_size = 1; // update size
        region_groups[region_index].my_size = group_index;
        //printf("Angle: %.2f\tPosition: %d, %d\n", source_angle, spawn_pixel.first, spawn_pixel.second);
        //printf("max_index: %d\t group_index: %d\n", max_index, group_index);
        //printf("Group: %d\t Region: %d, %d\n\n", region_index, region.first, region.second);
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
        double source_angle = random_source_angle_uniform(); // uniform random source angle
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

__device__ int evolve_rays(RayGroup* group, double* densities, double* doses, int N, int M)
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
                    ///////////////////////////spawn_secondary_rays(group, visited_pixel, r->get_current_energy(), N); // UNCOMMENT 
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

__device__ void evolve_to_completion(RayGroup* group, double* densities, double* doses, int N, int M)
{
    int rays_evolved = group->my_size;
    //  int raysevolved = ray_group->size();
    while (rays_evolved > 0)
    {
        rays_evolved = evolve_rays(group, densities, doses, N, M);
        // DEBUG(DB_GENERAL, std::cout << raysevolved << " rays evolved" << std::endl);
        // DEBUG(DB_GENERAL, std::cout << (ray_group->size() - prev_num_rays) << " rays added" << std::endl <<
        // std::endl)
    }

    return;
}
__device__ void run_serial(RegionGroup *region_groups, Region cur_region, double* densities, double* doses, int N, int M)
{
    printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
    // Each primary ray is done serially as its own individual ray group
    

    /// Sanity check to see if GPU received all the correct data. It does
    /*
    int num_regions = pow(ceilf(N / M), 2); // number of regions
    printf("N: %d\tM: %d\tnum_regions: %d\n", N, M, num_regions);
    int count = 0;
    for (int i = 0; i < num_regions; i++){
        RegionGroup cur_region_group = region_groups[i];
        int cur_region_group_size = cur_region_group.my_size;
        //printf("cur_region_group_size: %d\n", cur_region_group_size);
        for (int j = 0; j < cur_region_group_size; j++){
            RayGroup cur_ray_group = cur_region_group.my_ray_groups[j];
                Ray r = cur_ray_group.my_rays[0];
                Pixel position = r.get_current_pixel();
                Region region = r.get_current_region();
                printf("Primary %d: G: %d, P: %d, %d\tR: %d, %d\n", count, i, position.first, position.second, region.first, region.second);
                count++;
        }
    }
    return;
    */
    int Rx = floor((float) N / M);
    int group_index = cur_region.second * Rx + cur_region.first;
    RegionGroup cur_region_group = region_groups[group_index];
    int num_ray_groups = cur_region_group.my_size;

    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_index < num_ray_groups){
        RayGroup cur_ray_group = cur_region_group.my_ray_groups[thread_index];
        evolve_to_completion(&cur_ray_group, densities, doses, N, M);
    }
    
    

    //free(primary_ray_group.my_rays);

    return;
    
}

/// Kernel function
__global__ void run_rays(RegionGroup *region_groups, Region *schedule, double* densities, double* doses, int N, int M)
{
    int num_regions = pow(ceilf(N / M), 2);
    Region cur_region;
    for (int i = 0; i < num_regions; i++){
        cur_region = schedule[i];
        run_serial(region_groups, cur_region, densities, doses, N, M);
    }
    
}

// TODO could actually make N a command line argument, right?
int main(void)
{
    DEBUG(DB_GPU, std::cout << "Starting simulation, allocating grids" << std::endl);
    int N = 1000; // grid size in pixels per side
    int M = 250; // region size in pixels per side
    int num_primary_rays = 100; // number of primary rays to be scattered across the grid
    int rays_per_group = PARAM_KS + 1; // initial capacity of a ray group
    int ray_groups_per_region = 10; // intial capacity of a region group

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
    size_t heap_limit = pow(2, 26); // default 8MB, this sets to 64 MB
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, heap_limit);
 
    DEBUG(DB_GPU, std::cout << "Writing densities" << std::endl);
    write_to_csv_file(densities, N, "densities.csv");
    


    DEBUG(DB_GPU, std::cout << "Running rays on threads" << std::endl);
    int num_regions = pow(ceil(N / M), 2); // number of regions
    DEBUG(DB_GPU, std::cout << "Number of regions: " << num_regions << std::endl);

    // create all primary rays on host. One primary ray per (thread) ray group
    RegionGroup *region_groups = (RegionGroup *) malloc(num_regions * sizeof(RegionGroup));
    initialize_region_groups(region_groups, num_regions, ray_groups_per_region, rays_per_group);
    spawn_primary_rays(region_groups, num_primary_rays, N, M);
    
    // print all primary rays spawned
    int count = 0;
    for (int i = 0; i < num_regions; i++){
        RegionGroup cur_region_group = region_groups[i];
        int cur_region_group_size = cur_region_group.my_size;
        for (int j = 0; j < cur_region_group_size; j++){
            RayGroup cur_ray_group = cur_region_group.my_ray_groups[j];
                Ray r = cur_ray_group.my_rays[0];
                Pixel position = r.get_current_pixel();
                Region region = r.get_current_region();
                printf("Primary %d: G: %d, P: %d, %d\tR: %d, %d\n", count, i, position.first, position.second, region.first, region.second);
                count++;
        }
    }

    // copy all host data to GPU
    RegionGroup *region_groups_cuda;
    cudaMallocManaged(&region_groups_cuda, num_regions * sizeof(RegionGroup));
    cudaMemcpy(region_groups_cuda, region_groups, num_regions * sizeof(RegionGroup), cudaMemcpyHostToDevice);
    for (int i = 0; i < num_regions; i++){
        ray_groups_per_region = region_groups[i].my_size;
        RayGroup *ray_groups_cuda;
        cudaMallocManaged(&ray_groups_cuda, ray_groups_per_region * sizeof(RayGroup));

        cudaMemcpy(ray_groups_cuda, region_groups[i].my_ray_groups, 
            ray_groups_per_region * sizeof(RayGroup), cudaMemcpyHostToDevice);

        region_groups_cuda[i].my_ray_groups = ray_groups_cuda;
        region_groups_cuda[i].my_size = region_groups[i].my_size;
        region_groups_cuda[i].max_size = region_groups[i].max_size;

        for (int j = 0; j < ray_groups_per_region; j++){
            Ray *rays_cuda;
            cudaMallocManaged(&rays_cuda, rays_per_group * sizeof(Ray));

            cudaMemcpy(rays_cuda, region_groups[i].my_ray_groups[j].my_rays,
                rays_per_group * sizeof(Ray), cudaMemcpyHostToDevice);
            region_groups_cuda[i].my_ray_groups[j].my_rays = rays_cuda;
            region_groups_cuda[i].my_ray_groups[j].my_size = region_groups[i].my_ray_groups[j].my_size;
            region_groups_cuda[i].my_ray_groups[j].max_size = region_groups[i].my_ray_groups[j].max_size;
        }
    }

    //TODO: Free host data?
    printf("building schedule\n");
    Region *schedule;//
    cudaMallocManaged(&schedule, num_regions * sizeof(Region));
    // = (Region *) malloc(num_regions * sizeof(Region));
    build_schedule(schedule, N, M, true, true);
    printf("Schedule is: ");
    for (int i = 0; i < num_regions; i++){
        printf("%d, %d\n", schedule[i].first, schedule[i].second);
    }

    int grid_size = 4; // number of thread blocks
    int block_size = 4;   // TODO 1 thread per block, does this make sense?
    run_rays<<<grid_size, block_size>>>(region_groups_cuda, schedule, densities, doses, N, M);
    // TODO: Call run_rays again with new schedule (backward pass, etc)
          // Or schedule can be made within run_rays. Had device/host issue earlier but haven't tried recently

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

Region get_region(Pixel position, int N, int M){
    int px = position.first;
    int py = position.second;

    Region region;
    region.first = floor((float) px / M);
    region.second = floor((float) py / M);
    return region;
}

__host__ __device__ void build_schedule(Region *schedule, int N, int M, bool forward, bool lr_diag){
    Region r0 = {0,0};
    Region r1 = {0,1};
    Region r2 = {1,0};
    Region r3 = {0,2};
    Region r4 = {1,1};
    Region r5 = {2,0};
    Region r6 = {0,3};
    Region r7 = {1,2};
    Region r8 = {2,1};
    Region r9 = {3,0};
    Region r10 = {1,3};
    Region r11 = {2,2};
    Region r12 = {3,1};
    Region r13 = {2,3};
    Region r14 = {3,2};
    Region r15 = {3,3};
    schedule[0] = r0;
    schedule[1] = r1;
    schedule[2] = r2;
    schedule[3] = r3;
    schedule[4] = r4;
    schedule[5] = r5;
    schedule[6] = r6;
    schedule[7] = r7;
    schedule[8] = r8;
    schedule[9] = r9;
    schedule[10] = r10;
    schedule[11] = r11;
    schedule[12] = r12;
    schedule[13] = r13;
    schedule[14] = r14;
    schedule[15] = r15;
    /*
    Region cur_region = {0,0};
    int index = 0;
    while (cur_region.first != -1 || cur_region.second != -1){
        printf("cur_region: %d, %d\n", cur_region.first, cur_region.second);
        schedule[index++] = cur_region;
        cur_region = get_next_region(cur_region, N, M, forward, lr_diag);
    }
    */
}

__host__ __device__ Region get_next_region(Region cur_region, int N, int M, bool forward, bool lr_diag){
    int dx, dy; // which direction is "progress" from V to V', not from cur_region to next_region
    int edge_x, edge_y; // boundary condition. If at some boundary, go somewhere else
    int max_region = ceilf(N / M) - 1;
    int rx = cur_region.first, ry = cur_region.second;

    if (forward == true){   // forward pass
        if (lr_diag == true){ // V is top left node. V' is bottom right node. so diagonal between them is going left to right
            dx = 1, dy = 1; // right and down
            edge_x = 0, edge_y = 0;
        }
        else{
            dx = -1, dy = 1; // left and down
            edge_x = max_region, edge_y = 0;
        }
    }
    else{                   // backwards pass
        if (lr_diag == true){
            dx = -1, dy = -1; // left and up
            edge_x = max_region, edge_y = max_region;
        }
        else{
            dx = 1, dy = -1; // right and up
            edge_x = 0, edge_y = max_region;
        }
    }
    // TODO implement backwards and other diagonal. Right now, it's just forward, lr logic
    
    Region result;
    if (rx == edge_x && ry == edge_y){
        result.first = rx; result.second = ry + dy;
        return result;
    }
    else if (ry == edge_y){
        result.first = 0; result.second = rx + dy;
        return result;
    }
    else{
        result.first = rx + dx; result.second = ry - dy;
        if (result.first > max_region || result.second > max_region){
            result.first = max_region, result.second = max_region;
        }
        if(rx == max_region && ry == max_region){
            result.first = -1, result.second = -1;
        }
        return result;

    }
    //return result;


}

