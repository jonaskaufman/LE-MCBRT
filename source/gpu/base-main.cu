#include "../parameters.hpp"
#include "base-simulation.cuh"

#include <chrono>
#include <iostream>
#include <stdlib.h>

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cerr << "Invalid number of arguments [grid size] [number of primary rays] [density "
                     "initialization method] "
                     "[additional]"
                  << std::endl;
        std::cerr << "Initialization methods and their additional arguments are:" << std::endl;
        std::cerr << "C - Constant [density]" << std::endl;
        std::cerr << "R - Random" << std::endl;
        std::cerr << "G - Gaussian [maximum density] [spread]" << std::endl;
        std::cerr << "M - Multiple random Gaussians [number of Gaussians] [maximum density] [spread]" << std::endl;
        exit(1);
    }

    // Grid size
    const int N = atoi(argv[1]);
    if (N < 1)
    {
        std::cerr << "Grid size must be greater than 0" << std::endl;
        exit(1);
    }

    // Number of rays
    const int num_primary_rays = atoi(argv[2]);
    if (num_primary_rays < 1)
    {
        std::cerr << "Ray count must be greater than 0" << std::endl;
        exit(1);
    }

    // Set up grids
    // NOTE: all 2D arrays unwrapped as 1D arrays/vectors use linear indexing
    // of the form i,j -> i + j * edge_dimension
    DEBUG(DB_HOST, std::cout << "Starting simulation, allocating grids" << std::endl);
    double *densities, *doses;
    cudaMallocManaged(&densities, N * N * sizeof(double));
    cudaMallocManaged(&doses, N * N * sizeof(double));

    // Initialize densities and doses
    DEBUG(DB_HOST, std::cout << "Initializing doses to zero" << std::endl);
    initialize_doses(doses, N);

    DEBUG(DB_HOST, std::cout << "Initializing densities" << std::endl);
    // Initialize densities according to given method
    int init_method_arg = 3;
    switch (argv[init_method_arg][0])
    {
    case 'C': // constant
    {
        if (argc < (init_method_arg + 1 + 1))
        {
            std::cerr << "Not enough arguments for chosen initialization option" << std::endl;
            exit(1);
        }

        const double density = strtod(argv[init_method_arg + 1], NULL);
        if (density < 0 || density > 1)
        {
            std::cerr << "Density must be between 0 and 1" << std::endl;
            exit(1);
        }
        initialize_densities_constant(densities, N, density);
        break;
    }
    case 'R': // random
    {
        initialize_densities_random(densities, N);
        break;
    }

    case 'G': // centered Gaussian
    {
        if (argc < (init_method_arg + 1 + 2))
        {
            std::cerr << "Not enough arguments for chosen initialization option" << std::endl;
            exit(1);
        }

        const double max_density = strtod(argv[init_method_arg + 1], NULL);
        if (max_density < 0 || max_density > 1)
        {
            std::cerr << "Density must be between 0 and 1" << std::endl;
            exit(1);
        }
        const double spread = strtod(argv[init_method_arg + 2], NULL);

        initialize_densities_centered_gaussian(densities, N, max_density, spread);
        break;
    }

    case 'M': // multiple random Gaussians
    {
        if (argc < (init_method_arg + 1 + 3))
        {
            std::cerr << "Not enough arguments for chosen initialization option" << std::endl;
            exit(1);
        }

        const int number = atoi(argv[init_method_arg + 1]);
        if (number < 1)
        {
            std::cerr << "Number of Gaussians must be greater than 0" << std::endl;
            exit(1);
        }

        const double max_density = strtod(argv[init_method_arg + 2], NULL);
        if (max_density < 0 || max_density > 1)
        {
            std::cerr << "Density must be between 0 and 1" << std::endl;
            exit(1);
        }
        const double spread = strtod(argv[init_method_arg + 3], NULL);

        initialize_densities_random_gaussians(densities, N, number, max_density, spread);
        break;
    }

    default:
    {
        std::cerr << "Invalid initialization method given" << std::endl;
        exit(1);
        break;
    }
    }

    std::cout << "Base-GPU version" << std::endl << std::endl;

    // Print parameters
    std::cout << "Listing parameters" << std::endl;
    std::cout << "Grid side length N: " << N << std::endl;
    std::cout << "Source offset distance D / N: " << PARAM_D << std::endl;
    std::cout << "Source angle mean, std dev: " << PARAM_MEAN << ", " << PARAM_SIGMA << std::endl;
    std::cout << "Initial primary ray energy E0: " << PARAM_E0 << std::endl;
    std::cout << "Interaction probability scaling a: " << PARAM_A << std::endl;
    std::cout << "Primary ray energy deposit fraction f: " << PARAM_F << std::endl;
    std::cout << "Secondary ray energy deposit constant g: " << PARAM_G << std::endl;
    std::cout << "Number of secondary rays per primary k_s: " << PARAM_KS << std::endl;
    std::cout << std::endl;

    std::cout << "Threads per block: " << GPU_BLOCK_SIZE << std::endl;
    std::cout << std::endl;

    // Write densities
    std::cout << "Writing densities" << std::endl;
    std::cout << std::endl;
    write_to_csv_file(densities, N, "densities.csv");

    // Run the given number of rays on threads (includes spawning)
    std::cout << "Spawning and running " << num_primary_rays << " primary rays" << std::endl;

    cudaDeviceSetLimit(cudaLimitMallocHeapSize, GPU_HEAP_LIMIT); // increase device heap limit
    int grid_size = 1 + num_primary_rays / GPU_BLOCK_SIZE;       // number of thread blocks
    int block_size = GPU_BLOCK_SIZE;                             // threads per block

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    run_rays<<<grid_size, block_size>>>(num_primary_rays, densities, doses, N);
    cudaDeviceSynchronize(); // wait for GPU computation to finish

    // Stop timer
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;

    // Write result
    std::cout << std::endl;
    std::cout << "Writing doses" << std::endl;
    write_to_csv_file(doses, N, "doses.csv");

    // Free memory
    cudaFree(densities);
    cudaFree(doses);

    std::cout << "All done" << std::endl;
    return 0;
}

