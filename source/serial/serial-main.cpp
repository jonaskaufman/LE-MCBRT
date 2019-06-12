#include "../parameters.hpp"
#include "simulation.hpp"

#include <chrono>
#include <iostream>
//#include <stdlib.h>

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cerr << "Invalid number of arguments [grid size] [number of primary rays] [density initialization method] "
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

    // Create simulation
    Simulation s = Simulation(N);

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

        s.initialize_densities_constant(density);
        break;
    }
    case 'R': // random
    {
        s.initialize_densities_random();
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

        s.initialize_densities_centered_gaussian(max_density, spread);
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

        s.initialize_densities_random_gaussians(number, max_density, spread);
        break;
    }

    default:
    {
        std::cerr << "Invalid initialization method given" << std::endl;
        exit(1);
        break;
    }
    }

    std::cout << "Serial version" << std::endl << std::endl;

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

    // Write densities
    std::cout << "Writing densities" << std::endl;
    std::cout << std::endl;
    s.write_densities_to_file("densities.csv");

    // Run a given number of rays
    std::cout << "Spawning and running " << num_primary_rays << " primary rays" << std::endl;

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    int batch_size = 100; // run in batches to avoid segfaults
    s.run_serial(num_primary_rays, batch_size);

    // Stop timer
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;

    // Write result
    std::cout << std::endl;
    std::cout << "Writing doses" << std::endl;
    s.write_doses_to_file("doses.csv");

    std::cout << "All done" << std::endl;
    return 0;
}

