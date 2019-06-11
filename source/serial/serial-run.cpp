#include "../parameters.hpp"
#include "simulation.hpp"

#include <chrono>
#include <iostream>
#include <stdlib.h>

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
    const int ray_count = atoi(argv[2]);
    if (ray_count < 1)
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

    DEBUG(DB_GENERAL, std::cout << "Densities initialized" << std::endl << std::endl);

    // Write densities
    std::cout << "Writing densities" << std::endl;
    s.write_densities_to_file("densities.csv");

    auto start = std::chrono::high_resolution_clock::now();

    // Run a given number of rays
    std::cout << "Spawning and running " << ray_count << " primary rays" << std::endl;
    int batch_size = 10000;
    s.run_serial(ray_count, batch_size);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;

    // Write result
    std::cout << "Writing doses" << std::endl;
    s.write_doses_to_file("doses.csv");

    std::cout << "All done" << std::endl;
    return 0;
}

