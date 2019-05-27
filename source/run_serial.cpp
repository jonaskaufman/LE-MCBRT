#include "simulation.hpp"
#include "parameters.hpp"

#include <iostream>
#include <stdlib.h>

int main(int argc, char **argv) 
{
    if (argc != 4){
        std::cerr << "Invalid number of arguments. [grid size] [density (or -1 for random)] [# of primary rays]" << std::endl;
        exit(1);
    }

    const int N = atoi(argv[1]);
    
    if (N < 1) {
        std::cerr << "grid size must be greater than 0" << std::endl;
        exit(1);
    }
    DEBUG(DB_ARGPASS, std::cout << "N is " << N << std::endl);

    const double density = strtod(argv[2], NULL);
    if ( (density != -1) && (density < 0 || density > 1) ){
        std::cerr << "density must be between 0 and 1 (or -1 for random densities)";
    }
    DEBUG(DB_ARGPASS, std::cout << "density is " << density << std::endl);
    
    const int ray_count = atoi(argv[3]);
    if (ray_count < 1){
        std::cerr << "ray count must be greater than 0";
    }
    DEBUG(DB_ARGPASS, std::cout << "ray count is " << ray_count << std::endl);

    DEBUG(DB_ARGPASS, std::cout << "about to create simulation" << std::endl);

    // create simulation class
    SimulationSerial s = SimulationSerial(N, density, ray_count);
    DEBUG(DB_SIMCONST, std::cout << "created simulation" << std::endl);
    //DEBUG(DB_SIMCONST, s.print_m_densities());
    //DEBUG(DB_SIMCONST, std::cout << "******************" << std::endl);
    //DEBUG(DB_SIMCONST, s.print_m_doses());


    // initialize grid
    
    // run a given number of rays
    
    // write result
    s.write_to_file();
    

    return 0;
}

