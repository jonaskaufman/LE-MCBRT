

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "parameters.hpp"
#endif

#ifndef RAY_H
#define RAY_H
//#include "ray.hpp"
#endif

#include <vector>
#include <random>

#ifndef _IOSTREAM_H
#define _IOSTREAM
#include <iostream>
#endif



class SimulationSerial
{
    public: 
    // SimulationSerial() = delete;
    
    SimulationSerial(const int N, const double density); /// will this work? for static arrays?

    

    /// Label the pixel edges
    enum class PIXEL_EDGE {TOP, BOTTOM, LEFT, RIGHT};

    /// Initalize densities
    void initialize_densities_random();
    void initialize_densities_constant(const double density);

    /// Run simulation for a given number of primary rays
    void run(int num_primary_rays);

    /// Get array of doses
    double* get_doses();

    // Print data
    void print_m_densities();
    void print_m_doses();
    
    private:
    const int N;
    double** m_densities;         /// pixel density values
    double** m_doses;                   /// pixel dose values
    //std::vector<Ray::Ray> m_rays;                /// active rays

    /// Randomly sample source angle for primary rays
    double _random_source_angle();

    /// Generate new primary ray
    void _spawn_primary_ray();

    /// Evolve all rays by one step
    void _evolve_rays();

    /// Evolve all rays until complete
    void _evolve_to_completion();

    

};


