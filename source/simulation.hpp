#include parameters.hpp

class SimulationSerial
{
    public: 
    // SimulationSerial() = delete;

    // SimulationSerial(const int N); /// will this work? for static arrays?

    /// Label the pixel edges
    enum class PIXEL_EDGE {TOP, BOTTOM, LEFT, RIGHT};

    /// Initalize densities
    void initialize_densities_random();
    void initialize_densities_constant(const double density);

    /// Run simulation for a given number of primary rays
    void run(int num_primary_rays);

    /// Get array of doses
    double* get_doses();
    
    private:
    const double m_densities[N][N];         /// pixel density values
    double m_doses[N][N];                   /// pixel dose values
    std::vector<Ray> m_rays;                /// active rays

    /// Randomly sample source angle for primary rays
    double _random_source_angle();

    /// Generate new primary ray
    void _spawn_primary_ray();

    /// Evolve all rays by one step
    void _evolve_rays();

    /// Evolve all rays until complete
    void _evolve_to_completion();

}
