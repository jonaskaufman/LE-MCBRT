/// Global definitions (maybe make these constructor arguments?)
#define N 1000      /// pixels per side
#define E0 1.0      /// initial primary ray energy
#define d 1.0       /// source distance
#define sigma 1.0   /// source spread


class SimulationSerial
{
    public: 
    /// Label the pixel edges
    enum class PIXEL_EDGE {TOP, BOTTOM, LEFT, RIGHT};

    /// Initalize densities
    void initialize_densities_random();
    void initialize_densities_constant(const double density);

    /// Run simulation for a given number of primary rays
    void run(int num_primary_rays);

    /// Get array of doses
    double* get_doses;
    
    private:
    const double m_densities[N][N];         /// pixel density values
    double m_doses[N][N];                   /// pixel dose values
    std::vector<Ray> m_rays;                /// active rays

    /// Generate new primary ray
    void _spawn_primary_ray();

    /// Evolve all rays by one step
    void _evolve_rays();

    /// Evolve all rays until complete
    void _evolve_to_completion();

}
