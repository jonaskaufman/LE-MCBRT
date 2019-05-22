/// Global definitions (maybe make these constructor arguments?)
#define N 1000      /// pixels per side
#define E0 1.0      /// initial primary ray energy
#define d 1.0       /// source distance
#define sigma 1.0   /// source spread
#define a 1.0       /// primary ray interaction probability scaling
#define F 0.5       /// primary ray interaction energy fraction
#define G 0.1       /// secondary ray deposition constant


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


class Ray
{
    public:
    /// Evolve ray by one step, return a tuple of:
    ///     - i,j for pixel visited
    ///     - energy deposited to that pixel
    ///     - new rays generated
    std::tuple<std::pair<int>, double, std::vector<Ray>> evolve();

    private:
    bool m_active;                          /// whether ray is active
    const bool m_primary;                   /// whether ray is primary
    const couble m_angle;                   /// angle of ray
    std::pair<int> m_current_pixel;         /// current pixel
    Simulation::PIXEL_EDGE m_current_edge;  /// current edge
    double m_current_edge_dist;             /// current distance along edge (from TOP or LEFT)

    /// Determine whether primary ray interacts at current pixel
    bool _interact();
    
    /// Trace ray through current pixel and return:
    ///     - distance travelled in pixel
    ///     - i,j for next pixel
    std::pair<double, std::pair<int>> _trace(); 

}

