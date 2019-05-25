

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "parameters.hpp"
#endif

#ifndef SIMULATION_H
#define SIMULATION_H
#include "simulation.hpp"
#endif

#include <utility>

class Ray
{
    public:
    Ray() = delete;

    /// Constructor
    Ray(const bool primary, const double angle, std::pair<int, int> current_pixel, SimulationSerial::PIXEL_EDGE current_edge, double current_edge_dist);

    /// Evolve ray by one step, return a tuple of:
    ///     - i,j for pixel visited
    ///     - energy deposited to that pixel
    ///     - new rays generated
    std::tuple<std::pair<int, int>, double, std::vector<Ray>> evolve();

    private:
    bool m_active;                          /// whether ray is active
    const bool m_primary;                   /// whether ray is primary
    const double m_angle;                   /// angle of ray
    std::pair<int, int> m_current_pixel;         /// current pixel
    SimulationSerial::PIXEL_EDGE m_current_edge;  /// current edge
    double m_current_edge_dist;             /// current distance along edge (from TOP or LEFT)

    /// Determine whether primary ray interacts at current pixel
    bool _random_interact();
    
    /// Trace ray through current pixel and return:
    ///     - distance travelled in pixel
    ///     - i,j for next pixel
    std::pair<double, std::pair<int, int>> _trace(); 

    /// Generate secondary rays from primary ray
    std::vector<Ray> _spawn_secondary_rays();

};

