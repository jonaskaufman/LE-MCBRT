#ifndef RAY_H
#define RAY_H

#include "parameters.hpp"

#include <utility>
#include <vector>
#include <iostream>
#include <tuple>
#include <math.h>

typedef std::pair<int, int> PIXEL;
typedef std::pair<double, double> POS;


class Ray
{
    public:
    
    Ray();

    /// Constructor
    Ray(const bool primary, const double angle, PIXEL current_pixel, \
        PIXEL_EDGE current_edge, double current_edge_dist, double current_energy);

    /// Evolve ray by one step, return a tuple of:
    ///     - i,j for pixel visited
    ///     - energy deposited to that pixel
    ///     - new rays generated
    std::tuple<PIXEL, double, std::vector<Ray>> evolve();
    const double m_angle;
    private:
    bool m_active;                          /// whether ray is active
    const bool m_primary;                   /// whether ray is primary
                       /// angle of ray
    PIXEL m_current_pixel;                  /// current pixel
    PIXEL_EDGE m_current_edge;              /// current edge
    double m_current_edge_dist;             /// current distance along edge (from TOP or LEFT), 
                                                // becomes inactive when dist is out of bounds
    double m_current_energy;                        // energy remaining, becomes inactive when 0

    /// Determine whether primary ray interacts at current pixel
    bool _random_interact();
    
    /// Trace ray through current pixel and return:
    ///     - distance travelled in pixel
    ///     - i,j for next pixel
    std::pair<double, PIXEL> _trace(); 
    std::pair<double, PIXEL> _tb_trace();
    std::pair<double, PIXEL> _lr_trace();

    /// Generate secondary rays from primary ray
    std::vector<Ray> _spawn_secondary_rays();

};

#endif