#ifndef RAY_H
#define RAY_H

#include "parameters.hpp"

#include <utility>
#include <vector>
#include <iostream>
#include <tuple>
#include <math.h>
#include <string>

typedef std::pair<int, int> PIXEL;

class Ray
{
    public:
    
    Ray() = delete;
    ~Ray();

    /// Constructor
    Ray(const bool primary, const double angle, PIXEL current_pixel, \
        PIXEL_EDGE current_edge, double current_edge_dist, double current_energy);

    /// Evolve ray by one step, return a tuple of:
    ///     - i,j for pixel visited
    ///     - energy deposited to that pixel
    ///     - new rays generated
    std::tuple<PIXEL, double, std::vector<Ray>> evolve();

    /// Trace ray through current pixel and return:
    ///     - distance travelled in pixel
    ///     - i,j for next pixel
    std::pair<double, PIXEL> _trace(); 

    void deactivate();
    bool is_active();
    PIXEL get_current_pixel();
    double get_current_energy();
    void set_current_energy(double new_energy);
    bool is_primary();
    PIXEL_EDGE get_current_edge();
    double get_current_edge_dist();

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

    PIXEL _update_ray(int delta_x, int delta_y, double a, double b);
    std::string _get_edge_name(PIXEL_EDGE edge); // gets string name of edge
    
    

};

#endif