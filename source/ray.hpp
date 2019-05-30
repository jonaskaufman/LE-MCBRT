#ifndef RAY_H
#define RAY_H

#include "parameters.hpp"

#include <iostream>
#include <math.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

typedef std::pair<int, int> PIXEL;

class Ray
{
public:
    Ray() = delete;
    ~Ray();
    Ray(const bool primary,
        const double angle,
        PIXEL current_pixel,
        PIXEL_EDGE current_edge,
        double current_edge_dist,
        double current_energy); // TODO make private

    static Ray primary(const double angle, PIXEL spawn_pixel, PIXEL_EDGE spawn_edge, double spawn_edge_dist);
    static Ray secondary_from_center(const double angle, PIXEL spawn_pixel, double energy);

    /// Trace ray through current pixel and return:
    ///     - distance travelled in pixel
    ///     - i,j for next pixel
    std::pair<double, PIXEL> trace();

    /// Ray activation / correction
    void deactivate();
    bool is_active();
    void set_corrected();

    /// Access data
    PIXEL get_current_pixel();
    double get_current_energy();
    void set_current_energy(double new_energy);
    bool is_primary();
    PIXEL_EDGE get_current_edge();
    double get_current_edge_dist();

    const double m_angle; /// angle of ray, ideally should be private

private:
    bool m_active;              /// whether ray is active
    const bool m_primary;       /// whether ray is primary
    PIXEL m_current_pixel;      /// current pixel
    PIXEL_EDGE m_current_edge;  /// current edge
    double m_current_edge_dist; /// current distance along edge (from TOP or LEFT),
                                // becomes inactive when dist is out of bounds
    double m_current_energy;    /// energy remaining, becomes inactive when 0
    bool m_corrected = true;    /// whether ray is corrected TODO what is this?

    /// Updates meta-data about the ray, such as it's current pixel and edge distance
    PIXEL _update_ray(int delta_x, int delta_y, double a, double b);

    /// String name of given edge
    std::string _get_edge_name(PIXEL_EDGE edge);
};

#endif
