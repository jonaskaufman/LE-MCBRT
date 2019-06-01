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

    /// Named constructor for primary rays
    static Ray primary(const double angle, PIXEL spawn_pixel, PIXEL_EDGE spawn_edge, double spawn_edge_dist);

    /// Named constructor for secondary rays, originating from center of given pixel
    static Ray secondary_from_center(const double angle, PIXEL spawn_pixel, double energy);

    /// Trace ray through current pixel, return distance travelled in pixel and pixel visited TODO switch order?
    std::pair<double, PIXEL> trace();

    /// Check whether ray is primary or secondary
    bool is_primary();

    /// Check or set ray activation
    bool is_active();
    void deactivate();

    /// Access ray position data TODO check if all of these are used
    PIXEL get_current_pixel();
    PIXEL_EDGE get_current_edge();
    double get_current_edge_dist();

    /// Access or set ray energy
    double get_current_energy();
    void set_current_energy(double new_energy);

private:
    /// Full parameterized constructor
    Ray(const bool primary,
        const double angle,
        PIXEL current_pixel,
        PIXEL_EDGE current_edge,
        double current_edge_dist,
        double current_energy);

    bool m_active;              /// whether ray is active
    const bool m_primary;       /// whether ray is primary
    const double m_angle;       /// ray angle
    PIXEL m_current_pixel;      /// current pixel
    PIXEL_EDGE m_current_edge;  /// current edge
    double m_current_edge_dist; /// current distance along edge (from TOP or LEFT of current pixel),
    double m_current_energy;    /// energy remaining

    /// String name of a given edge
    std::string _get_edge_name(PIXEL_EDGE edge);

    /// Local x,y coordinates relative to TOP LEFT of current pixel
    std::pair<double, double> _get_local_pixel_coordinates();
};

#endif
