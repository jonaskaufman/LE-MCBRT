#ifndef RAY_H
#define RAY_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "parameters.hpp"

#include <iostream>
#include <math.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

struct Pixel
{
    int first;
    int second;
};

struct Region{
    int first;
    int second;
};

struct Coords
{
    double x;
    double y;
};

struct TraceHistory
{
    Pixel visited;
    double distance;
};

class Ray
{
public:
    CUDA_CALLABLE_MEMBER Ray() = delete;

    /// Named constructor for primary rays
    CUDA_CALLABLE_MEMBER static Ray
    primary(const double angle, Pixel spawn_pixel, PIXEL_EDGE spawn_edge, double spawn_edge_dist, Region region);

    /// Named constructor for secondary rays, originating from center of given pixel
    CUDA_CALLABLE_MEMBER static Ray secondary_from_center(const double angle, Pixel spawn_pixel, double energy);

    /// Check whether ray is primary or secondary
    CUDA_CALLABLE_MEMBER bool is_primary();

    /// Check or set ray activation
    CUDA_CALLABLE_MEMBER bool is_active();
    CUDA_CALLABLE_MEMBER void deactivate();

    /// Access ray position
    CUDA_CALLABLE_MEMBER Pixel get_current_pixel();

    /// Access or set ray energy
    CUDA_CALLABLE_MEMBER double get_current_energy();
    CUDA_CALLABLE_MEMBER void set_current_energy(double new_energy);

    CUDA_CALLABLE_MEMBER Region get_current_region();

    /// Trace ray through current pixel, return distance travelled in pixel and pixel visited
    CUDA_CALLABLE_MEMBER TraceHistory trace();

    /// String name of a given edge
    static std::string get_edge_name(PIXEL_EDGE edge);



private:
    /// Full parameterized constructor
    CUDA_CALLABLE_MEMBER Ray(const bool& primary,
        const double& angle,
        const Pixel& pixel,
        const PIXEL_EDGE& edge,
        const double& edge_dist,
        const Region& region,
        const double& energy);

        CUDA_CALLABLE_MEMBER Ray(const bool& primary,
            const double& angle,
            const Pixel& pixel,
            const PIXEL_EDGE& edge,
            const double& edge_dist,
            const double& energy);


    bool m_active;              /// whether ray is active
    bool m_primary;             /// whether ray is primary (had to make non const for operator=)
    double m_angle;             /// ray angle (had to make non const for operator=)
    Pixel m_current_pixel;      /// current pixel
    PIXEL_EDGE m_current_edge;  /// current edge
    double m_current_edge_dist; /// current distance along edge from TOP/LEFT of current pixel
    double m_current_energy;    /// energy remaining
    Region m_current_region;     /// current region in the grid

    /// Local x,y coordinates relative to TOP LEFT of current pixel
    CUDA_CALLABLE_MEMBER Coords _get_local_pixel_coordinates();

    
};

#endif
