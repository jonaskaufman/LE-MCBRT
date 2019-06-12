#include "ray.hpp"

#include <assert.h>
#include <iostream>
#include <math.h>

Ray::Ray(const bool& primary,
         const double& angle,
         const PIXEL& pixel,
         const PIXEL_EDGE& edge,
         const double& edge_dist,
         const double& energy)
    : m_active(true),
      m_primary(primary),
      m_angle(angle),
      m_current_pixel(pixel),
      m_current_edge(edge),
      m_current_edge_dist(edge_dist),
      m_current_energy(energy)
{
}

Ray Ray::primary(const double angle, PIXEL spawn_pixel, PIXEL_EDGE spawn_edge, double spawn_edge_dist)
{
    return Ray(true, angle, spawn_pixel, spawn_edge, spawn_edge_dist, PARAM_E0);
}

Ray Ray::secondary_from_center(double angle, PIXEL spawn_pixel, double energy)
{
    DEBUG(DB_INIT_SEC, std::cout << "New secondary ray spawning from center of pixel " << spawn_pixel.first << ","
                                 << spawn_pixel.second << " with angle " << angle << " and energy " << energy
                                 << std::endl);
    PIXEL_EDGE new_edge;
    int flip;
    double offset;
    int i_adjust = 0;
    int j_adjust = 0;

    // Starting from center of current pixel,
    // figure out which edge the ray is going to hit
    if (angle > 7 * M_PI / 4 || angle < M_PI / 4) // SOUTH
    {
        new_edge = PIXEL_EDGE::TOP; // bottom edge to top edge
        flip = 1;
        offset = 0;
        j_adjust = 1;
    }
    else if (angle < 3 * M_PI / 4) // EAST
    {
        new_edge = PIXEL_EDGE::LEFT; // right edge to left edge
        flip = -1;
        offset = M_PI / 2;
        i_adjust = 1;
    }
    else if (angle < 5 * M_PI / 4) // NORTH
    {
        new_edge = PIXEL_EDGE::BOTTOM; // top edge to bottom edge
        flip = -1;
        offset = M_PI;
        j_adjust = -1;
    }
    else // if (angle < 2 * M_PI) // WEST
    {
        new_edge = PIXEL_EDGE::RIGHT; // left edge to right edge
        flip = 1;
        offset = 3 * M_PI / 2;
        i_adjust = -1;
    }

    // Update pixel and edge dist accordingly
    PIXEL new_pixel(spawn_pixel.first + i_adjust, spawn_pixel.second + j_adjust);
    double new_edge_dist = 0.5 + 0.5 * flip * tan(angle - offset);
    DEBUG(DB_INIT_SEC, std::cout << "Moved to " << _get_edge_name(new_edge) << " edge of pixel " << new_pixel.first
                                 << "," << new_pixel.second << ", new edge dist is " << new_edge_dist << std::endl);
    assert(new_edge_dist > 0 && new_edge_dist <= 1);
    return Ray(false, angle, new_pixel, new_edge, new_edge_dist, energy);
}

bool Ray::is_primary() { return m_primary; }

bool Ray::is_active() { return m_active; }

void Ray::deactivate()
{
    m_active = false;
    return;
}

PIXEL Ray::get_current_pixel() { return m_current_pixel; }

double Ray::get_current_energy() { return m_current_energy; }

void Ray::set_current_energy(double new_energy)
{
    m_current_energy = new_energy;
    return;
}

/*  Extended tracing description:
 *
 *  Simulates a ray moving from one pixel to another.
 *  Returns distance traveled and pixel traveresed.
 *
 *  Tracing is done by considering the two possible edges that a ray may hit,
 *  based on the quadrant in which its angle lies.
 *
 *  For example, a ray moving SW from the TOP will hit either the LEFT or BOTTOM.
 *  The correct edge can be determined by forming triangles with each candidate.
 *
 *  Notation:
 *  Angles defined with zero pointing SOUTH, increasing COUNTERCLOCKWISE
 *  Pixels and edge distances are measured from NW (TOP LEFT) corner
 *
 *  Special cases:
 *  Rays that pass exactly through the corner of a pixel are not considered,
 *  nor are those that travel perfectly vertically or horizontally.
 *  With randomized rays, these cases will never occur in practice.
 */
std::pair<PIXEL, double> Ray::trace()
{

    // Get local coordinates of ray origin position
    // i.e. horizontal and vertical distance from top-left of current pixel
    std::pair<double, double> local_coords = _get_local_pixel_coordinates();
    double pixel_dist_x = local_coords.first;
    double pixel_dist_y = local_coords.second;

    DEBUG(DB_TRACE, std::cout << "Beginning trace, ray is on " << _get_edge_name(m_current_edge) << " edge of pixel "
                              << m_current_pixel.first << "," << m_current_pixel.second << std::endl);
    DEBUG(DB_TRACE, std::cout << "x,y distance from top left of current pixel is " << pixel_dist_x << ","
                              << pixel_dist_y << std::endl);

    // Variables numbered with 1 are for case of moving to an adjacent horizontal pixel
    // Variables numbered with 2 are for case of moving to an adjacent vertical pixel
    // a and b are the perpendicular sides of the triangle, where a is horizonatl and b is vertical
    // c is the hypotenuse

    double a1;    // Case 1 horizontal
    double b2;    // Case 2 vertical
    double alpha; // The angle between b and c

    // Keeping track of this to eventually determine how to get new edge dist from triangle
    int dir_vert;
    int dir_horiz;

    // Find the two candidate edges that the ray can go to, based on the quadrant of its angle
    // Also determine the correct angle alpha of the two candidate triangles (same for both)
    PIXEL_EDGE candidate1, candidate2;
    if (m_angle < M_PI / 2) // going SE
    {
        // If in this quadrant and not on the correct edge, something went wrong
        assert(m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::LEFT);
        candidate1 = PIXEL_EDGE::RIGHT;
        candidate2 = PIXEL_EDGE::BOTTOM;
        alpha = m_angle;
        dir_vert = 1;
        dir_horiz = 1;
        DEBUG(DB_TRACE, std::cout << "Going SE with angle " << m_angle << std::endl);
    }
    else if (m_angle < M_PI) // going NE
    {
        // If in this quadrant and not on the correct edge, something went wrong
        assert(m_current_edge == PIXEL_EDGE::BOTTOM || m_current_edge == PIXEL_EDGE::LEFT);
        candidate1 = PIXEL_EDGE::RIGHT;
        candidate2 = PIXEL_EDGE::TOP;
        alpha = M_PI - m_angle;
        dir_vert = -1;
        dir_horiz = 1;
        DEBUG(DB_TRACE, std::cout << "Going NE with angle " << m_angle << std::endl);
    }
    else if (m_angle < 3 * M_PI / 2) // going NW
    {
        // If in this quadrant and not on the correct edge, something went wrong
        assert(m_current_edge == PIXEL_EDGE::BOTTOM || m_current_edge == PIXEL_EDGE::RIGHT);
        candidate1 = PIXEL_EDGE::LEFT;
        candidate2 = PIXEL_EDGE::TOP;
        alpha = m_angle - M_PI;
        dir_vert = -1;
        dir_horiz = -1;
        DEBUG(DB_TRACE, std::cout << "Going NW with angle " << m_angle << std::endl);
    }
    else // going SW
    {
        // If in this quadrant and not on the correct edge, something went wrong
        assert(m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::RIGHT);
        candidate1 = PIXEL_EDGE::LEFT;
        candidate2 = PIXEL_EDGE::BOTTOM;
        alpha = 2 * M_PI - m_angle;
        dir_vert = 1;
        dir_horiz = -1;
        DEBUG(DB_TRACE, std::cout << "Going SW with angle " << m_angle << std::endl);
    }

    // Determine known side lengths of candidate triangles
    if (candidate1 == PIXEL_EDGE::LEFT)
    {
        a1 = pixel_dist_x;
    }
    else // RIGHT
    {
        a1 = 1 - pixel_dist_x;
    }

    if (candidate2 == PIXEL_EDGE::TOP)
    {
        b2 = pixel_dist_y;
    }
    else // BOTTOM
    {
        b2 = 1 - pixel_dist_y;
    }

    // Solve the triangles
    double b1 = a1 / tan(alpha);
    double c1 = a1 / sin(alpha);
    double a2 = b2 * tan(alpha);
    double c2 = a2 / sin(alpha);

    // Find shortest hypotenuse and figure out ray adjustments accordingly
    double dist_traveled;
    PIXEL old_pixel = m_current_pixel;
    PIXEL_EDGE new_edge;
    double new_edge_dist;

    int delta_x = 0; // horizontal pixel adjustment
    int delta_y = 0; // vertical pixel adjustment
    if (c1 < c2)
    {
        // Candidate 1: Traveled horizontally
        dist_traveled = c1;

        if (candidate1 == PIXEL_EDGE::LEFT)
        {
            // Travel to the left pixel
            delta_x = -1;
            new_edge = PIXEL_EDGE::RIGHT;
        }
        else
        {
            // Travel to the right pixel
            delta_x = 1;
            new_edge = PIXEL_EDGE::LEFT;
        }
        new_edge_dist = pixel_dist_y + dir_vert * b1;
    }
    else
    {
        // Candidate 2: Travelled vertically
        dist_traveled = c2;

        if (candidate2 == PIXEL_EDGE::TOP)
        {
            // Travel to the top pixel
            delta_y = -1;
            new_edge = PIXEL_EDGE::BOTTOM;
        }
        else
        {
            // Travel to the bottom pixel
            delta_y = 1;
            new_edge = PIXEL_EDGE::TOP;
        }
        new_edge_dist = pixel_dist_x + dir_horiz * a2;
    }

    // Edge distance must be well-defined
    assert(new_edge_dist > 0 && new_edge_dist <= 1);

    PIXEL new_pixel(m_current_pixel.first + delta_x, m_current_pixel.second + delta_y);
    DEBUG(DB_TRACE, std::cout << "Moving to " << _get_edge_name(new_edge) << " edge of pixel " << new_pixel.first << ","
                              << new_pixel.second << ", new edge dist is " << new_edge_dist << std::endl);
    m_current_pixel = new_pixel;
    m_current_edge = new_edge;
    m_current_edge_dist = new_edge_dist;

    // All done
    DEBUG(DB_TRACE, std::cout << "Done" << std::endl << std::endl);
    std::pair<PIXEL, double> result(old_pixel, dist_traveled);
    return result;
}

std::pair<double, double> Ray::_get_local_pixel_coordinates()
{
    double pixel_dist_x, pixel_dist_y;
    if (m_current_edge == PIXEL_EDGE::TOP)
    {
        pixel_dist_x = m_current_edge_dist;
        pixel_dist_y = 0;
    }
    else if (m_current_edge == PIXEL_EDGE::BOTTOM)
    {
        pixel_dist_x = m_current_edge_dist;
        pixel_dist_y = 1;
    }
    else if (m_current_edge == PIXEL_EDGE::LEFT)
    {
        pixel_dist_x = 0;
        pixel_dist_y = m_current_edge_dist;
    }
    else // if (m_current_edge == PIXEL_EDGE::RIGHT)
    {
        pixel_dist_x = 1;
        pixel_dist_y = m_current_edge_dist;
    }

    std::pair<double, double> result(pixel_dist_x, pixel_dist_y);
    return result;
}

std::string Ray::_get_edge_name(PIXEL_EDGE edge)
{
    if (edge == PIXEL_EDGE::TOP)
    {
        return "top";
    }
    else if (edge == PIXEL_EDGE::RIGHT)
    {
        return "right";
    }
    else if (edge == PIXEL_EDGE::BOTTOM)
    {
        return "bottom";
    }
    else
    {
        return "left";
    }
}

