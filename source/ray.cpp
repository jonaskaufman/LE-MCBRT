#include "ray.hpp"

Ray::Ray(const bool primary,
         const double angle,
         PIXEL current_pixel,
         PIXEL_EDGE current_edge,
         double current_edge_dist,
         double current_energy)
    : m_primary(primary), m_angle(angle)
{

    m_active = true;
    m_current_pixel = current_pixel;
    m_current_edge = current_edge;
    m_current_edge_dist = current_edge_dist;
    m_current_energy = current_energy;
} // TODO initialization list

Ray Ray::primary(const double angle, PIXEL spawn_pixel, PIXEL_EDGE spawn_edge, double spawn_edge_dist)
{
    return Ray(true, angle, spawn_pixel, spawn_edge, spawn_edge_dist, PARAM_E0);
}

Ray Ray::secondary_from_center(double angle, PIXEL spawn_pixel, double energy)
{
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

    assert(new_edge_dist > 0 && new_edge_dist <= 1);

    DEBUG(DB_TRACE, printf("New secondary ray edge dist %.5f\n", new_edge_dist));

    // TODO need to check edges for out of bounds, could do in simulation class
    return Ray(false, angle, new_pixel, new_edge, new_edge_dist, energy);
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
std::pair<double, PIXEL> Ray::trace()
{

    // Get local coordinates of ray origin position
    // i.e. horizontal and vertical distance from top-left of current pixel
    std::pair<double, double> local_coords = _get_local_pixel_coordinates();
    double pixel_dist_x = local_coords.first;
    double pixel_dist_y = local_coords.second;

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
        DEBUG(DB_TRACE, printf("Going SE. angle, alpha: %.5f %.5f\n", m_angle, alpha));
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
        DEBUG(DB_TRACE, printf("Going NE. angle, alpha: %.5f %.5f\n", m_angle, alpha));
    }
    else if (m_angle < 3 * M_PI / 2) // going NW
    {
        DEBUG(DB_TRACE, std::cout << "Going NW but edge is: " << _get_edge_name(m_current_edge) << std::endl);
        // If in this quadrant and not on the correct edge, something went wrong
        assert(m_current_edge == PIXEL_EDGE::BOTTOM || m_current_edge == PIXEL_EDGE::RIGHT);
        candidate1 = PIXEL_EDGE::LEFT;
        candidate2 = PIXEL_EDGE::TOP;
        alpha = m_angle - M_PI;
        dir_vert = -1;
        dir_horiz = -1;
        DEBUG(DB_TRACE, printf("Going NW. angle, alpha: %.5f %.5f\n", m_angle, alpha));
    }
    else // if (m_angle < 2 * M_PI) // going SW
    {
        // If in this quadrant and not on the correct edge, something went wrong
        assert(m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::RIGHT);
        candidate1 = PIXEL_EDGE::LEFT;
        candidate2 = PIXEL_EDGE::BOTTOM;
        alpha = 2 * M_PI - m_angle;
        dir_vert = 1;
        dir_horiz = -1;
        DEBUG(DB_TRACE, printf("Going SW. angle, alpha: %.5f %.5f\n", m_angle, alpha));
    }
    // else
    // {
    //     std::cout << "unexpected angle to trace" << std::endl;
    // }

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

    DEBUG(DB_TRACE, printf("pixel_dist_x: %.5f\n", pixel_dist_x));
    DEBUG(DB_TRACE, printf("pixel_dist_y: %.5f\n", pixel_dist_y));
    DEBUG(DB_TRACE, printf("old edge distance: %.5f\n", m_current_edge_dist));

    int delta_x = 0; // horizontal pixel adjustment
    int delta_y = 0; // vertical pixel adjustment
    if (c1 < c2)
    {
        // Candidate 1: Traveled horizontally
        dist_traveled = c1;
        DEBUG(DB_TRACE, printf("a1, b1: %.5f %.5f\n", a1, b1));
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
        DEBUG(DB_TRACE, printf("a2, b2: %.5f %.5f\n", a2, b2));
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

    DEBUG(DB_TRACE, printf("new_edge_dist: %.5f\n", new_edge_dist));
    // Edge distance must be well-defined
    assert(new_edge_dist > 0 && new_edge_dist <= 1);

    // Update ray
    DEBUG(DB_TRACE, printf("delta x, y: %d %d\n", delta_x, delta_y));

    PIXEL new_pixel(m_current_pixel.first + delta_x, m_current_pixel.second + delta_y);
    m_current_pixel = new_pixel;
    m_current_edge = new_edge;
    m_current_edge_dist = new_edge_dist;

    DEBUG(DB_TRACE, printf("Traveled: %.5f\n", dist_traveled));
    // All done
    std::pair<double, PIXEL> result(dist_traveled, old_pixel);
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

//Ray::~Ray() { DEBUG(DB_TRACE, "Ray is getting deleted"); }

void Ray::deactivate() { m_active = false; }

bool Ray::is_active() { return m_active; }

PIXEL Ray::get_current_pixel() { return m_current_pixel; }

double Ray::get_current_energy() { return m_current_energy; }

void Ray::set_current_energy(double new_energy) { m_current_energy = new_energy; }

bool Ray::is_primary() { return m_primary; }

PIXEL_EDGE Ray::get_current_edge() { return m_current_edge; }

double Ray::get_current_edge_dist() { return m_current_edge_dist; }

