#include "ray.hpp"
/*
Ray::Ray() : m_primary(false), m_angle(0.0){
    std::cout << "default" << std::endl;
}
*/

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
}

Ray Ray::primary(const double angle, PIXEL spawn_pixel, PIXEL_EDGE spawn_edge, double spawn_edge_dist)
{
    return Ray(true, angle, spawn_pixel, spawn_edge, spawn_edge_dist, E0);
}


Ray Ray::secondary_from_center(double angle, PIXEL spawn_pixel, double energy)
{
    PIXEL_EDGE new_edge;
    int flip;
    double offset;
    int i_adjust = 0;
    int j_adjust = 0;

    // TODO Possible to change these to single comparison? by going in correct order 
    // Starting from center of current pixel,
    // figure out which edge the ray is going to hit
    if (angle > 7 * M_PI / 4 && angle < M_PI / 4)   // SOUTH
    {
        new_edge = PIXEL_EDGE::TOP; // bottom edge to top edge
        flip = 1;
        offset = 0;
        j_adjust = -1;
    }
    else if (angle > M_PI / 4 && angle < 3 * M_PI / 4)  // EAST
    {
        new_edge = PIXEL_EDGE::LEFT; // right edge to left edge
        flip = -1;
        offset = M_PI / 2;
        i_adjust = 1;
    }
    else if (angle > 3 * M_PI / 4 && angle < 5 * M_PI / 4)  // NORTH
    {
        new_edge = PIXEL_EDGE::BOTTOM; // top edge to bottom edge
        flip = -1;
        offset = M_PI;
        j_adjust = 1;
    }
    else // if (angle > 5 * M_PI / 4 && angle < 2 * M_PI) // WEST
    {
        new_edge = PIXEL_EDGE::RIGHT; // left edge to right edge
        flip = 1;
        offset = 3 * M_PI / 2;
        i_adjust = -1;
    }

    // Update pixel and edge dist accordingly
    PIXEL new_pixel(spawn_pixel.first + i_adjust, spawn_pixel.second + j_adjust);
    double new_edge_dist = 0.5 + 2 * flip * tan(angle - offset);

    // TODO need to check edges for out of bounds, could do in simulation class
    return Ray(false, angle, new_pixel, new_edge, new_edge_dist, energy);
}


/* Simulates a ray moving from one pixel to another.
    Returns distance traveled and pixel traveresed
.
   This is done by taking it's source angle and deciding if it's moving SW, NW, NE, or SE.
   Recall that the source angle is normalized to be between 0 and 2 pi, so direction is
   based off unit circle quadrants (rotated so that 0 rad is going South).

   Next part is better understood with example. Say ray is moving SW from the top.
   Then it will either cross the left edge, or the bottom edge. So both edges are considered,
   triangles leading to both are formed and solved. The one with shortest hyponetuse is correct.
*/


// Break into helper functions wherever possible!
std::pair<double, PIXEL> Ray::trace()
{

    double pixel_x, pixel_y; // horizontal and vertical displacement from a top-left edge of current pixel

    double a1; // When solving any triangle, horizonal legs are labeled as a1, a2
               // all variables labeled _1 are considering net horizontal movement

    double b2; // When solving any triangle, vertical legs are labeled as b1, b2
               // all variables labeled _2 are considering net vertical movement
    
    double theta;
    int delta_x1, delta_x2, delta_y1, delta_y2; // TODO initialize to zero
    double edge_dist1, edge_dist2;

    // Get local coordinates of ray origin
    if (m_current_edge == PIXEL_EDGE::TOP)
    {
        pixel_x = m_current_edge_dist;
        pixel_y = 0;
    }
    else if (m_current_edge == PIXEL_EDGE::BOTTOM)
    {
        pixel_x = m_current_edge_dist;
        pixel_y = 1;
    }
    else if (m_current_edge == PIXEL_EDGE::LEFT)
    {
        pixel_x = 0;
        pixel_y = m_current_edge_dist;
    }
    else
    {
        pixel_x = 1;
        pixel_y = m_current_edge_dist;
    }

    // For given quadrant, find relevant triangle sides etc. in each case
    if (m_angle < M_PI / 2) // going SE
    { 
        // If in this quadrant and not on the correct edge, something went wrong    
        assert( m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::LEFT );
        theta = m_angle;
     
        // Case 1: Ray hits bottom wall
        b1 = 1 - pixel_y;
        delta_y1 = 1;

        // Case 2: Ray hits right wall
        a2 = 1 - pixel_x;
        delta_x2 = 1;

    }
    else if (m_angle < M_PI) // going NE
    {
        // If in this quadrant and not on the correct edge, something went wrong        
        assert( m_current_edge == PIXEL_EDGE::BOTTOM || m_current_edge == PIXEL_EDGE::LEFT );
        theta = M_PI - m_angle;
        
        // Case 1: Ray hits top wall
        b1 = pixel_y;
        delta_y1 = -1;

        // Case 2: Ray hits right wall
        a2 = 1 - pixel_x;
        delta_x2 = 1;
    }
    else if (m_angle < 3 * M_PI / 2) // going NW
    {
        // If in this quadrant and not on the correct edge, something went wrong
        assert( m_current_edge == PIXEL_EDGE::BOTTOM || m_current_edge == PIXEL_EDGE::RIGHT );
        theta = m_angle - M_PI;
        
        // Case 1: Ray hits top wall
        b1 = pixel_y;        
        delta_y1 = -1;

        // Case 2: Ray hits left wall
        a2 = pixel_x;
        delta_x2 = -1;
    
    }
    else // if (m_angle < 2 * M_PI) // going SW
    {
        // If in this quadrant and not on the correct edge, something went wrong
        assert( m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::RIGHT );
        
        //TODO: what is theta?
        
        // Case 1: Ray hits bottom wall
        b1 = 1 - pixel_y;
        delta_y1 = 1;
        
        // Case 2: Ray hits left wall;
        a2 = pixel_x;
        delta_x2 = -1;

    }
   // else    
   // {
   //     std::cout << "unexpected angle to trace" << std::endl;
   // }

    // Solve the triangles
    

    // c1 = a1 / sin theta
    // c2 = b2 / cos theta
            
    double b1 = a1 * tan(alpha1);
    double c1 = sqrt(a1 * a1 + b1 * b1);
    DEBUG(DB_TRACE, printf("a1, b1, c1: %.2f, %.2f, %.2f\n", a1, b1, c1));

    double a2 = b2 * tan(alpha2);
    double c2 = sqrt(a2 * a2 + b2 * b2);
    DEBUG(DB_TRACE, printf("a2, b2, c2: %.2f, %.2f, %.2f\n", a2, b2, c2));

    double dist_traveled;
    PIXEL old_pixel;


    // TODO get rid of update ray stuff
    // shortest hypotenuse determines next pixel and edge
    if (c1 < c2)
    {
        dist_traveled = c1;
        old_pixel = _update_ray(delta_x1, delta_y1, a1, b1);
    }
    else
    {
        dist_traveled = c2;
        old_pixel = _update_ray(delta_x2, delta_y2, a2, b2);
    }
    
    if (pixel_x < 0 || pixel_y < 0)
    {
        printf("error. pixel x and y should be non-negative\n");
        exit(1);
    }
    if (pixel_x > 1 || pixel_y > 1)
    {
        printf("error. pixel x and y should be smaller than 1\n");
        exit(1);
    }

    std::pair<double, PIXEL> result(dist_traveled, old_pixel);
    return result;
}


/*
 Updates meta-data for the current ray. This includes current pixel, current edge,
   and current edge distance. Returns pixel it traversed

   Here, delta_x and delta_y describe pixel offset. So (-1, 0) implies movement to the left neighbor pixel
   a, b are legs of the triangle describing the rays movement. So (0.3, 0.4) implies the ray moved 0.3 horizontal units
   and 0.4 vertical units a, b are always positive

PIXEL Ray::_update_ray(int delta_x, int delta_y, double a, double b)
{
    // get next pixel
    int next_x = m_current_pixel.first + delta_x;
    int next_y = m_current_pixel.second + delta_y;

    PIXEL old_pixel(m_current_pixel.first, m_current_pixel.second);
    PIXEL next_pixel(next_x, next_y);

    std::string old_edge = _get_edge_name(m_current_edge);
    DEBUG(DB_TRACE, printf("old edge distance: %.2f\n", m_current_edge_dist));
    // update edge
    if (delta_x != 0)
    { // moving horizontally to left/right edge
        if (m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::BOTTOM)
        {

            if (m_angle > M_PI / 2 && m_angle < 3 * M_PI / 2)
            { // going up
                m_current_edge_dist = -1.0 * b;
                if (m_current_edge_dist < next_pixel.second)
                {
                    m_current_edge_dist += 1;
                }
                
                if (m_current_pixel.second == 0){ // don't go negative for the 0th row. Just fractional
                    m_current_edge_dist += 1;
                }
                
            }
            else
            { // going down
                m_current_edge_dist = b;
            }
        }

        else
        { // already on left/right, ending up on left/right (with some up/down displacement)
            if (m_angle > M_PI / 2 && m_angle < 3 * M_PI / 2)
            { // going up
                m_current_edge_dist -= b;
                DEBUG(DB_TRACE, printf("decreased edge distance by b: %.2f\n", b));
            }
            else
            { // going down
                m_current_edge_dist += b;
                DEBUG(DB_TRACE, printf("increased edge distance by b: %.2f\n", b));
            }
        }

        if (delta_x == 1)
        {
            m_current_edge = PIXEL_EDGE::RIGHT; // set right edge
        }
        else
        {
            m_current_edge = PIXEL_EDGE::LEFT; // set left edge
        }
    }
    else
    { // moving vertically
        if (m_current_edge == PIXEL_EDGE::LEFT || m_current_edge == PIXEL_EDGE::RIGHT)
        {
            if (m_angle < M_PI)
            { // going left
                m_current_edge_dist = - a;

               //  if (m_corrected == false && m_primary == false && m_current_edge_dist < next_x)
               // { // correct for secondary ray going opposite primary
               //     m_current_edge_dist += 1;
              //      DEBUG(DB_TRACE, printf("adjusted current edge distance by +1\n"));
              //  }
            }
            else
            { // going right
                m_current_edge_dist = m_current_pixel.first + a;
            }
        }

        else
        { // already on top/bottom, ending up on top/bottom (with some left/right displacement)
            if (m_angle < M_PI)
            { // going left
                m_current_edge_dist -= a;
                DEBUG(DB_TRACE, printf("decreased edge distance by a: %.2f\n", a));
            }
            else
            { // going right
                m_current_edge_dist += a;
                DEBUG(DB_TRACE, printf("increased edge distance by a: %.2f\n", a));
            }
        }

        if (delta_y == 1)
        {
            //if (m_current_edge == PIXEL_EDGE::LEFT && m_angle < M_PI && m_corrected == true)
           // {
            ///    m_current_edge_dist += 1; // special case where the formula gives incorrect result b/c pixel is same
             //   DEBUG(DB_TRACE, printf("special case, increased current edge distance by 1\n"));
           // }
            //m_current_edge = PIXEL_EDGE::BOTTOM; // set bottom edge
        }
        else
        {

        ////    if (m_current_edge == PIXEL_EDGE::LEFT && m_angle < M_PI && m_corrected == true)
        //    {
        //        m_current_edge_dist += 1; // special case where the formula gives incorrect result b/c pixel is same
        //        DEBUG(DB_TRACE, printf("special case, increased current edge distance by 1\n"));
         //   }

            m_current_edge = PIXEL_EDGE::TOP; // set top edge
        }
 
    }
    std::string new_edge = _get_edge_name(m_current_edge);
    DEBUG(DB_TRACE, std::cout << "went from " << old_edge << " edge to " << new_edge << " edge" << std::endl);
    DEBUG(DB_TRACE, printf("new edge distance is %.2f\n", m_current_edge_dist));
    
    DEBUG(DB_TRACE, printf("Going from %d, %d to %d, %d\n", m_current_pixel.first, m_current_pixel.second,
                           next_pixel.first, next_pixel.second));

    // Segfault 11 occurs after this print statement

    // update current pixel
    m_current_pixel = next_pixel;

    return old_pixel;
}
*/

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

Ray::~Ray() { DEBUG(DB_TRACE, "Ray is getting deleted"); }

void Ray::deactivate() { m_active = false; }

bool Ray::is_active() { return m_active; }

PIXEL Ray::get_current_pixel() { return m_current_pixel; }

double Ray::get_current_energy() { return m_current_energy; }

void Ray::set_current_energy(double new_energy) { m_current_energy = new_energy; }

bool Ray::is_primary() { return m_primary; }

PIXEL_EDGE Ray::get_current_edge() { return m_current_edge; }

double Ray::get_current_edge_dist() { return m_current_edge_dist; }

