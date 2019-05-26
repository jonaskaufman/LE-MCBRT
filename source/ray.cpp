#include "ray.hpp"
/*
Ray::Ray() : m_primary(false), m_angle(0.0){
    std::cout << "default" << std::endl;
}
*/

Ray::Ray(const bool primary, const double angle, std::pair<int, int> current_pixel, \
        PIXEL_EDGE current_edge, double current_edge_dist, double current_energy) : m_primary(primary), m_angle(angle)
{
    m_active = true;
    m_current_pixel = current_pixel;
    m_current_edge = current_edge;
    m_current_edge_dist = current_edge_dist;
    m_current_energy = current_energy;
    
}
/* Moves ray one pixel. 
   Obsolete???????
   

*/
std::tuple<PIXEL, double, std::vector<Ray>> Ray::evolve(){
    std::pair<double, PIXEL> trace = _trace();
    double distance_traveled = trace.first;
    PIXEL next_pixel = trace.second;
    
    std::tuple<PIXEL, double, std::vector<Ray>> result;

    std::vector<Ray> new_rays;
    result = std::make_tuple(next_pixel, distance_traveled, new_rays);
    return result;
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
std::pair<double, PIXEL> Ray::_trace()
{
    
    double pixel_x, pixel_y; // horizontal and vertical displacement from a top-left edge of current pixel

    double a1;  // When solving any triangle, horizonal legs are labeled as a1, a2
                // all variables labeled _1 are considering net horizontal movement

    double b2;  // When solving any triangle, vertical legs are labeled as b1, b2
                // all variables labeled _2 are considering net vertical movement
                

    if (m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::BOTTOM){ // on top/bottom edge
        pixel_x = m_current_edge_dist - m_current_pixel.first;
        pixel_y = m_current_pixel.second;
    }
    else{                                                                            // on left/right edge
        pixel_x = m_current_pixel.first;
        pixel_y = m_current_edge_dist - m_current_pixel.second;
    }
    
    DEBUG(DB_TRACE, printf("pixel x,y: %.2f, %.2f\n", pixel_x, pixel_y));
    double alpha1, alpha2;
    int delta_x1, delta_x2, delta_y1, delta_y2;

    if (m_angle < M_PI / 2){ // going SW
        DEBUG(DB_TRACE, printf("moving SW. Angle: %.2f\n", m_angle * 180 / M_PI))
        alpha1 = M_PI / 2 - m_angle; // hits left wall
        alpha2 = m_angle; // hits bottom wall
        delta_x1 = -1, delta_y1 = 0;
        delta_x2 = 0, delta_y2 = 1;

        if (m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::BOTTOM){ // on top/bottom edge
            a1 = pixel_x;
            b2 = 1;
        }
        else{                                                                           // on left/right edge
            a1 = 1;
            b2 = 1 - pixel_y;
        }
    }
    else if (m_angle < M_PI){ // going NW
        DEBUG(DB_TRACE, printf("moving NW. Angle: %.2f\n", m_angle * 180 / M_PI))
        alpha1 = m_angle - M_PI / 2; // hits left wall
        alpha2 = M_PI - m_angle; // hits top wall
        delta_x1 = -1; delta_y1 = 0;
        delta_x2 = 0; delta_y2 = -1;

        if (m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::BOTTOM){ // on top/bottom edge
            a1 = pixel_x;
            b2 = 1;
        }
        else{                                                                           // on left/right edge
            a1 = 1;
            b2 = pixel_y;
        }
    }
    else if (m_angle < 3 * M_PI / 2){ // going NE
        DEBUG(DB_TRACE, printf("moving NE. Angle: %.2f\n", m_angle * 180 / M_PI))
        alpha1 = 3 * M_PI / 2 - m_angle; // hits right wall
        alpha2 = m_angle - M_PI; // hits top wall
        delta_x1 = 1; delta_y1 = 0;
        delta_x2 = 0; delta_y2 = -1;

        if (m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::BOTTOM){ // on top/bottom edge
            a1 = 1 - pixel_x;
            b2 = 1;
        }
        else{                                                                           // on left/right edge
            a1 = 1;
            b2 = pixel_y;
        }

    }
    else if (m_angle < 2 * PI){ // going SE
        DEBUG(DB_TRACE, printf("moving SE. Angle: %.2f\n", m_angle * 180 / M_PI))
        alpha1 = m_angle - 3 * M_PI / 2; // hits right wall
        alpha2 = 2 * M_PI - m_angle; // hits bottom wall
        delta_x1 = 1, delta_y1 = 0;
        delta_x2 = 0; delta_y2 = 1;

        if (m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::BOTTOM){ // on top/bottom edge
            a1 = 1 - pixel_x;
            b2 = 1;
        }
        else{                                                                           // on left/right edge
            a1 = 1;
            b2 = 1 - pixel_y;
        }
    }
    else{
        std::cout << "unexpected angle to trace" << std::endl;
    }

    // solve the triangles
    double b1 = a1 * tan(alpha1);
    double c1 = sqrt(a1 * a1 + b1 * b1);
    DEBUG(DB_TRACE, printf("a1, b1, c1: %.2f, %.2f, %.2f\n", a1, b1, c1));

    double a2 = b2 * tan(alpha2);
    double c2 = sqrt(a2 * a2 + b2 * b2);
    DEBUG(DB_TRACE, printf("a2, b2, c2: %.2f, %.2f, %.2f\n", a2, b2, c2));

    double dist_traveled;
    PIXEL old_pixel;

    // shortest hypotenuse determines next pixel and edge
    if (c1 < c2){
        dist_traveled = c1;
        old_pixel = _update_ray(delta_x1, delta_y1, a1, b1);
    } 
    else{ 
        dist_traveled = c2;
        old_pixel = _update_ray(delta_x2, delta_y2, a2, b2);
    }

    std::pair<double, PIXEL> result(dist_traveled, old_pixel);
    return result;
}

/* Updates meta-data for the current ray. This includes current pixel, current edge,
   and current edge distance. Returns pixel it traversed

   Here, delta_x and delta_y describe pixel offset. So (-1, 0) implies movement to the left neighbor pixel
   a, b are legs of the triangle desribing the rays movement. So (-0.3, 0.4) implies the ray moved 0.3 units left and 0.4 down.
*/
PIXEL Ray::_update_ray(int delta_x, int delta_y, double a, double b){
    // get next pixel
    int next_x = m_current_pixel.first + delta_x;
    int next_y = m_current_pixel.second + delta_y;

    PIXEL old_pixel(m_current_pixel.first, m_current_pixel.second);
    PIXEL next_pixel(next_x, next_y);
    DEBUG(DB_TRACE, printf("Going from %d, %d to %d, %d\n", \
        m_current_pixel.first, m_current_pixel.second, next_pixel.first, next_pixel.second)); 
    
    std::string old_edge = _get_edge_name(m_current_edge);
    DEBUG(DB_TRACE, printf("old edge distance: %.2f\n", m_current_edge_dist));
    // update edge
    if (delta_x != 0){                      // moving horizontally to left/right edge
        if (m_current_edge == PIXEL_EDGE::TOP){
            m_current_edge_dist = m_current_pixel.second + b;
            
        }
        else if (m_current_edge == PIXEL_EDGE::BOTTOM){
            m_current_edge_dist = m_current_pixel.second - b;
        }
   
        else{      // already on left/right, ending up on left/right (with some up/down displacement)                                              
            if (m_angle > M_PI / 2 && m_angle < 3 * M_PI / 2 ){ // going up
                m_current_edge_dist -= b;
                DEBUG(DB_TRACE, printf("decreased edge distance by b: %.2f\n", b));
            }
            else{                                               // going down
                m_current_edge_dist += b;
                DEBUG(DB_TRACE, printf("increased edge distance by b: %.2f\n", b));
            }
        }
        
        if (delta_x == 1){
            m_current_edge = PIXEL_EDGE::RIGHT; // set right edge
        }
        else{
            m_current_edge = PIXEL_EDGE::LEFT; // set left edge
        }
    }
    else{                                   // moving vertically
        if (m_current_edge == PIXEL_EDGE::LEFT){
            m_current_edge_dist = m_current_pixel.first + a;
        }
        else if (m_current_edge == PIXEL_EDGE::RIGHT){
            m_current_edge_dist = m_current_pixel.second - a;
        }
   
        else{       // already on top/bottom, ending up on top/bottom (with some left/right displacement)                                             
            if (m_angle < M_PI ){                           // going left
                m_current_edge_dist -= a;
                DEBUG(DB_TRACE, printf("decreased edge distance by a: %.2f\n", a));
            }
            else{                                           // going right
                m_current_edge_dist += a;
                DEBUG(DB_TRACE, printf("increased edge distance by a: %.2f\n", a));
            }
        }
        
        if (delta_y == 1){
            m_current_edge = PIXEL_EDGE::BOTTOM; // set bottom edge
        }
        else{
            m_current_edge = PIXEL_EDGE::TOP; // set top edge
        }
    }
    std::string new_edge = _get_edge_name(m_current_edge);
    DEBUG(DB_TRACE, std::cout << "went from " << old_edge << " edge to " << new_edge << " edge" << std::endl);
    DEBUG(DB_TRACE, printf("new edge distance is %.2f\n", m_current_edge_dist));
    
    // update current pixel
    m_current_pixel = next_pixel;

    return old_pixel;
}

std::string Ray::_get_edge_name(PIXEL_EDGE edge){
    if (edge == PIXEL_EDGE::TOP){
        return "top";
    }
    else if (edge == PIXEL_EDGE::RIGHT){
        return "right";
    }
    else if (edge == PIXEL_EDGE::BOTTOM){
        return "bottom";
    }
    else{
        return "left";
    }
}

Ray::~Ray(){
    DEBUG(DB_TRACE, "Ray is getting deleted");
}

void Ray::deactivate(){
    m_active = false;
}

bool Ray::is_active(){
    return m_active;
}

PIXEL Ray::get_current_pixel(){
    return m_current_pixel;
}

double Ray::get_current_energy(){
    return m_current_energy;
}

void Ray::set_current_energy(double new_energy){
    m_current_energy = new_energy;
}

bool Ray::is_primary(){
    return m_primary;
}

PIXEL_EDGE Ray::get_current_edge(){
    return m_current_edge;
}

double Ray::get_current_edge_dist(){
    return m_current_edge_dist;
}