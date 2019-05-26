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

std::tuple<PIXEL, double, std::vector<Ray>> Ray::evolve(){
    std::pair<double, PIXEL> trace = _trace();
    double distance_traveled = trace.first;
    PIXEL next_pixel = trace.second;
    DEBUG(DB_TRACE, printf("Going from %d, %d to %d, %d\n", \
        m_current_pixel.first, m_current_pixel.second, next_pixel.first, next_pixel.second)); 
    std::tuple<PIXEL, double, std::vector<Ray>> result;

    std::vector<Ray> new_rays;
    result = std::make_tuple(next_pixel, distance_traveled, new_rays);
    return result;
}


std::pair<double, PIXEL> Ray::_trace()
{
    double pixel_x, pixel_y;
    double a1, b2; // a is horizontal triangle legs, b is verticle triangle legs

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

    double b1 = a1 * tan(alpha1);
    double c1 = sqrt(a1 * a1 + b1 * b1);
    DEBUG(DB_TRACE, printf("a1, b1, c1: %.2f, %.2f, %.2f\n", a1, b1, c1));

    double a2 = b2 * tan(alpha2);
    double c2 = sqrt(a2 * a2 + b2 * b2);
    DEBUG(DB_TRACE, printf("a2, b2, c2: %.2f, %.2f, %.2f\n", a2, b2, c2));

    

    if (c1 < c2){
        int next_x = m_current_pixel.first + delta_x1;
        int next_y = m_current_pixel.second + delta_y1;
        PIXEL next_pixel(next_x, next_y);
        m_current_pixel = next_pixel;
        m_current_edge_dist = 0; // update edge
        m_current_edge = PIXEL_EDGE::TOP; // TO DO

        std::pair<double, PIXEL> result(c1, next_pixel);
        return result;
    } 
    else{ 
        int next_x = m_current_pixel.first + delta_x2;
        int next_y = m_current_pixel.second + delta_y2;
        PIXEL next_pixel(m_current_pixel.first + delta_x2, m_current_pixel.second + delta_y2);
        m_current_pixel = next_pixel;
        m_current_edge_dist = 0; // update edge
        m_current_edge = PIXEL_EDGE::TOP; // TO DO

        std::pair<double, PIXEL> result(c2, next_pixel);
        return result;
    }
}

