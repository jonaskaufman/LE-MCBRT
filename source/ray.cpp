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

}


std::pair<double, PIXEL> Ray::_trace()
{
    if (m_current_edge == PIXEL_EDGE::TOP || m_current_edge == PIXEL_EDGE::BOTTOM){
        return _tb_trace();
    }
    else{
        return _lr_trace();
    }


}

std::pair<double, PIXEL> Ray::_lr_trace(){
    

}

std::pair<double, PIXEL> Ray::_tb_trace(){
    double pixel_x = m_current_edge_dist - m_current_pixel.first;
    double pixel_y = m_current_pixel.second;
    DEBUG(DB_TRACE, printf("pixel x,y: %.2f, %.2f\n", pixel_x, pixel_y));

    if (m_angle < M_PI / 2){ // going SW
        // hits left wall
        DEBUG(DB_TRACE, printf("moving SW. Angle: %.2f\n", m_angle * 180 / M_PI))
        double alpha1 = M_PI / 2 - m_angle;
        double a1 = pixel_x;
        double b1 = a1 * tan(alpha1);
        double c1 = sqrt(a1 * a1 + b1 * b1);
        DEBUG(DB_TRACE, printf("a1, b1, c1: %.2f, %.2f, %.2f", a1, b1, c1));

        // hits bottom wall
        double alpha2 = m_angle;
        double b2 = 1;
        double a2 = b2 * tan(alpha2);
        double c2 = sqrt(a2 * a2 + b2 * b2);
        DEBUG(DB_TRACE, printf("a2, b2, c2: %.2f, %.2f, %.2f", a2, b2, c2));

        if (c1 < c2){ // hits left wall
            PIXEL next_pixel(m_current_pixel.first - 1, m_current_pixel.second);
            std::pair<double, PIXEL> result(c1, next_pixel);
            return result;
        } 
        else{       // hits bottom wall
            PIXEL next_pixel(m_current_pixel.first, m_current_pixel.second + 1);
            std::pair<double, PIXEL> result(c2, next_pixel);
            return result;
        }
    }
    else if (m_angle < M_PI){ // going NW
        DEBUG(DB_TRACE, printf("moving NW. Angle: %.2f\n", m_angle * 180 / M_PI))
        // hits left wall
        double alpha1 = m_angle - M_PI / 2;
        double a1 = pixel_x;
        double b1 = a1 * tan(alpha1);
        double c1 = sqrt(a1 * a1 + b1 * b1);
        DEBUG(DB_TRACE, printf("a1, b1, c1: %.2f, %.2f, %.2f", a1, b1, c1));

        // hits top wall
        double alpha2 = M_PI - m_angle;
        double b2 = 1;
        double a2 = b2 * tan(alpha2);
        double c2 = sqrt(a2 * a2 + b2 * b2);
        DEBUG(DB_TRACE, printf("a2, b2, c2: %.2f, %.2f, %.2f", a2, b2, c2));

        if (c1 < c2){ // hits left wall
            PIXEL next_pixel(m_current_pixel.first - 1, m_current_pixel.second);
            std::pair<double, PIXEL> result(c1, next_pixel);
            return result;
        } 
        else{       // hits top wall
            PIXEL next_pixel(m_current_pixel.first, m_current_pixel.second - 1);
            std::pair<double, PIXEL> result(c2, next_pixel);
            return result;
        }
    }
    else if (m_angle < 3 * M_PI / 2){ // going NE
        DEBUG(DB_TRACE, printf("moving NE. Angle: %.2f\n", m_angle * 180 / M_PI))
        // hits right wall
        double alpha1 = 3 * M_PI / 2 - m_angle;
        double a1 = pixel_x;
        double b1 = a1 * tan(alpha1);
        double c1 = sqrt(a1 * a1 + b1 * b1);
        DEBUG(DB_TRACE, printf("a1, b1, c1: %.2f, %.2f, %.2f", a1, b1, c1));

        // hits top wall
        double alpha2 = m_angle - M_PI;
        double b2 = 1;
        double a2 = b2 * tan(alpha2);
        double c2 = sqrt(a2 * a2 + b2 * b2);
        DEBUG(DB_TRACE, printf("a2, b2, c2: %.2f, %.2f, %.2f", a2, b2, c2));

        if (c1 < c2){ // hits right wall
            PIXEL next_pixel(m_current_pixel.first + 1, m_current_pixel.second);
            std::pair<double, PIXEL> result(c1, next_pixel);
            return result;
        } 
        else{       // hits top wall
            PIXEL next_pixel(m_current_pixel.first, m_current_pixel.second - 1);
            std::pair<double, PIXEL> result(c2, next_pixel);
            return result;
        }
    }
    else if (m_angle < 2 * PI){ // going SE
        DEBUG(DB_TRACE, printf("moving SE. Angle: %.2f\n", m_angle * 180 / M_PI))
        // hits right wall
        double alpha1 = m_angle - 3 * M_PI / 2;
        double a1 = pixel_x;
        double b1 = a1 * tan(alpha1);
        double c1 = sqrt(a1 * a1 + b1 * b1);
        DEBUG(DB_TRACE, printf("a1, b1, c1: %.2f, %.2f, %.2f", a1, b1, c1));

        // hits bottom wall
        double alpha2 = 2 * M_PI - m_angle;
        double b2 = 1;
        double a2 = b2 * tan(alpha2);
        double c2 = sqrt(a2 * a2 + b2 * b2);
        DEBUG(DB_TRACE, printf("a2, b2, c2: %.2f, %.2f, %.2f", a2, b2, c2));

        if (c1 < c2){ // hits right wall
            PIXEL next_pixel(m_current_pixel.first + 1, m_current_pixel.second);
            std::pair<double, PIXEL> result(c1, next_pixel);
            return result;
        } 
        else{       // hits bottom wall
            PIXEL next_pixel(m_current_pixel.first, m_current_pixel.second + 1);
            std::pair<double, PIXEL> result(c2, next_pixel);
            return result;
        }
    }
    else{
        std::cout << "unexpected angle to trace";
    }
}

