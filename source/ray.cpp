#include "ray.hpp"

    Ray::Ray(const bool primary, const double angle, std::pair<int, int> current_pixel, \
            SimulationSerial::PIXEL_EDGE current_edge, double current_edge_dist) : m_primary(primary), m_angle(angle)
    {
        
    }




    std::pair<double, std::pair<int, int>> Ray::_trace()
    {
        /// figure out two candiate edges



    }
