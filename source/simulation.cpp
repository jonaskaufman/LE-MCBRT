#include "simulation.hpp"

SimulationSerial::SimulationSerial(const int N) : m_N(N)
{
    for (int i = 0; i < m_N; i++)
    {
        std::vector<double> temp;
        for (int j = 0; j < m_N; j++)
        {
            temp.push_back(0);
        }
        m_densities.push_back(temp);
        m_doses.push_back(temp);
    }
}


void SimulationSerial::initialize_densities_constant(const double density)
{
    for (int i = 0; i < m_N; i++)
    {
        for (int j = 0; j < m_N; j++)
        {
            m_densities[i][j] = density;
        }
    }
}


void SimulationSerial::initialize_densities_random()
{
    for (int i = 0; i < m_N; i++)
    {
        for (int j = 0; j < m_N; j++)
        {
            m_densities[i][j] = uniform_dist(random_engine);
            DEBUG(DB_SIMCONST,
                  std::cout << "m_densities[" << i << "][" << j << "]: " << m_densities[i][j] << std::endl);
        }
    }
}


void SimulationSerial::run(int num_primary_rays)
{
    // a thought: would it help with memory issues to spawn a certain number at a time, then evolve those?
    // would require clearing the ray vector in _evolve_to_completion
    for (int i = 0; i < num_primary_rays; i++)
    {
        _spawn_primary_ray();
    }
    _evolve_to_completion();

    DEBUG(DB_INITPRIM, printf("spawned %lu/%d primary rays\n", m_rays.size(), num_primary_rays));
    return;
}


double SimulationSerial::_random_source_angle(bool normal)
{
    if (normal)
    { // normal
        return normal_dist(random_engine);
    }
    else
    { // uniform between 0 and 2 pi
        return uniform_angle_dist(random_engine);
    }
}


void SimulationSerial::_spawn_primary_ray()
{
    // Randomly select source angle, check if it is valid
    double source_angle;
    bool valid = false;
    while (!valid)
    {
        source_angle = _random_source_angle(true); // normal distribution
        // source_angle = _normalize_angle(source_angle);    // TODO: is this necessary?

        DEBUG(DB_INITPRIM, printf("Angle: %.2f\t", source_angle));

        if (source_angle >= M_PI / 4 && source_angle <= 3 * M_PI / 4)
        {
            DEBUG(DB_INITPRIM, printf("ray spawned backwards from source")); 
        }
        else
        {
            valid = true;
        }
    }    

    // Calculate initial ray position  
    double horiz_dist_from_center = D * tan(source_angle); // horizontal distance from center of top edge
    //DEBUG(DB_INITPRIM, printf("horiz_dist: %.2f\t", horiz_dist));
    int middle_pixel = m_N / 2;
    double horiz_dist_from_left = middle_pixel + horiz_dist_from_center;
     
    // If ray does not miss grid entirely, spawn it
    if (horiz_dist_from_left < 0 || horiz_dist_from_left >= m_N )
    {
        DEBUG(DB_INITPRIM, printf("primay ray missed the grid\n")); 
    }
    else
    {   
        double horiz_dist_from_left_rounded = floor(horiz_dist_from_left);
        PIXEL pixel(horiz_dist_from_left_rounded, 0); // starts from top of grid 
        double edge_dist = horiz_dist_from_left - horiz_dist_from_left_rounded;
        m_rays.push_back( Ray::primary(source_angle, pixel, PIXEL_EDGE::TOP, edge_dist) );
    }

    return;
}

void SimulationSerial::_spawn_secondary_rays(PIXEL spawn_pixel, double total_energy)
{
    for (int i = 0; i < KS; i++)
    {
        DEBUG(DB_SECONDARY, printf("Secondary %d:\n", i));
        double source_angle = _random_source_angle(false); // random source angle (normal dist=false)
        double partial_energy = total_energy / KS;
        m_rays.push_back( Ray::secondary_from_center(source_angle, spawn_pixel, partial_energy) );
    }

    return;
}

/* Continue to evolve rays as long as there are rays to evolve.
   Couldn't get vector to delete the rays so rays are deactivated instead.
   If a ray is deactivated, for loop goes to next ray
   For loop returns number of rays it evolved.
*/
void SimulationSerial::_evolve_to_completion()
{
    int rays_evolved = m_rays.size();
    while (rays_evolved > 0)
    {
        rays_evolved = _evolve_rays();
        DEBUG(DB_SECONDARY, printf("rays_evolved: %d\n\n", rays_evolved));
        DEBUG(DB_SECONDARY, std::this_thread::sleep_for(std::chrono::milliseconds(0)));
    }
    // rays.clear()
}


/* Evolve all active rays.
   Return number of rays evolved
*/
int SimulationSerial::_evolve_rays()
{

    // First parallelization taks: put this on GPU
    int rays_evolved = 0;
    for (int i = 0; i < m_rays.size(); i++)
    {
        Ray* r = &m_rays[i];

        if (r->is_active() == false)
        {
            DEBUG(DB_SECONDARY, printf("ray %d is not active\n", i));
            continue;
        }
        DEBUG(DB_SECONDARY, printf("ray %d is active\n", i));
        if (r->is_primary() == true)
        {
            DEBUG(DB_SECONDARY, printf("ray %d is primary\n", i));
        }
        else
        {
            DEBUG(DB_SECONDARY, printf("ray %d is secondary\n", i));
        }

        rays_evolved++;

        /// trace ray
        std::pair<double, PIXEL> rtrace = r->trace();
        double distance_traveled = rtrace.first;
        PIXEL visited_pixel = rtrace.second;          // pixel visited TODO spawn secondary rays from here
        PIXEL current_pixel = r->get_current_pixel(); // updated new pixel? what if it goes off the edge? TODO

        std::vector<Ray> new_rays; // TODO needed if primary?

        if (r->is_primary() == false)
        { // secondary rays deposit energy
            _deposit_energy(r, visited_pixel, distance_traveled);
            DEBUG(DB_SECONDARY, printf("Remaining energy: %.2f\n", r->get_current_energy()));
        }
        else if (_random_interact(r, visited_pixel, distance_traveled))
        { // primary rays check for interaction

            DEBUG(DB_TRACE, printf("Deactivated vector %d because of interaction\n", i));
            //new_rays = _spawn_secondary_rays(r);
            
            
//    PIXEL current_pixel = primary->get_current_pixel();
//    PIXEL_EDGE current_edge = primary->get_current_edge();
//    double energy_remaining = primary->get_current_energy();
//    double partial_energy = energy_remaining / KS;
//    double current_edge_dist = primary->get_current_edge_dist();
//    double current_angle = primary->m_angle;
            _spawn_secondary_rays(visited_pixel, m_rays[i].get_current_energy());
            m_rays[i].deactivate(); // deactivate primary ray. Not sure why r->deactivate wouldn't work TODO check this
        }

        // TODO feel like this should be handled in ray, but I guess that class has no knowledge of simulation size
        // seems fine then
        if (current_pixel.first < 0 || current_pixel.first >= m_N || current_pixel.second < 0 ||
            current_pixel.second >= m_N)
        {
            m_rays[i].deactivate();
            DEBUG(DB_TRACE, printf("Deactivated vector %d because out-of-bounds\n", i));
        }

        DEBUG(
            DB_TRACE,
            printf("Ray %d finished. Visited pixel %d, %d. Distance traveled is %.2f. Generated %lu secondary rays\n\n",
                   i, visited_pixel.first, visited_pixel.second, distance_traveled, new_rays.size()));

        // DEBUG(DB_SECONDARY, std::this_thread::sleep_for(std::chrono::milliseconds(1000)));
    }
    // DEBUG(DB_TRACE, printf("%lu vectors remaining\n", m_rays.size()));
    return rays_evolved;
}


/* Determine whether primary ray interacts at current pixel
   If it does, deposit fraction of energy there
   */
bool SimulationSerial::_random_interact(Ray* r, PIXEL visited, double distance)
{
    int i = visited.first, j = visited.second;
    double density = m_densities[i][j];
    double l_ep = density * distance;
    double probability = std::exp(-A / l_ep);
    double rand = uniform_dist(random_engine);

    bool result = rand < probability;
    std::string result_string = result ? "TRUE" : "FALSE";

    DEBUG(DB_INTERACT, printf("l_ep: %.2f\tprobability: %.2f\trand: %.2f\t", l_ep, probability, rand));
    DEBUG(DB_INTERACT, std::cout << "interact? " << result_string << std::endl);

    if (result)
    {
        double energy_deposited = F * l_ep * r->get_current_energy();
        m_doses[i][j] += energy_deposited;
        r->set_current_energy(r->get_current_energy() - energy_deposited);
        DEBUG(DB_INTERACT, printf("deposited %.2f energy into pixel %d, %d\n", energy_deposited, i, j));
    }

    return result;
}


// Also just take a pixel and a partial energy
/// Generate secondary rays from primary ray
/*
std::vector<Ray> SimulationSerial::_spawn_secondary_rays(Ray* primary)
{
    std::vector<Ray> result;

    PIXEL current_pixel = primary->get_current_pixel();
    PIXEL_EDGE current_edge = primary->get_current_edge();
    double energy_remaining = primary->get_current_energy();
    double partial_energy = energy_remaining / KS;
    double current_edge_dist = primary->get_current_edge_dist();
    double current_angle = primary->m_angle;

    for (int i = 0; i < KS; i++)
    {
        DEBUG(DB_SECONDARY, printf("Secondary %d:\n", i));
        double source_angle = _random_source_angle(false); // random source angle (normal dist=false)

        std::pair<int, int> adjustments = _fix_position(current_edge, current_angle, source_angle);

        // TODO: make adjust pixel function?
        // More generally why does secondary ray origin depend on primary ray angle?
        // I say we just spawn secondary rays from center of pixel where interaction occured
        PIXEL adjusted_pixel(current_pixel.first + adjustments.first, current_pixel.second + adjustments.second);

        Ray secondary_ray = Ray(false, source_angle, adjusted_pixel, current_edge, current_edge_dist, partial_energy);

        // secondary_ray.fix_position()

        if (adjustments.first != 0 || adjustments.second != 0)
        {
            secondary_ray.set_corrected();
        }

        // why are both of these being updated?
        result.push_back(secondary_ray);
        m_rays.push_back(secondary_ray);
    }

    DEBUG(DB_SECONDARY, printf("spawned %d secondary rays\n", KS));
    return result;
}
*/

void SimulationSerial::print_densities()
{
    std::cout << "m_densities: " << std::endl;
    for (int i = 0; i < m_N; i++)
    {
        for (int j = 0; j < m_N; j++)
        {
            printf("%.2f ", m_densities[i][j]);
        }
        std::cout << std::endl;
    }
}


void SimulationSerial::print_doses()
{
    std::cout << "m_doses: " << std::endl;
    for (int i = 0; i < m_N; i++)
    {
        for (int j = 0; j < m_N; j++)
        {
            printf("%.2f ", m_doses[i][j]);
        }
        std::cout << std::endl;
    }
}


double SimulationSerial::_normalize_angle(double angle)
{
    while (angle < 0)
    {
        angle += 2 * M_PI;
    }
    while (angle >= 2 * M_PI)
    {
        angle -= 2 * M_PI;
    }
    return angle;
}


/* Secondary rays deposit energy into the pixel they visited
   Rethink this???? TODO: what do you mean?

   Currently, each secondary ray has initial energy of E0 / KS.
   So deposit random fraction of this initial energy
*/
void SimulationSerial::_deposit_energy(Ray* r, PIXEL visited, double distance)
{
    int i = visited.first, j = visited.second;
    double density = m_densities[i][j];
    double l_ep = density * distance;

    double current_energy = r->get_current_energy();
    double energy_deposited = std::min(G * l_ep * E0 / KS, current_energy);
    m_doses[i][j] += energy_deposited;
    r->set_current_energy(current_energy - energy_deposited);
    if (r->get_current_energy() < MIN_ENERGY)
    {
        r->deactivate();
        DEBUG(DB_SECONDARY, printf("secondary ray ran out of energy.\n"));
    }
    DEBUG(DB_SECONDARY, printf("deposited %.2f energy into pixel %d, %d\n", energy_deposited, i, j));
}

/*
PIXEL SimulationSerial::_fix_position(PIXEL_EDGE edge, double current_angle, double new_angle)
{
    std::pair<int, int> result(0, 0);
    // DEBUG(DB_SECONDARY, printf("fixing position. Current angle: %.2f\t New angle: %.2f\n", current_angle * 180 /
    // M_PI, new_angle * 180 / M_PI));
    if (edge == PIXEL_EDGE::RIGHT)
    {
        if (current_angle > M_PI && new_angle < M_PI)
        {
            result.first = -1;
            DEBUG(DB_SECONDARY,
                  printf("Spawned secondary going left but primary was going right. Adjusting pixel location\n"));
        }
    }
    else if (edge == PIXEL_EDGE::LEFT)
    {
        if (current_angle < M_PI && new_angle > M_PI)
        {
            result.first = 1;
            DEBUG(DB_SECONDARY,
                  printf("Spawned secondary going right but primary was going left. Adjusting pixel location\n"));
        }
    }
    else if (edge == PIXEL_EDGE::TOP)
    {
        if ((current_angle > M_PI / 2 && current_angle < 3 * M_PI / 2) &&
            (new_angle < M_PI / 2 || new_angle > 3 * M_PI / 2))
        {
            result.second = 1;
            DEBUG(DB_SECONDARY,
                  printf("Spawned secondary going down but primary was going up. Adjusting pixel location\n"));
        }
    }
    else
    {
        if ((current_angle < M_PI / 2 || current_angle > 3 * M_PI / 2) &&
            (new_angle > M_PI / 2 && new_angle < 3 * M_PI / 2))
        {
            result.second = -1;
            DEBUG(DB_SECONDARY,
                  printf("Spawned secondary going up but primary was going down. Adjusting pixel location\n"));
        }
    }

    return result;
}
*/

void SimulationSerial::write_to_file()
{
    std::ofstream densities;
    densities.open("densities.csv");
    for (int j = 0; j < m_densities.size(); j++)
    {
        for (int i = 0; i < m_densities.size() - 1; i++)
        {
            densities << m_densities[i][j] << ",";
        }
        densities << m_densities[m_densities.size() - 1][j] << "\n";
    }
    densities.close();

    std::ofstream doses;
    doses.open("doses.csv");
    for (int j = 0; j < m_doses.size(); j++)
    {
        for (int i = 0; i < m_doses.size() - 1; i++)
        {
            doses << m_doses[i][j] << ",";
        }
        doses << m_doses[m_doses.size() - 1][j] << "\n";
    }
    doses.close();
}

