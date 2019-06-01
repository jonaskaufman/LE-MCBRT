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
    return;
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
    return;
}

void SimulationSerial::initialize_densities_random()
{
    for (int i = 0; i < m_N; i++)
    {
        for (int j = 0; j < m_N; j++)
        {
            m_densities[i][j] = uniform_dist(random_engine);
        }
    }
    return;
}

void SimulationSerial::initialize_densities_centered_gaussian(const double max_density, const double spread)
{
    for (int i = 0; i < m_N; i++)
    {
        for (int j = 0; j < m_N; j++)
        {
            int mid = m_N / 2; // middle pixel
            int x = i - mid;
            int y = j - mid;
            double std_dev = spread * m_N;
            m_densities[i][j] = max_density * exp(-(x * x + y * y) / (2 * std_dev * std_dev));
        }
    }
    return;
}

void SimulationSerial::initialize_densities_random_gaussians(const int n_gaussians,
                                                             const double max_density,
                                                             const double spread)
{
    double std_dev = spread * m_N;
    double highest = 0;

    // Add the Gaussians
    for (int k = 0; k < n_gaussians; k++)
    {
        int mid_x = floor(uniform_dist(random_engine) * m_N);
        int mid_y = floor(uniform_dist(random_engine) * m_N);
        for (int i = 0; i < m_N; i++)
        {
            for (int j = 0; j < m_N; j++)
            {
                int x = i - mid_x;
                int y = j - mid_y;
                m_densities[i][j] += max_density * exp(-(x * x + y * y) / (2 * std_dev * std_dev));
                highest = fmax(highest, m_densities[i][j]);
            }
        }
    }

    // Normalize the resulting density distribution
    for (int i = 0; i < m_N; i++)
    {
        for (int j = 0; j < m_N; j++)
        {
            m_densities[i][j] = max_density * m_densities[i][j] / highest;
        }
    }
    return;
}

void SimulationSerial::run(int num_primary_rays)
{
    // TODO: Would it be helpful to add batching? i.e. spawn R rays, evolve them, repeat
    for (int i = 0; i < num_primary_rays; i++)
    {
        _spawn_primary_ray();
    }
    DEBUG(DB_INITPRIM, std::cout << "Spawned " << m_rays.size() << " primary rays, evolving..." << std::endl)
    _evolve_to_completion();
    return;
}

double SimulationSerial::_random_source_angle(bool normal)
{
    if (normal)
    { // Normal distribution
        double angle = normal_dist(random_engine);

        // Normalize angle
        if (angle < 0)
        {
            angle += 2 * M_PI;
        }
        else if (angle >= 2 * M_PI)
        {
            angle -= 2 * M_PI;
        }
        return angle;
    }
    else
    { // Uniform distribution between 0 and 2 pi
        return uniform_angle_dist(random_engine);
    }
}

void SimulationSerial::_spawn_primary_ray()
{
    // Randomly select source angle from normal distribution
    double source_angle = _random_source_angle(true);

    // Calculate initial ray position
    double horiz_dist_from_center = PARAM_D * tan(source_angle); // horizontal distance from center of top edge
    int middle_pixel = m_N / 2;
    double horiz_dist_from_left = middle_pixel + horiz_dist_from_center;

    // If ray does not miss grid entirely, spawn it
    if (horiz_dist_from_left < 0 || horiz_dist_from_left >= m_N ||
        (source_angle >= M_PI / 2 && source_angle <= 3 * M_PI / 2))
    {
        DEBUG(DB_INITPRIM, std::cout << "New primary ray missed the grid, not adding" << std::endl);
    }
    else
    {
        double horiz_dist_from_left_rounded = floor(horiz_dist_from_left);
        PIXEL pixel(horiz_dist_from_left_rounded, 0); // always starts from top of grid
        double edge_dist = horiz_dist_from_left - horiz_dist_from_left_rounded;
        m_rays.push_back(Ray::primary(source_angle, pixel, PIXEL_EDGE::TOP, edge_dist));
        DEBUG(DB_INITPRIM, std::cout << "New primary ray added at pixel " << pixel.first << "," << pixel.second
                                     << " with angle " << source_angle << std::endl);
    }
    return;
}

void SimulationSerial::_spawn_secondary_rays(PIXEL spawn_pixel, double total_energy)
{
    DEBUG(DB_INITSECOND, std::cout << "Spawning " << PARAM_KS << " secondary rays from pixel " << spawn_pixel.first
                                   << "," << spawn_pixel.second << "..." << std::endl);
    for (int i = 0; i < PARAM_KS; i++)
    {
        double source_angle = _random_source_angle(false); // uniform random source angle
        double partial_energy = total_energy / PARAM_KS;
        Ray new_ray = Ray::secondary_from_center(source_angle, spawn_pixel, partial_energy);
        PIXEL current_pixel = new_ray.get_current_pixel();
        if (_out_of_bounds(current_pixel))
        {
            DEBUG(DB_INITSECOND, std::cout << "New secondary ray is out of bounds, not adding" << std::endl);
        }
        else
        {
            m_rays.push_back(new_ray);
            DEBUG(DB_INITSECOND, std::cout << "New secondary ray added at pixel " << current_pixel.first << ","
                                           << current_pixel.second << " with angle " << source_angle << std::endl);
        }
    }
    return;
}

bool SimulationSerial::_out_of_bounds(const PIXEL& current_pixel)
{
    return (current_pixel.first < 0 || current_pixel.first >= m_N || current_pixel.second < 0 ||
            current_pixel.second >= m_N);
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
    }

    // Clear out the ray vector
    m_rays.clear();
    return;
}

/* Evolve all active rays.
   Return number of rays evolved
*/
int SimulationSerial::_evolve_rays()
{

    // First parallelization task: put this on GPU
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
            // new_rays = _spawn_secondary_rays(r);

            //    PIXEL current_pixel = primary->get_current_pixel();
            //    PIXEL_EDGE current_edge = primary->get_current_edge();
            //    double energy_remaining = primary->get_current_energy();
            //    double partial_energy = energy_remaining / PARAM_KS;
            //    double current_edge_dist = primary->get_current_edge_dist();
            //    double current_angle = primary->m_angle;
            DEBUG(DB_TRACE, printf("Spawning secondaries"));
            _spawn_secondary_rays(visited_pixel, m_rays[i].get_current_energy());
            m_rays[i].deactivate(); 
        }
 
        if (_out_of_bounds(current_pixel))
        {
            DEBUG(DB_TRACE, std::cout << "Ray " << i << " is out-of-bounds, deactivating" << std::endl);
            m_rays[i].deactivate();
        }

        DEBUG(
            DB_TRACE,
            printf("Ray %d finished. Visited pixel %d, %d. Distance traveled is %.2f. Generated %lu secondary rays\n\n",
                   i, visited_pixel.first, visited_pixel.second, distance_traveled, new_rays.size()));
    }
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
    double probability = std::exp(-PARAM_A / l_ep);
    double rand = uniform_dist(random_engine);

    bool result = rand < probability;
    std::string result_string = result ? "TRUE" : "FALSE";

    DEBUG(DB_INTERACT, printf("l_ep: %.2f\tprobability: %.2f\trand: %.2f\t", l_ep, probability, rand));
    DEBUG(DB_INTERACT, std::cout << "interact? " << result_string << std::endl);

    if (result)
    {
        double energy_deposited = PARAM_F * l_ep * r->get_current_energy();
        m_doses[i][j] += energy_deposited;
        r->set_current_energy(r->get_current_energy() - energy_deposited);
        DEBUG(DB_INTERACT, printf("deposited %.2f energy into pixel %d, %d\n", energy_deposited, i, j));
    }

    return result;
}

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
    return;
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
    return;
}

/* Secondary rays deposit energy into the pixel they visited
   Rethink this???? TODO: what do you mean?

   Currently, each secondary ray has initial energy of E0 / PARAM_KS.
   So deposit random fraction of this initial energy
*/
void SimulationSerial::_deposit_energy(Ray* r, PIXEL visited, double distance)
{
    int i = visited.first, j = visited.second;
    double density = m_densities[i][j];
    double l_ep = density * distance;

    double current_energy = r->get_current_energy();
    double energy_deposited =
        std::min(PARAM_G * l_ep * PARAM_E0 / PARAM_KS, current_energy); // TODO shouldnt this be current energy?
    m_doses[i][j] += energy_deposited;
    r->set_current_energy(current_energy - energy_deposited);
    if (r->get_current_energy() < PARAM_MINERGY)
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
    return;
}

