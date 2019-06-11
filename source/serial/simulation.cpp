#include "simulation.hpp"

Simulation::Simulation(int N) : m_N(N)
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

void Simulation::initialize_densities_random()
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

void Simulation::initialize_densities_constant(double density)
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

void Simulation::initialize_densities_centered_gaussian(double max_density, double spread)
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

void Simulation::initialize_densities_random_gaussians(int n_gaussians,
                                                       double max_density,
                                                       double spread)
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

void Simulation::run_serial(int num_primary_rays, int batch_size)
{
    int num_batches = 1 + num_primary_rays / batch_size;
    int rays_added = 0;
    for (int b = 0; b < num_batches; b++)
    {
        DEBUG(DB_INIT_PRI, std::cout << "Spawning " << num_primary_rays << " primary rays..." << std::endl);
        for (int i = 0; i < batch_size; i++)
        {
            if (rays_added < num_primary_rays)
            {
                _spawn_primary_ray();
                rays_added++;
            }
        }
        DEBUG(DB_INIT_PRI, std::cout << "Done, " << m_rays.size() << " rays added" << std::endl << std::endl);
        DEBUG(DB_GENERAL, std::cout << "Evolving rays..." << std::endl)
        _evolve_to_completion();
        DEBUG(DB_GENERAL, std::cout << "Done" << std::endl);
    }
    return;
}

void Simulation::write_densities_to_file(const std::string& filename)
{
    std::ofstream output;
    output.open(filename);
    for (int j = 0; j < m_densities.size(); j++)
    {
        for (int i = 0; i < m_densities.size() - 1; i++)
        {
            output << m_densities[i][j] << ",";
        }
        output << m_densities[m_densities.size() - 1][j] << "\n";
    }
    output.close();
    return;
}

void Simulation::write_doses_to_file(const std::string& filename)
{
    std::ofstream output;
    output.open(filename);
    for (int j = 0; j < m_doses.size(); j++)
    {
        for (int i = 0; i < m_doses.size() - 1; i++)
        {
            output << m_doses[i][j] << ",";
        }
        output << m_doses[m_doses.size() - 1][j] << "\n";
    }
    output.close();
    return;
}

double Simulation::_random_source_angle(bool normal)
{
    if (normal)
    { // Normal distribution
        double angle = normal_dist(random_engine);

        // Normalize angle
        while (angle < 0.0)
        {
            angle += 2 * M_PI;
        }
        while (angle >= 2 * M_PI)
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

bool Simulation::_out_of_bounds(const PIXEL& current_pixel)
{
    return (current_pixel.first < 0 || current_pixel.first >= m_N || current_pixel.second < 0 ||
            current_pixel.second >= m_N);
}

void Simulation::_spawn_primary_ray()
{
    // Randomly select source angle from normal distribution
    double source_angle = _random_source_angle(true);

    // Calculate initial ray position
    double horiz_dist_from_center = PARAM_D * m_N * tan(source_angle); // horizontal distance from center of top edge
    int middle_pixel = m_N / 2;
    double horiz_dist_from_left = middle_pixel + horiz_dist_from_center;

    // If ray does not miss grid entirely, spawn it
    if (horiz_dist_from_left < 0 || horiz_dist_from_left >= m_N ||
        (source_angle >= M_PI / 2 && source_angle <= 3 * M_PI / 2))
    {
        DEBUG(DB_INIT_PRI, std::cout << "New primary ray missed the grid, not adding" << std::endl);
    }
    else
    {
        double horiz_dist_from_left_rounded = floor(horiz_dist_from_left);
        PIXEL pixel(horiz_dist_from_left_rounded, 0); // always starts from top of grid
        double edge_dist = horiz_dist_from_left - horiz_dist_from_left_rounded;
        m_rays.push_back(Ray::primary(source_angle, pixel, PIXEL_EDGE::TOP, edge_dist));
        DEBUG(DB_INIT_PRI, std::cout << "New primary ray added at pixel " << pixel.first << "," << pixel.second
                                     << " with angle " << source_angle << std::endl);
    }
    return;
}

void Simulation::_spawn_secondary_rays(PIXEL spawn_pixel, double total_energy)
{
    DEBUG(DB_INIT_SEC, std::cout << "Spawning " << PARAM_KS << " secondary rays from pixel " << spawn_pixel.first << ","
                                 << spawn_pixel.second << "..." << std::endl);
    for (int i = 0; i < PARAM_KS; i++)
    {
        double source_angle = _random_source_angle(false); // uniform random source angle
        double partial_energy = total_energy / PARAM_KS;
        Ray new_ray = Ray::secondary_from_center(source_angle, spawn_pixel, partial_energy);
        PIXEL current_pixel = new_ray.get_current_pixel();
        if (_out_of_bounds(current_pixel))
        {
            DEBUG(DB_INIT_SEC, std::cout << "Ray is out of bounds, not adding" << std::endl);
        }
        else
        {
            m_rays.push_back(new_ray);
            DEBUG(DB_INIT_SEC, std::cout << "Ray is in bounds, adding" << std::endl);
        }
    }
    DEBUG(DB_INIT_SEC, std::cout << "Done" << std::endl << std::endl);
    return;
}

bool Simulation::_random_interact(PIXEL target_pixel, double distance)
{
    int i = target_pixel.first, j = target_pixel.second;
    double density = m_densities[i][j];
    double l_ep = density * distance; // effective path length travelled in pixel
    double probability = std::exp(-PARAM_A / l_ep);
    double rand = uniform_dist(random_engine);
    return (rand < probability);
}

void Simulation::_transfer_energy(Ray* ray, PIXEL target_pixel, double unscaled_energy)
{
    int i = target_pixel.first, j = target_pixel.second;
    double density = m_densities[i][j];
    double transfer_energy = unscaled_energy * density; // scale energy by pixel density
    double current_ray_energy = ray->get_current_energy();

    // Ray cannot transfer more energy that it has
    transfer_energy = fmin(transfer_energy, current_ray_energy);

    // Remove energy from ray and add it to pixel dose
    ray->set_current_energy(current_ray_energy - transfer_energy);
    m_doses[i][j] += transfer_energy;

    return;
}

int Simulation::_evolve_rays()
{

    // First parallelization task: put this on GPU
    int rays_evolved = 0;
    for (int i = 0; i < m_rays.size(); i++)
    {
        Ray* r = &m_rays[i];

        // Only evolve active rays
        if (r->is_active())
        {
            // Trace ray
            DEBUG(DB_TRACE, std::cout << "Tracing ray " << i << std::endl);
            std::pair<PIXEL, double> rtrace = r->trace();
            PIXEL visited_pixel = rtrace.first;
            double travel_distance = rtrace.second; // distance traveled in visited pixel
            rays_evolved++;

            if (r->is_primary()) // primary ray
            {
                DEBUG(DB_EVOLVE_PRI, std::cout << "Primary ray " << i << "  visited pixel " << visited_pixel.first
                                               << "," << visited_pixel.second << " with travel dist " << travel_distance
                                               << std::endl);
                if (_random_interact(visited_pixel, travel_distance))
                {
                    DEBUG(DB_EVOLVE_PRI, std::cout << "Primary ray " << i << " interacted" << std::endl);
                    // Deposit energy to pixel
                    DEBUG(DB_EVOLVE_PRI, std::cout << "Depositing energy to pixel" << std::endl);
                    DEBUG(DB_EVOLVE_PRI, std::cout << "Starting energy " << r->get_current_energy() << std::endl);
                    double energy_to_deposit = PARAM_F * travel_distance * r->get_current_energy();
                    _transfer_energy(r, visited_pixel, energy_to_deposit);
                    DEBUG(DB_EVOLVE_PRI, std::cout << "Energy after deposit " << r->get_current_energy() << std::endl);
                    DEBUG(DB_EVOLVE_PRI, std::cout << "Will spawn secondary rays next" << std::endl << std::endl);
                    // Spawn secondary rays, transferring remaining energy to them
                    _spawn_secondary_rays(visited_pixel, r->get_current_energy());
                    r->set_current_energy(0);
                }
                else
                {
                    DEBUG(DB_EVOLVE_PRI, std::cout << "No interaction" << std::endl << std::endl);
                }
            }
            else // secondary ray
            {
                DEBUG(DB_EVOLVE_SEC, std::cout << "Secondary ray " << i << " visited pixel " << visited_pixel.first
                                               << "," << visited_pixel.second << " with travel dist " << travel_distance
                                               << std::endl);
                double energy_to_deposit = PARAM_G * travel_distance;
                DEBUG(DB_EVOLVE_SEC, std::cout << "Depositing energy to pixel" << std::endl);
                DEBUG(DB_EVOLVE_SEC, std::cout << "Starting energy " << r->get_current_energy() << std::endl);
                DEBUG(DB_EVOLVE_SEC, std::cout << "Unscaled energy to deposit " << energy_to_deposit << std::endl);
                _transfer_energy(r, visited_pixel, energy_to_deposit);
                DEBUG(DB_EVOLVE_SEC, std::cout << "Energy after deposit " << r->get_current_energy() << std::endl
                                               << std::endl);
            }

            // Deactivate ray if out of energy or outside of the grid bounds
            if (r->get_current_energy() < PARAM_EPSILON || _out_of_bounds(r->get_current_pixel()))
            {
                DEBUG(DB_EVOLVE_SEC, std::cout << "Ray " << i << " is out of energy or bounds, deactivating"
                                               << std::endl
                                               << std::endl);
                r->deactivate();
            }
        }
    }
    return rays_evolved;
}

void Simulation::_evolve_to_completion()
{
    int rays_evolved = m_rays.size();
    while (rays_evolved > 0)
    {
        int prev_num_rays = m_rays.size();
        rays_evolved = _evolve_rays();
        DEBUG(DB_GENERAL, std::cout << rays_evolved << " rays evolved" << std::endl);
        DEBUG(DB_GENERAL, std::cout << (m_rays.size() - prev_num_rays) << " rays added" << std::endl << std::endl)
    }

    // Clear out the ray vector
    m_rays.clear();
    return;
}

