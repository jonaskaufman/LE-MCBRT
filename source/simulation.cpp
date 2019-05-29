#include "simulation.hpp"

SimulationSerial::SimulationSerial(const int N, const double density, const int ray_count) : N(N)
{

   for (int i = 0; i < N; i++){
      std::vector<double> temp;
      for (int j = 0; j < N; j++){
         temp.push_back(0);
      }
      m_densities.push_back(temp);
      m_doses.push_back(temp);
   }
   

   if (density == -1) {
      DEBUG(DB_SIMCONST, std::cout << "density is -1" << std::endl);
      initialize_densities_random();
   }
   else{
      DEBUG(DB_SIMCONST, std::cout << "density is not -1" << std::endl);
      initialize_densities_constant(density);
   }
   
   for (int i = 0; i < ray_count; i++){
      _spawn_primary_ray();
   }
   DEBUG(DB_INITPRIM, printf("spawned %lu/%d primary rays\n", m_rays.size(), ray_count));
   
   _evolve_to_completion();
   //_evolve_rays(); // test one step of evolution
}

void SimulationSerial::initialize_densities_constant(const double density)
{
   for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         m_densities[i][j] = density;
      }
   }
}

void SimulationSerial::initialize_densities_random()
{
   /*
   std::random_device rd;
   std::mt19937 mt(rd());
   std::uniform_real_distribution<double> uni_dist(0.0, 1.0);
   */
   for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         m_densities[i][j] = uniform_dist(random_engine);
         DEBUG(DB_SIMCONST, std::cout << "m_densities[" << i << "][" << j << "]: " << m_densities[i][j] << std::endl);
      }
   }
}

double SimulationSerial::_random_source_angle(bool normal){
   if (normal){ // normal
      return normal_dist(random_engine);
   }
   else{ // uniform between 0 and 2 pi
      return uniform_angle_dist(random_engine);
   }
}


void SimulationSerial::_spawn_primary_ray(){
   double source_angle = _random_source_angle(true); // random source angle (normal dist)
   source_angle = _normalize_angle(source_angle);

   DEBUG(DB_INITPRIM, printf("Angle: %.2f\t", source_angle));

   double horiz_dist = D * tan(source_angle); // horizontal distance from center of top edge
   DEBUG(DB_INITPRIM, printf("horiz_dist: %.2f\t", horiz_dist));

   int middle_pixel = N / 2;
   
   double x_coord = middle_pixel + horiz_dist; // x coordinate on top edge that primary ray originates from
   
   if (x_coord < 0 || x_coord >= N){ // primary ray missed the grid
      DEBUG(DB_INITPRIM, printf("\n"));
      return;
   }
   DEBUG(DB_INITPRIM, printf("x_coord: %.2f\n", x_coord));

   PIXEL current_pixel (floor(x_coord), 0); // current pixel rounds down x_coord and starts at top
   
   

   Ray ray = Ray(true, source_angle, current_pixel, PIXEL_EDGE::TOP, x_coord, E0); // passing E0 energy but primary won't deposit it
   m_rays.push_back(ray);
   
}

/* Continue to evolve rays as long as there are rays to evolve.
   Couldn't get vector to delete the rays so rays are deactivated instead.
   If a ray is deactivated, for loop goes to next ray
   For loop returns number of rays it evolved.
*/
void SimulationSerial::_evolve_to_completion(){
   int rays_evolved = m_rays.size();
   while (rays_evolved > 0){
      rays_evolved = _evolve_rays();
      DEBUG(DB_SECONDARY, printf("rays_evolved: %d\n\n", rays_evolved));
      DEBUG(DB_SECONDARY, std::this_thread::sleep_for(std::chrono::milliseconds(10)));
   }
}
/* Evolve all active rays.
   Return number of rays evolved
*/
int SimulationSerial::_evolve_rays(){
   std::vector<int> to_remove;
   int rays_evolved = 0;
   for (int i = 0; i < m_rays.size(); i++){
      Ray *r = &m_rays[i];
      
      if (r->is_active() == false){
         DEBUG(DB_SECONDARY, printf("ray %d is not active\n", i));
         continue;
      }
      DEBUG(DB_SECONDARY, printf("ray %d is active\n", i));
      if (r->is_primary() == true){
         DEBUG(DB_SECONDARY, printf("ray %d is primary\n", i));
      }
      else{
         DEBUG(DB_SECONDARY, printf("ray %d is secondary\n", i));
      }

      rays_evolved++;

      std::pair<double, PIXEL> r_trace = r->_trace();
      double distance_traveled = r_trace.first;
      PIXEL visited_pixel = r_trace.second;
      PIXEL current_pixel = r->get_current_pixel();

      std::vector<Ray> new_rays;

      if (r->is_primary() == false){ // secondary rays deposit energy
         _deposit_energy(r, visited_pixel, distance_traveled);
         DEBUG(DB_SECONDARY, printf("Remaining energy: %.2f\n", r->get_current_energy()));
      }
      else if (_random_interact(r, visited_pixel, distance_traveled)){ // primary rays check for interaction
         //to_remove.push_back(i);
         //m_rays.erase(m_rays.begin() + i);
         
         DEBUG(DB_TRACE, printf("Deactivated vector %d because of interaction\n", i));
         new_rays = _spawn_secondary_rays(r);
         m_rays[i].deactivate(); // deactivate primary ray. Not sure why r->deactivate wouldn't work
         
      }
      
      if (current_pixel.first < 0 || current_pixel.first >= N || \
          current_pixel.second < 0 || current_pixel.second >= N){
         m_rays[i].deactivate();
         DEBUG(DB_TRACE, printf("Deactivated vector %d because out-of-bounds\n", i));
      }

      
      
      DEBUG(DB_TRACE, printf("Ray %d finished. Visited pixel %d, %d. Distance traveled is %.2f. Generated %lu secondary rays\n\n", \
         i, visited_pixel.first, visited_pixel.second, distance_traveled, new_rays.size()));

      //DEBUG(DB_SECONDARY, std::this_thread::sleep_for(std::chrono::milliseconds(1000)));
      
   }
   //DEBUG(DB_TRACE, printf("%lu vectors remaining\n", m_rays.size()));
   return rays_evolved;
   /*
   for (int i = 0; i < to_remove.size(); i++){
      int remove = to_remove[i] - i;
      m_rays.erase(m_rays.begin() + remove);
   }*/
}


/* Determine whether primary ray interacts at current pixel
   If it does, deposit fraction of energy there
   */
bool SimulationSerial::_random_interact(Ray *r, PIXEL visited, double distance){
   int i = visited.first, j = visited.second;
   double density = m_densities[i][j];
   double l_ep = density * distance;
   double probability = std::exp(-A / l_ep);
   double rand = uniform_dist(random_engine);

   bool result = rand < probability;
   std::string result_string = result ? "TRUE" : "FALSE";

   DEBUG(DB_INTERACT, printf("l_ep: %.2f\tprobability: %.2f\trand: %.2f\t", l_ep, probability, rand));
   DEBUG(DB_INTERACT, std::cout << "interact? " << result_string << std::endl);

   if (result){
      double energy_deposited = F * l_ep * r->get_current_energy();
      m_doses[i][j] += energy_deposited;
      r->set_current_energy(r->get_current_energy() - energy_deposited);
      DEBUG(DB_INTERACT, printf("deposited %.2f energy into pixel %d, %d\n", energy_deposited, i, j));
   }
   
   return result;
}

/// Generate secondary rays from primary ray
std::vector<Ray> SimulationSerial::_spawn_secondary_rays(Ray *primary){
   std::vector<Ray> result;
   

   PIXEL current_pixel = primary->get_current_pixel();
   PIXEL_EDGE current_edge = primary->get_current_edge();
   double energy_remaining = primary->get_current_energy();
   double partial_energy = energy_remaining / KS;
   double current_edge_dist = primary->get_current_edge_dist();

   for (int i = 0; i < KS; i++){
      DEBUG(DB_SECONDARY, printf("Secondary %d:\n", i));
      double source_angle = _random_source_angle(false); // random source angle (normal dist=false)
      std::pair<int, int> adjustments = _fix_position(current_edge, primary->m_angle, source_angle);
      
      PIXEL adjusted_pixel(current_pixel.first + adjustments.first, current_pixel.second + adjustments.second);
   
      Ray secondary_ray = Ray(false, source_angle, adjusted_pixel, current_edge, current_edge_dist, partial_energy);
      if (adjustments.first != 0 || adjustments.second != 0){
         secondary_ray.set_corrected();
      }
      result.push_back(secondary_ray);
      m_rays.push_back(secondary_ray);
   }


   DEBUG(DB_SECONDARY, printf("spawned %d secondary rays\n", KS));
   return result;
}

void SimulationSerial::print_m_densities(){
   std::cout << "m_densities: " << std::endl;
   for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         printf("%.2f ", m_densities[i][j]);
      }
      std::cout << std::endl;
   }
}

void SimulationSerial::print_m_doses(){
   std::cout << "m_doses: " << std::endl;
   for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         printf("%.2f ", m_doses[i][j]);
      }
      std::cout << std::endl;
   }
}

double SimulationSerial::_normalize_angle(double angle){
   while (angle < 0){
      angle += 2 * M_PI;
   }
   while (angle >= 2 * M_PI){
      angle -= 2 * M_PI;
   }
   return angle;
}

/* Secondary rays deposit energy into the pixel they visited
   Rethink this????
   
   Currently, each secondary ray has initial energy of E0 / KS.
   So deposit random fraction of this initial energy
*/
void SimulationSerial::_deposit_energy(Ray *r, PIXEL visited, double distance){
   int i = visited.first, j = visited.second;
   double density = m_densities[i][j];
   double l_ep = density * distance;
   double initial_energy = E0 / KS;
   double current_energy = r->get_current_energy();
   double energy_deposited = std::min(G * l_ep * E0 / KS, current_energy);
   m_doses[i][j] += energy_deposited;
   r->set_current_energy(current_energy - energy_deposited);
   if (r->get_current_energy() < MIN_ENERGY){
      r->deactivate();
      DEBUG(DB_SECONDARY, printf("secondary ray ran out of energy.\n"));
   }
   DEBUG(DB_SECONDARY, printf("deposited %.2f energy into pixel %d, %d\n", energy_deposited, i, j));
}

/* Fixes position discrepency when spawning secondary rays that are going in the opposite direction of the primary ray
   Returns corrections to current pixel
*/
std::pair<int,int> SimulationSerial::_fix_position(PIXEL_EDGE edge, double current_angle, double new_angle){
   std::pair<int,int> result(0,0);
   if (edge == PIXEL_EDGE::RIGHT){
      if (current_angle > M_PI && new_angle < M_PI){
         result.first = -1;
         DEBUG(DB_SECONDARY, printf("Spawned secondary going left but primary was going right. Adjusting pixel location\n"));
      }
   }
   else if (edge == PIXEL_EDGE::LEFT){
      if (current_angle < M_PI && new_angle > M_PI){
         result.first = 1;
         DEBUG(DB_SECONDARY, printf("Spawned secondary going right but primary was going left. Adjusting pixel location\n"));
      }
   }
   else if (edge == PIXEL_EDGE::TOP){
      if ( (current_angle > M_PI / 2 && current_angle < 3 * M_PI / 2) && \
           (new_angle < M_PI / 2 || new_angle > 3 * M_PI / 2) ){
         result.second = 1;
         DEBUG(DB_SECONDARY, printf("Spawned secondary going down but primary was going up. Adjusting pixel location\n"));
      }
   }
   else{
      if ( (current_angle < M_PI / 2 || current_angle > 3 * M_PI / 2) && \
           (new_angle > M_PI / 2 && new_angle < 3 * M_PI / 2) ){
         result.second = -1;
         DEBUG(DB_SECONDARY, printf("Spawned secondary going up but primary was going down. Adjusting pixel location\n"));
      }
   }

   return result;
}

void SimulationSerial::write_to_file(){
   std::ofstream densities;
   densities.open("densities.csv");
   for (int i = 0; i < m_densities.size(); i++){
      for (int j = 0; j < m_densities.size() - 1; j++){
         densities << m_densities[i][j] << ",";
      }
      densities << m_densities[i][m_densities.size() - 1] << "\n";
   }
   densities.close();

   std::ofstream doses;
   doses.open("doses.csv");
   for (int i = 0; i < m_doses.size(); i++){
      for (int j = 0; j < m_doses.size() - 1; j++){
         doses << m_doses[i][j] << ",";
      }
      doses << m_doses[i][m_doses.size() - 1] << "\n";
   }
   doses.close();


}
