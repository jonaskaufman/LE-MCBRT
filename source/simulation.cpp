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
   

   //_evolve_rays();
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

double SimulationSerial::_random_source_angle(){
   return normal_dist(random_engine);
}


void SimulationSerial::_spawn_primary_ray(){
   double source_angle = _random_source_angle(); // random source angle (normal dist)
   while (source_angle < 0){
      source_angle += 2 * M_PI;
   }
   while (source_angle >= 2 * M_PI){
      source_angle -= 2 * M_PI;
   }
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

void SimulationSerial::_evolve_to_completion(){
   while (m_rays.size() > 0){
      _evolve_rays();
   }
}

void SimulationSerial::_evolve_rays(){
   for (Ray r: m_rays){
      std::tuple<PIXEL, double, std::vector<Ray>> r_tuple = r.evolve();
      //r.evolve();
   }
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