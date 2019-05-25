#include "simulation.hpp"


SimulationSerial::SimulationSerial(const int N, const double density) : N(N)
{
   m_doses = new double *[N];
   m_densities = new double *[N];
   

   if (density == -1) {
      initialize_densities_random();
   }
   else{
      initialize_densities_constant(density);
   }

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
   for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         std::random_device rd;
         std::mt19937 mt(rd());
         std::uniform_real_distribution<double> uni_dist(0.0, 1.0);
         m_densities[i][j] = uni_dist(mt);
      }
   }
}

void SimulationSerial::print_m_densities(){
   for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         std::cout << m_densities[i][j] << " ";
      }
      std::cout << std::endl;
   }
}

void SimulationSerial::print_m_doses(){
   for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         std::cout << m_doses[i][j] << " ";
      }
      std::cout << std::endl;
   }
}