#include <random>
#include <omp.h>
#include <cmath>
#include "random_tools.hpp"

const double pi = 3.141592654;
std::default_random_engine generator[4];
std::uniform_real_distribution<double> distribution(0,1);

double roll(){
  return distribution(generator[omp_get_thread_num()]);
}

Vec2 box_muller(double sigma){
  double u1 = roll(), u2 = roll();
  return { 
    sigma * sqrt(-2 * log(u1)) * cos(2 * pi * u2),
    sigma * sqrt(-2 * log(u1)) * sin(2 * pi * u2)
  };
}
