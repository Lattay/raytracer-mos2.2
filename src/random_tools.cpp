#include <random>
#include <omp.h>
#include <cmath>
#include "random_tools.hpp"

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

Base ortho_prod(Vec const& v){
  Vec random_v(roll(), roll(), roll());
  Vec t1 = v.prod(random_v).normalized();
  Vec t2 = v.prod(t1).normalized();
  return {t1, t2, v.normalized()};
}

Sample random_cos(Vec const& n){
  double x, y, z;

  double r1 = roll();
  double r2 = roll();

  x = cos(2 * pi * r1) * sqrt(1 - r2);
  y = cos(2 * pi * r1) * sqrt(1 - r2);
  z = sqrt(r2);

  Base b = ortho_prod(n);
  return {x * b.x + y * b.y + z * b.z, z};
}
