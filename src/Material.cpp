#include <omp.h>
#include <random>
#include "Material.hpp"

std::default_random_engine generator[4];
std::uniform_real_distribution<double> distribution(0,1);

const double pi = 3.14159265;

static double roll(){
  return distribution(generator[omp_get_thread_num()]);
}

Vec Reflective::reflex_dir(Vec const& source, Vec const& n) const{
  return source - 2.0 * n.dot(source) * n;
}

Vec Diffuse::reflex_dir(Vec const& source, Vec const& n) const{
  double x, y, z;

  double r1 = roll();
  double r2 = roll();

  x = cos(2 * pi * r1) * sqrt(1 - r2);
  y = cos(2 * pi * r1) * sqrt(1 - r2);
  z = sqrt(r2);

  Vec random_v(roll(), roll(), roll());
  Vec t1 = n.prod(random_v).normalized();
  Vec t2 = n.prod(t1).normalized();

  return x * t1 + y * t2 + z * n;
}

Vec Transparent::reflex_dir(Vec const& source, Vec const& n) const{
  return source;
}
