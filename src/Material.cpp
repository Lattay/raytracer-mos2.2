#include "Material.hpp"
#include "random_tools.hpp"

typedef struct {Vec x; Vec y; Vec z;} Base;

const double pi = 3.14159265;

Vec Reflective::reflex_dir(Vec const& source, Vec const& n) const{
  return source - 2.0 * n.dot(source) * n;
}

Base ortho_prod(Vec const& v){
  Vec random_v(roll(), roll(), roll());
  Vec t1 = v.prod(random_v).normalized();
  Vec t2 = v.prod(t1).normalized();
  Base b{t1, t2, v.normalized()};
  return b;
}

Base ortho_min(Vec const& v){
  double x = v.x(), y = v.y(), z = v.z();
  Vec t1, t2;

  if(std::abs(x) > std::abs(y) && std::abs(x) > std::abs(z)){
    t1 = Vec(z, y, -x).normalized();
    t2 = Vec(y, -x, z).normalized();
  } else if(std::abs(y) > std::abs(x) && std::abs(y) > std::abs(z)){
    t1 = Vec(x, z, -y).normalized();
    t2 = Vec(-y, x, z).normalized();
  } else {
    t1 = Vec(-z, y, x).normalized();
    t2 = Vec(x, -z, y).normalized();
  }

  Base b{t1, t2, v.normalized()};
  return b;
}

Base ortho_mixt(Vec v){
  double x = v.x(), y = v.y(), z = v.z();
  Vec t1;
  if(std::abs(x) > std::abs(y) && std::abs(x) > std::abs(z)){
    t1 = Vec(z, y, -x).normalized();
  } else if(std::abs(y) > std::abs(x) && std::abs(y) > std::abs(z)){
    t1 = Vec(x, z, -y).normalized();
  } else {
    t1 = Vec(-z, y, x).normalized();
  }
  Base b{t1, v.prod(t1).normalized(), v.normalized()};
  return b;
}

Vec Diffuse::reflex_dir(Vec const& source, Vec const& n) const{
  double x, y, z;

  double r1 = roll();
  double r2 = roll();

  x = cos(2 * pi * r1) * sqrt(1 - r2);
  y = cos(2 * pi * r1) * sqrt(1 - r2);
  z = sqrt(r2);

  Base b = ortho_prod(n);
  Vec v1(b.x), v2(b.y), v3(b.z);

  return x * v1 + y * v2 + z * v3;
}

Vec Transparent::reflex_dir(Vec const& source, Vec const& n) const{
  return source;
}
