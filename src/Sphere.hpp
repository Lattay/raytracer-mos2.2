#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "Vec.hpp"
#include "Ray.hpp"

class Sphere {

  private:
    Vec m_origin;
    double m_radius;

  public:
    Sphere(Vec origin, double radius):
      m_origin(origin), m_radius(radius) {};

    double r() const{return m_radius;};
    Vec origin() const{return m_origin;};
    bool intersect(Ray const& r) const;
    Vec intersection(Ray const& r) const;
};

#endif
