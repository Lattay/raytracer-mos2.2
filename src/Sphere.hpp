#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "Vec.hpp"
#include "Ray.hpp"

class Sphere {

  private:
    Vec m_origin;
    Vec m_albedo;
    double m_radius;

  public:
    Sphere(Vec origin, double radius):
      m_origin(origin), m_albedo(Vec(255, 255, 255)), m_radius(radius) {};
    Sphere(Vec origin, double radius, Vec albedo):
      m_origin(origin), m_albedo(albedo), m_radius(radius) {};

    double r() const{return m_radius;};
    Vec origin() const{return m_origin;};
    Vec albedo() const{return m_albedo;};
    bool intersect(Ray const& r) const;
    Vec intersection(Ray const& r) const;
};

#endif
