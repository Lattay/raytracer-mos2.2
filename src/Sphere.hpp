#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "Vec.hpp"
#include "Ray.hpp"
#include "Material.hpp"
#include "Intersection.hpp"

class Sphere {

  private:
    Vec m_origin;
    Material m_material;
    double m_radius;

  public:
    Sphere(Vec origin, double radius):
      m_origin(origin), m_material(white_emit), m_radius(radius) {};
    Sphere(Vec origin, double radius, Material material):
      m_origin(origin), m_material(material), m_radius(radius) {};

    double r() const{return m_radius;};
    Vec origin() const{return m_origin;};
    Material material() const{return m_material;};
    bool intersect(Ray const& r) const;
    Intersection intersection(Ray const& r) const;
};

#endif
