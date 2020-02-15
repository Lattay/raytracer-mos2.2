#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "Vec.hpp"
#include "Ray.hpp"
#include "Material.hpp"
#include "Intersection.hpp"
#include "Object.hpp"

class Sphere: public Object {

  private:
    Vec m_origin;
    Material const& m_material;
    double m_radius;

  public:
    Sphere(Vec origin, double radius):
      m_origin(origin), m_material(white_emit), m_radius(radius) {};
    Sphere(Vec origin, double radius, Material const& material):
      m_origin(origin), m_material(material), m_radius(radius) {};

    double r() const{return m_radius;};
    Vec origin() const{return m_origin;};

    virtual Material const& material() const{return m_material;};
    // virtual bool intersect(Ray const& r) const;
    virtual Intersection intersection(Ray const& r) const;
};

#endif
