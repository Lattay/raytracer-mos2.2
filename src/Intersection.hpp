#ifndef INTERSECTION_HPP
#define INTERSECTION_HPP
#include "Vec.hpp"
#include "Ray.hpp"
#include "Material.hpp"

class Intersection {
  private:
    bool m_valid;
    Vec m_position;
    Vec m_normal;
    Material m_material;

  public:
    Intersection():
      m_valid(false), m_position(), m_normal(), m_material() {};
    Intersection(Vec pos, Vec norm, Material mat):
      m_valid(true), m_position(pos), m_normal(norm), m_material(mat) {};

    bool valid() const{return m_valid;};
    Vec position() const{return m_position;};
    Vec normal() const{return m_normal;};
    Material material() const{return m_material;};
};
#endif
