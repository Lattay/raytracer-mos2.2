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
    bool m_mesh;
    bool m_own_mat = false;
    Material const* m_material;

  public:
    Intersection():
      m_valid(false), m_position(), m_normal(), m_mesh(false), m_material(&white_emit) {};
    Intersection(Vec pos, Vec norm, Material const& mat):
      m_valid(true), m_position(pos), m_normal(norm.normalized()), m_mesh(false), m_own_mat(false), m_material(&mat) {};
    Intersection(Vec pos, Vec norm, Material* mat):
      m_valid(true), m_position(pos), m_normal(norm.normalized()), m_mesh(false), m_own_mat(true), m_material(mat) {};

    ~Intersection(){
      /*
      if(m_own_mat)
        delete m_material;
      */
    }

    bool valid() const{return m_valid;};
    Vec position() const{return m_position;};
    Vec normal() const{return m_normal;};
    Material const& material() const{return *m_material;};
    bool is_mesh(bool yes) { return m_mesh = yes; };
    bool is_mesh() { return m_mesh; };
};

const Intersection invalid_intersect();
#endif
