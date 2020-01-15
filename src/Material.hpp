#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "Vec.hpp"

class Material {
  private:
    Vec m_color;
    bool m_mirror;

  public:
    Material():
      m_color(Vec(1, 1, 1)) {};
    Material(Vec color):
      m_color(color), m_mirror(false) {};
    Material(bool mirror):
      m_color(Vec(1, 1, 1)), m_mirror(mirror) {};
    Vec color() const{return m_color;};
    bool mirror() const{return m_mirror;};
};

const Material white_emit = Material();
#endif
