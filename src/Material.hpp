#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "Vec.hpp"

class Material {
  private:
    Vec m_color;

  public:
    Material():
      m_color(Vec(1, 1, 1)) {};
    Material(Vec color):
      m_color(color) {};
    Vec color() const{return m_color;};
};

const Material white_emit = Material();
#endif
