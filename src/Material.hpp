#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "Vec.hpp"

typedef enum { MIRROR, TRANSPARENT} SFX;

class Material {
  private:
    Vec m_color;
    bool m_mirror;
    bool m_transparent;

  public:
    Material(): m_color(Vec(1, 1, 1)), m_mirror(false), m_transparent(false) {};
    Material(Vec color):
      m_color(color), m_mirror(false), m_transparent(false) {};
    Material(SFX effect):
      m_color(Vec(1, 1, 1)), m_mirror(effect == MIRROR), m_transparent(effect == TRANSPARENT) {};
    Material(Vec color, SFX effect):
      m_color(color), m_mirror(effect == MIRROR), m_transparent(effect == TRANSPARENT) {};
    Vec color() const{return m_color;};
    bool mirror() const{return m_mirror;};
    bool transparent() const{return m_transparent;};
    double index() const{return 1.5;};
};

const Material white_emit = Material();
#endif
