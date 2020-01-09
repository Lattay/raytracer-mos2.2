#ifndef LIGHT_HPP
#define LIGHT_HPP
#include "Vec.hpp"

class Light {
  private:
    Vec m_source;
    Vec m_color;
    double m_intensity;

  public:
    Light(Vec v, double i):
      m_source(v), m_color(Vec(1, 1, 1)), m_intensity(i) {};
    Light(Vec v, double i, Vec color):
      m_source(v), m_color(color), m_intensity(i) {};
    Vec source() const{return m_source;};
    double intensity() const{return m_intensity;};
};

#endif
