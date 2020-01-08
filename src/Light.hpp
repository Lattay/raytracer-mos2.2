#ifndef LIGHT_HPP
#define LIGHT_HPP
#include "Vec.hpp"

class Light {
  private:
    Vec m_source;
    double m_intensity;

  public:
    Light(Vec v, double i):
      m_source(v), m_intensity(i) {};
    Vec source() const{return m_source;};
    double intensity() const{return m_intensity;};
};

#endif
