#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "Vec.hpp"

typedef enum { MIRROR, TRANSPARENT} SFX;

class Material {

  public:
    Vec color() const{return Vec(0, 0, 0);};
    bool direct_lighting() const{return false;};
    bool transparent() const{return false;};
    double index() const{return 1.0;};
    virtual Vec reflex_dir(Vec const& source, Vec const& n) const;
};

class Diffuse: public Material {
  private:
    Vec m_color;
  public:
    Diffuse(Vec color): m_color(color){};
    bool direct_lighting() const{return true;};
    Vec color() const{return m_color;};
    Vec reflex_dir(Vec const& source, Vec const& n) const;
};

class Reflective: public Material {
  public:
    Reflective(){};
    Vec reflex_dir(Vec const& source, Vec const& n) const;
};

class Transparent: public Material {
  private:
    double m_index;
  public:
    Transparent(double index): m_index(index) {};
    double index() const{return m_index;};
    bool transparent() const{return true;};
    Vec reflex_dir(Vec const& source, Vec const& n) const;
};

const Diffuse white_emit(Vec(1, 1, 1));

#endif
