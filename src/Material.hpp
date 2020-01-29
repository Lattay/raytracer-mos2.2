#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "Vec.hpp"

typedef enum { MIRROR, TRANSPARENT} SFX;

class Material {

  public:
    virtual Vec color() const{return Vec(0, 0, 0);};
    virtual bool direct_lighting() const{return false;};
    virtual bool transparent() const{return false;};
    virtual double index() const{return 1.0;};
    virtual Vec reflex_dir(Vec const& source, Vec const& n) const = 0;
};

class Diffuse: public Material {
  private:
    Vec m_color;
  public:
    Diffuse(Vec color): m_color(color){};
    virtual bool direct_lighting() const{return true;};
    virtual Vec color() const{return m_color;};
    virtual Vec reflex_dir(Vec const& source, Vec const& n) const;
};

class Reflective: public Material {
  public:
    Reflective(){};
    virtual Vec reflex_dir(Vec const& source, Vec const& n) const;
};

class Transparent: public Material {
  private:
    double m_index;
  public:
    Transparent(double index): m_index(index) {};
    virtual double index() const{return m_index;};
    virtual bool transparent() const{return true;};
    virtual Vec reflex_dir(Vec const& source, Vec const& n) const;
};

const Diffuse white_emit(Vec(1, 1, 1));

#endif
