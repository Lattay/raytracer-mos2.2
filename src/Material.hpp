#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include <vector>
#include "Vec.hpp"
#include "random_tools.hpp"

class Material {

  public:
    virtual ~Material() {};
    virtual Vec color() const{return Vec(0, 0, 0);};
    virtual bool direct_lighting() const{return false;};
    virtual bool transparent() const{return false;};
    virtual double index() const{return 1.0;};
    virtual Sample reflex_dir(Vec const& source, Vec const& n) const = 0;
};

class Diffuse: public Material {
  private:
    Vec m_color;
  public:
    Diffuse(Vec color): m_color(0.9 * color){};
    ~Diffuse() {};
    virtual bool direct_lighting() const{return true;};
    virtual Vec color() const{return m_color;};
    virtual Sample reflex_dir(Vec const& source, Vec const& n) const;
};

class Reflective: public Material {
  private:
    double m_loss;
  public:
    Reflective(double loss): m_loss(loss){};
    ~Reflective() {};
    virtual Sample reflex_dir(Vec const& source, Vec const& n) const;
    virtual Vec color() const{return {1, 1, 1};};
};

class Transparent: public Material {
  private:
    double m_index;
  public:
    Transparent(double index): m_index(index) {};
    ~Transparent() {};
    virtual double index() const{return m_index;};
    virtual bool transparent() const{return true;};
    virtual Sample reflex_dir(Vec const& source, Vec const& n) const;
};

const Diffuse white_emit(Vec(1, 1, 1));

#endif
