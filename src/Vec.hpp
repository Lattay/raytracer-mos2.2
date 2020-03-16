#ifndef VEC_HPP
#define VEC_HPP
#include <cmath>
#include <iostream>
#include <cassert>

typedef struct _b Base;

class Vec {
  private:
    double m_v[3];

  public:

    double x() const{return m_v[0];}
    double y() const{return m_v[1];}
    double z() const{return m_v[2];}

    double r() const{return m_v[0];}
    double g() const{return m_v[1];}
    double b() const{return m_v[2];}

    double &operator[](int i){
      assert(i < 3 && i >= 0);
      return m_v[i];
    }

    Vec():
      m_v{0.0, 0.0, 0.0} {};
    Vec(double x, double y, double z):
      m_v{x, y, z} {};

    double norm_sq() const;
    double norm() const;
    double dot(Vec const& b) const;
    Vec prod(Vec const& b) const;

    Vec normalized() const;

    Vec rotate(Base b) const;
};

typedef struct _b { Vec x, y, z; } Base;

Vec operator*(Vec const& v, double f);
Vec operator*(double f, Vec const& v);
Vec operator*(Vec const& a, Vec const& b);
Vec operator/(Vec const& a, double f);
Vec operator+(Vec const& a, Vec const& b);
Vec operator-(Vec const& a, Vec const& b);
Vec operator-(Vec const& a);
std::ostream& operator<<(std::ostream& a, Vec const& b);

Base base_from(Vec const& z);

#endif
