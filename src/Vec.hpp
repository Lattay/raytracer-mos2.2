#ifndef VEC_HPP
#define VEC_HPP
#include <cmath>
#include <iostream>
#include <cassert>

class Vec {
  private:
    double m_x, m_y, m_z;

  public:

    double x() const{return m_x;}
    double y() const{return m_y;}
    double z() const{return m_z;}

    double r() const{return m_x;}
    double g() const{return m_y;}
    double b() const{return m_z;}

    double &operator[](int i){
      assert(i > 2 || i < 0);
      switch(i){
        case 0: return m_x;
        case 1: return m_y;
        case 2: return m_z;
        default: return m_x; // unreachable
      }
    }

    Vec():
      m_x(0.0), m_y(0.0), m_z(0.0) {}
    Vec(double x, double y, double z):
      m_x(x),m_y(y),m_z(z) {}

    double norm_sq() const;
    double norm() const;
    double dot(Vec const& b) const;
    Vec prod(Vec const& b) const;

    Vec normalized() const;
};

typedef struct {Vec x; Vec y; Vec z;} Base;

Vec operator*(Vec const& v, double f);
Vec operator*(double f, Vec const& v);
Vec operator*(Vec const& a, Vec const& b);
Vec operator/(Vec const& a, double f);
Vec operator+(Vec const& a, Vec const& b);
Vec operator-(Vec const& a, Vec const& b);
Vec operator-(Vec const& a);
std::ostream& operator<<(std::ostream& a, Vec const& b);

#endif
