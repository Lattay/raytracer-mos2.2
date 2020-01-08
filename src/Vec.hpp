#ifndef VEC_HPP
#define VEC_HPP
#include <cmath>

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

    Vec():
      m_x(), m_y(), m_z() {}
    Vec(double x, double y, double z):
      m_x(x),m_y(y),m_z(z) {}

    double norm_sq() const;
    double norm() const;
    double dot(Vec const& b) const;
    Vec prod(Vec const& b) const;

    Vec normalized() const;
};

Vec operator*(Vec const& v, double f);
Vec operator*(double f, Vec const& v);
Vec operator/(Vec const& a, double f);
Vec operator+(Vec const& a, Vec const& b);
Vec operator-(Vec const& a, Vec const& b);
Vec operator-(Vec const& a);

#endif
