#include "Vec.hpp"


Vec operator*(Vec const& a, double f){
  return Vec(a.x()*f, a.y()*f, a.z()*f);
}

Vec operator*(double f, Vec const& a){
  return Vec(a.x()*f, a.y()*f, a.z()*f);
}

Vec operator/(Vec const& a, double f){
  return Vec(a.x()/f, a.y()/f, a.z()/f);
}

Vec operator+(Vec const& a, Vec const& b){
  return Vec(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
}

Vec operator-(Vec const& a, Vec const& b){
  return Vec(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
}

Vec operator-(Vec const& a){
  return Vec(-a.x(), -a.y(), -a.z());
}

double Vec::dot(Vec const& b) const{
  return m_x * b.x() + m_y * b.y() + m_z * b.z();
}

double Vec::norm_sq() const{
  return m_x * m_x + m_y * m_y + m_z * m_z;
}

double Vec::norm() const{
  return sqrt(norm_sq());
}

Vec Vec::prod(Vec const& b) const{
  return Vec(
      m_y * b.z() - m_z * b.y(),
      m_z * b.x() - m_x * b.z(),
      m_x * b.y() - m_y * b.x()
  );
}

Vec Vec::normalized() const{
  return (*this)/norm();
}
