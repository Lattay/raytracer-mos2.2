#include "Vec.hpp"


Vec operator*(Vec const& a, double f){
  return Vec(a.x()*f, a.y()*f, a.z()*f);
}

Vec operator*(double f, Vec const& a){
  return Vec(a.x()*f, a.y()*f, a.z()*f);
}

Vec operator*(Vec const& a, Vec const& b){
  return Vec(a.x()*b.x(), a.y()*b.y(), a.z()*b.z());
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

std::ostream& operator<<(std::ostream& a, Vec const& b){
  a << "(" << b.x() << ", " << b.y() << ", " << b.z() << ")";
  return a;
}

double Vec::dot(Vec const& b) const{
  return m_v[0] * b.m_v[0] + m_v[1] * b.m_v[1] + m_v[2] * b.m_v[2];
}

double Vec::norm_sq() const{
  return m_v[0] * m_v[0] + m_v[1] * m_v[1] + m_v[2] * m_v[2];
}

double Vec::norm() const{
  return sqrt(norm_sq());
}

Vec Vec::prod(Vec const& b) const{
  return Vec(
      m_v[1] * b.m_v[2] - m_v[2] * b.m_v[1],
      m_v[2] * b.m_v[0] - m_v[0] * b.m_v[2],
      m_v[0] * b.m_v[1] - m_v[1] * b.m_v[0]
  );
}

Vec Vec::normalized() const{
  return (*this)/norm();
}

Vec Vec::rotate(Base b) const{
  double x = b.x.m_v[0] * m_v[0] + b.y.m_v[0] * m_v[1] + b.z.m_v[0] * m_v[2];
  double y = b.x.m_v[1] * m_v[0] + b.y.m_v[1] * m_v[1] + b.z.m_v[1] * m_v[2];
  double z = b.x.m_v[2] * m_v[0] + b.y.m_v[2] * m_v[1] + b.z.m_v[2] * m_v[2];
  return Vec(x, y, z);
}

/* Create a base where xOz is the same plane as in canonical base
 * This way, when the camera looks in Oz+ direction whith 0x+ as vertical
 * the view is not rotated
 */
Base base_from(Vec const& z){
  Vec x = Vec(0, 1, 0).prod(z);
  return {x.normalized(), z.prod(x).normalized(), z.normalized()};
}
