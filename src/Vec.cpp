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

Vec Vec::rotate(Base b) const{
  double x = b.x.x() * m_x + b.y.x() * m_y + b.z.x() * m_z;
  double y = b.x.y() * m_x + b.y.y() * m_y + b.z.y() * m_z;
  double z = b.x.z() * m_x + b.y.z() * m_y + b.z.z() * m_z;
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
