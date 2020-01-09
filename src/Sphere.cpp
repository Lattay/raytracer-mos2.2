#include "Sphere.hpp"
#include <iostream>

bool Sphere::intersect(Ray const& r) const{
  // solve t^2 + b*t + c
  double b = 2 * r.direction().dot(r.origin() - m_origin);
  double c = (r.origin() - m_origin).norm_sq() - m_radius * m_radius;
  double delta = b * b - 4 * c;
  if(delta < 0){
    return false;
  }

  double t1 = (-b - sqrt(delta)) / 2.0;
  double t2 = (-b + sqrt(delta)) / 2.0;
  return (t1 >= 0) || (t2 >= 0);
}

Intersection Sphere::intersection(Ray const& r) const{
  // solve t^2 + b*t + c
  double b = 2 * r.direction().dot(r.origin() - m_origin);
  double c = (r.origin() - m_origin).norm_sq() - m_radius * m_radius;
  double delta = b * b - 4 * c;
  if(delta < 0){
    return Intersection();
  }

  double t1 = (-b - sqrt(delta)) / 2.0;
  double t2 = (-b + sqrt(delta)) / 2.0;

  if(t2 < 0){
    return Intersection();
  }

  Vec p = r.origin() + (t1 > 0.0 ? t1 : t2) * r.direction();

  return Intersection(p, (p - m_origin).normalized(), m_material);
}
