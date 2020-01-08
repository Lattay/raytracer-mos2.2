#include "Scene.hpp"

Vec Scene::get_color(Ray const& ray) const{
  Vec color(0.0, 0.0, 0.0);
  double min_sq_dist = 1e10;
  Vec closest_p;
  int closest_s;
  bool intersect = false;

  for(long unsigned int i = 0; i < m_spheres.size(); ++i){
    Sphere const& s = *m_spheres[i];
    if(s.intersect(ray)){
      intersect = true;
      Vec p = s.intersection(ray);
      double d = (p - ray.origin()).norm_sq();
      if(d < min_sq_dist){
        min_sq_dist = d;
        closest_p = p;
        closest_s = i;
      }
    }
  }

  if(intersect){
    Sphere sphere = *m_spheres[closest_s];
    Vec n = (closest_p - sphere.origin()).normalized();
    Vec vl = m_light.source() - closest_p;

    color = sphere.albedo() * m_light.intensity() * std::max(0.0, n.dot(vl.normalized())) / vl.norm_sq();
  }
  return color;
}

int Scene::add_sphere(Sphere const& s){
  m_spheres.push_back(&s);
  return m_spheres.size() - 1;
}

int Scene::add_new_sphere(Sphere s){
  m_spheres.push_back(new Sphere(s));
  return m_spheres.size() - 1;
}
