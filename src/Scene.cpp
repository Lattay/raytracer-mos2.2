#include "Scene.hpp"

Intersection Scene::intersection(Ray const& ray) const{
  double min_sq_dist = 1e10;
  Intersection final_intersect = Intersection();

  for(long unsigned int i = 0; i < m_spheres.size(); ++i){
    Sphere const& s = *m_spheres[i];
    Intersection intersect = s.intersection(ray);
    if(intersect.valid()){
      double d = (intersect.position() - ray.origin()).norm_sq();
      if(d < min_sq_dist){
        final_intersect = intersect;
        min_sq_dist = d;
      }
    }
  }
  return final_intersect;
}

int Scene::add_sphere(Sphere const& s){
  m_spheres.push_back(&s);
  return m_spheres.size() - 1;
}

int Scene::add_new_sphere(Sphere s){
  m_spheres.push_back(new Sphere(s));
  return m_spheres.size() - 1;
}
