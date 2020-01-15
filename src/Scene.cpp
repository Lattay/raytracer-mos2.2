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

Vec Scene::get_color(Ray const& ray, Light const& source) const{
  return get_color(ray, source, 5);
}

Vec Scene::get_color(Ray const& ray, Light const& source, int k) const{

  Intersection inter = intersection(ray);

  if(!inter.valid()){
    return Vec(0, 0, 0);
  }

  Vec vl = source.source() - inter.position();

  Ray r_light = Ray(inter.position() + 0.001 * inter.normal(), vl.normalized());

  double d2 = (inter.position() - source.source()).norm_sq();

  for(long unsigned int i = 0; i < m_spheres.size(); ++i){
    Sphere const& s = *m_spheres[i];
    Intersection i_light = s.intersection(r_light);
    if(i_light.valid()){
      if(d2 > (i_light.position() - inter.position()).norm_sq())
        return Vec(0, 0, 0);
    }
  }
  
  return inter.material().color() * source.intensity() * std::max(0.0, inter.normal().dot(vl.normalized())) / vl.norm_sq();
}

int Scene::add_sphere(Sphere const& s){
  m_spheres.push_back(&s);
  return m_spheres.size() - 1;
}

int Scene::add_new_sphere(Sphere s){
  m_spheres.push_back(new Sphere(s));
  return m_spheres.size() - 1;
}
