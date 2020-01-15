#include <iostream>
#include "Scene.hpp"

const Vec black(0, 0, 0);

Intersection Scene::intersection(Ray const& ray) const{

  double min_sq_dist = 1e10;
  // Intersection empty_intersection = Intersection();
  Intersection final_intersect;

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
  const double epsilon = 0.001;

  Intersection inter = intersection(ray);

  if(!inter.valid()){
    return black;
  }


  if(inter.material().mirror()){
    if(k < 0){
      return black;
    } else {
      Vec const& n = inter.normal();
      Ray new_ray(inter.position() + epsilon * n, ray.direction() - 2.0 * n.dot(ray.direction()) * n);

      return get_color(new_ray, source, k - 1);
    }
  } else if(inter.material().transparent() && k >= 0){
    if(k < 0){
      return black;
    } else {
      Vec n = inter.normal();
      Vec i = ray.direction();
      // Vec t = (i + n).normalized();

      double n2_n1 = inter.material().index();
      double i_dot_n = i.dot(n);
      double cos2_theta = 1 - 1.0/(n2_n1 * n2_n1) * (1 - i_dot_n*i_dot_n);

      Ray new_ray;
      if(cos2_theta > 0){
        // transmis
        Ray refracted_ray = Ray(inter.position() - epsilon * n,
            - sqrt(cos2_theta) * n + (i - i_dot_n * n) / n2_n1);

        Intersection inter2 = intersection(refracted_ray);

        if(!inter2.valid()){
          // error
          return black;
        } else {
          i = refracted_ray.direction();
          n = inter2.normal();

          double i_dot_n = i.dot(n);
          double cos2_theta = 1 - n2_n1 * n2_n1 * (1 - i_dot_n*i_dot_n);

          new_ray = Ray(inter2.position() - epsilon * n,
              - sqrt(cos2_theta) * n + (i - i_dot_n * n) * n2_n1);
        }
      } else {
        // reflechis
        new_ray.set_origin(inter.position() + epsilon * n);
        new_ray.set_direction(ray.direction() - 2.0 * n.dot(ray.direction()) * n);
      }

    return get_color(new_ray, source, k - 1);
    }
  } else {
    Vec vl = source.source() - inter.position();

    Ray r_light = Ray(inter.position() + epsilon * inter.normal(), vl.normalized());

    double d2 = (inter.position() - source.source()).norm_sq();

    for(long unsigned int i = 0; i < m_spheres.size(); ++i){
      Sphere const& s = *m_spheres[i];
      Intersection i_light = s.intersection(r_light);
      if(i_light.valid()){
        if(d2 > (i_light.position() - inter.position()).norm_sq())
          return black;
      }
    }
    return inter.material().color() * source.intensity() * std::max(0.0, inter.normal().dot(vl.normalized())) / vl.norm_sq();
  }
}

int Scene::add_sphere(Sphere const& s){
  m_spheres.push_back(&s);
  return m_spheres.size() - 1;
}

int Scene::add_new_sphere(Sphere s){
  m_spheres.push_back(new Sphere(s));
  return m_spheres.size() - 1;
}
