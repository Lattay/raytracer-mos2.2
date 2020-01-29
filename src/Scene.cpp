#include <iostream>
#include <random>
#include <omp.h>
#include "Scene.hpp"

const Vec black(0, 0, 0);

static std::default_random_engine generator[4];
static std::uniform_real_distribution<double> distribution(0,1);

static double roll(){
  return distribution(generator[omp_get_thread_num()]);
}

static double shlick(double n1, double n2, double cos_theta){
  double r0 = (n1 - n2)/(n1 + n2);
  return r0 * r0 + (1 + r0 * r0) * pow(1 - cos_theta, 5);
}

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
  return get_color(ray, source, 5, false);
}

Vec Scene::get_color(Ray const& ray, Light const& source, int k, bool inside) const{
  const double epsilon = 0.001;

  Intersection inter = intersection(ray);

  if(!inter.valid()){
    return black;
  }


  if(inter.material().transparent()){
    // purely transparent/refractive material
    
    if(k < 0){
      return black;
    } else {

      bool changing = true;

      Vec n = (inside ? -1 : 1) * inter.normal();
      Vec i = ray.direction();

      double i_dot_n = i.dot(n);

      Vec sin_theta1_t = i - i_dot_n * n;

      double sin2_theta1 = sin_theta1_t.norm_sq();

      double n1_n2 = inside ? inter.material().index() : 1.0/inter.material().index();

      double cos2_theta2 = 1 - n1_n2*n1_n2 * sin2_theta1;
      double cos_theta1 = sqrt(1 - sin2_theta1);
      double reflexion_probability = shlick(
          inside ? inter.material().index() : 1.0,
          !inside ? inter.material().index() : 1.0,
          cos_theta1);

      Ray new_ray;

      if(cos2_theta2 < 0 || roll() > reflexion_probability){
        // transmis
        new_ray = Ray(inter.position() - epsilon * n, - sqrt(cos2_theta2) * n + n1_n2 * sin_theta1_t);
      } else {
        // reflechis
        changing = false;
        new_ray = Ray(inter.position() + epsilon * n, i - 2.0 * i_dot_n * n);
      }

      return get_color(new_ray, source, k - 1, changing ? !inside : inside);
    }
  } else {

    Vec indirect(0, 0, 0), direct(0, 0, 0);

    if(k < 0){
      indirect = black;
    } else {
      Vec const& n = inter.normal();
      Ray new_ray(inter.position() + epsilon * n, inter.material().reflex_dir(ray.direction(), n));

      indirect = get_color(new_ray, source, k - 1, false);
    }

    if(inter.material().direct_lighting()){
      // Diffusive material have direct lighting
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
      direct = inter.material().color() * source.intensity() * std::max(0.0, inter.normal().dot(vl.normalized())) / vl.norm_sq();

    }
    return direct + indirect;
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
