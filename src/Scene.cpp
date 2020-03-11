#include <iostream>
#include "Scene.hpp"
#include "random_tools.hpp"

// 1/(2pi) = 0.159154
// 1/pi = 0.318310
#define ALPHA 0.318310

const Vec black(0, 0, 0);

static double shlick(double n1, double n2, double cos_theta){
  double r0 = (n1 - n2)/(n1 + n2);
  return r0 * r0 + (1 + r0 * r0) * pow(1 - cos_theta, 5);
}

Intersection Scene::intersection(Ray const& ray) const{

  double min_sq_dist = 1e10;
  // Intersection empty_intersection = Intersection();
  Intersection final_intersect;

  for(long unsigned int i = 0; i < m_objects.size(); ++i){
    Object const& o = *m_objects[i];
    Intersection intersect = o.intersection(ray);
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

Vec Scene::get_color(Ray const& ray, Light const& source, int depth) const{
  return get_color(ray, source, depth, false);
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

    Vec indirect(black), direct(black);

    if(k >= 0){
      Vec const& n = inter.normal();
      Sample s = inter.material().reflex_dir(ray.direction(), n);
      Ray new_ray(inter.position() + epsilon * n, s.dir);

      Vec c = get_color(new_ray, source, k - 1, false);

      indirect = inter.material().color() * c;
    }

    if(inter.material().direct_lighting()){
      // Diffusive material have direct lighting
      Vec vl = source.source() - inter.position();
      Sample s(random_cos(-vl));
      Vec p_emit(source.source() + s.dir*source.radius());

      Vec vl2 = p_emit - inter.position();

      double d2 = vl2.norm_sq();
      vl2 = vl2.normalized();

      Ray r_light = Ray(inter.position() + epsilon * inter.normal(), vl2);

      bool shadowed = false;
      for(long unsigned int i = 0; i < m_objects.size(); ++i){
        Object const& o = *m_objects[i];
        Intersection i_light = o.intersection(r_light);
        if(i_light.valid()){
          if(d2 > (i_light.position() - inter.position()).norm_sq()){
            shadowed = true;
            break;
          }
        }
      }
      if(!shadowed){
        double intensity = source.intensity()
          * std::max(0.0, inter.normal().dot(vl2))
          * std::max(0.0, s.dir.normalized().dot(-vl2));
        direct = inter.material().color() * intensity  / (pi * d2 * s.proba);
      }
    }
    return direct + indirect;
  }
}

int Scene::add_sphere(Sphere const& s){
  m_objects.push_back(&s);
  return m_objects.size() - 1;
}

int Scene::add_new_sphere(Sphere const& s){
  m_objects.push_back(new Sphere(s));
  return m_objects.size() - 1;
}

int Scene::add_mesh(Mesh const& mesh){
  m_objects.push_back(&mesh);
  return m_objects.size() - 1;
}

int Scene::add_new_mesh(Mesh const& mesh){
  m_objects.push_back(new Mesh(mesh));
  return m_objects.size() - 1;
}
