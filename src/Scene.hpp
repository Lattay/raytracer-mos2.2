#ifndef SCENE_HPP
#define SCENE_HPP
#include <vector>
#include "Vec.hpp"
#include "Ray.hpp"
#include "Sphere.hpp"
#include "Light.hpp"

class Scene {

  private:
    Light m_light;
    std::vector<Sphere const*> m_spheres;

  public:
    Scene(Light light):
      m_light(light), m_spheres() {};
    Vec get_color(Ray const& ray) const;
    int add_sphere(Sphere const& s);
};

#endif
