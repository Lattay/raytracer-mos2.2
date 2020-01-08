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

    // add a sphere to the scene with a reference
    int add_sphere(Sphere const& s);
    // add a sphere to the scene by allocating it
    int add_new_sphere(Sphere s);
};

#endif
