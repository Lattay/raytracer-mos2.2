#ifndef SCENE_HPP
#define SCENE_HPP
#include <vector>
#include "Vec.hpp"
#include "Ray.hpp"
#include "Object.hpp"
#include "Sphere.hpp"
#include "Mesh.hpp"
#include "Light.hpp"
#include "Intersection.hpp"

class Scene {

  private:
    std::vector<Object const*> m_objects;

  public:
    Scene():
      m_objects() {};
    Intersection intersection(Ray const& ray) const;
    Vec get_color(Ray const& ray, Light const& sourcr) const;
    Vec get_color(Ray const& ray, Light const& source, int k, bool inside) const;

    // add a sphere to the scene with a reference
    int add_sphere(Sphere const& s);
    // add a sphere to the scene by allocating it
    int add_new_sphere(Sphere s);
    int add_mesh(MeshBox* mesh);
};

#endif
