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
    Vec get_color(Ray const& ray, Light const& source, int depth) const;
    Vec get_color(Ray const& ray, Light const& source, int k, bool inside) const;

    int add_mesh(Mesh const& mesh);
    int add_sphere(Sphere const& s);
    // the following allocate the object themselves
    int add_new_mesh(Mesh const& mesh);
    int add_new_sphere(Sphere const& s);
};

#endif
