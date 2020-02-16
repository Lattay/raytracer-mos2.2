#ifndef OBJECT_HPP

#define OBJECT_HPP
#include "Vec.hpp"
#include "Ray.hpp"
#include "Material.hpp"
#include "Intersection.hpp"

class Object {
  public:
    virtual Intersection intersection(Ray const& r) const = 0;
};

#endif
