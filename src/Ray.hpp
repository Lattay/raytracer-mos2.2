#ifndef RAY_HPP
#define RAY_HPP

#include <cmath>
#include "Vec.hpp"

class Ray{

  private:
    Vec m_origin;
    Vec m_direction;

  public:
    Ray(Vec origin, Vec dir):
       m_origin(origin), m_direction(dir){};

    Vec direction() const {return m_direction;};
    Vec origin() const {return m_origin;};
};

#endif
