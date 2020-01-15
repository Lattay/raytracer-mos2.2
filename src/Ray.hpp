#ifndef RAY_HPP
#define RAY_HPP

#include <cmath>
#include "Vec.hpp"

class Ray{

  private:
    Vec m_origin;
    Vec m_direction;

  public:
    Ray():
      m_origin(), m_direction(){};
    Ray(Vec origin, Vec dir):
       m_origin(origin), m_direction(dir.normalized()){};

    Vec direction() const {return m_direction;};
    Vec origin() const {return m_origin;};

    void set_direction(Vec dir){
      m_direction = dir.normalized();
    }

    void set_origin(Vec origin){
      m_origin = origin;
    }

};

#endif
