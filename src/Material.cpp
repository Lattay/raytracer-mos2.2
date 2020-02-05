#include "Material.hpp"
#include "random_tools.hpp"

Sample Diffuse::reflex_dir(Vec const& source, Vec const& n) const{
  return random_cos(n);
}

Sample Transparent::reflex_dir(Vec const& source, Vec const& n) const{
  return {source, 1.0};
}

Sample Reflective::reflex_dir(Vec const& source, Vec const& n) const{
  Sample s{source - 2.0 * n.dot(source) * n, 1.0};
  return s;
}
