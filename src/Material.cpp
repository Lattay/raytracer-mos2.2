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

void Texture::load(){
  FILE* f;
  f = fopen(m_filename, "rb");
  unsigned char info[54];
  fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

  int w = *(int*)&info[18]; // extract image height and width from header
  int h = *(int*)&info[22];

  int size = 3 * w * h;
  m_raw.resize(size); // allocate 3 bytes per pixel
  fread(&m_raw[0], sizeof(unsigned char), size, f); // read the rest of the data at once
  fclose(f);

  for (int i = 0; i < size; i += 3) {
    std::swap(m_raw[i], m_raw[i + 2]);
  }
}

Diffuse Texture::get_diffuse(Vec const& i, Vec const& j, Vec const& k, Vec const& p) const{
  // FIXME
  Vec color(1, 1, 1);
  return Diffuse(color);
}
