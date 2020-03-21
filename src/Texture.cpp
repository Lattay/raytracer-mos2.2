#include "../vendor/stb_image.h"
#include "Texture.hpp"

Texture::~Texture(){
  if(m_data != nullptr)
    delete m_data;
}

Texture::Texture(const char* filename){
  m_data = stbi_load(filename, &m_w, &m_h, &m_n, 0);
  std::cout << "Loaded texture at " << filename
    << " (" << m_w << ", " << m_h << ")" << std::endl;
}

Vec Texture::get_color(Vec pos) const{
  int x = pos.x()*m_w;
  int y = (1.0 - pos.y())*m_h;

  Vec color;
  color[0] = m_data[m_n*(x + y * m_w)];
  color[1] = m_data[m_n*(x + y * m_w) + 1];
  color[2] = m_data[m_n*(x + y * m_w) + 2];
  return color/255.0;
}
