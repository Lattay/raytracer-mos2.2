#ifndef TEXTURE_H
#define TEXTURE_H
#include "Material.hpp"

class Texture {
  private:
    unsigned char* m_data = nullptr;
    int m_w, m_h, m_n;
  public:
    Texture() {};
    Texture(const char* filename);
    ~Texture();

    Vec get_color(Vec pos) const;
};

#endif
