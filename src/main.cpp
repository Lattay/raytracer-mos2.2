#include "main.hpp" 
 
const double pi = 3.141592654;
 
int main() {
  int W = 512;
  int H = 512;

  double fov = pi/3.0;

  Vec origin(0, 0, 55);
  Vec light(-10, 20, 40);
  Vec c(0, 0, 0);
  Sphere s(c, 10);
  Sphere s2(-c, 10);

  double I0 = 255;
  double rho = 1;

  std::vector<unsigned char> image(W*H * 3, 0);
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      double x, y, z;
      x = j - W/2.0 + 0.5;
      y = i - H/2 + 0.5;
      z = -W/(2*tan(fov/2.0));

      Vec dir = Vec(x, y, z);
      Ray r(origin, dir.normalized());

      if(s.intersect(r)){
        // Vec p = s.intersection(r);

        // Vec n = (p - c).normalized();
        // Vec vl = p - light;

        // int I = ceil(I0 * rho/pi * std::max(0.0, n.dot(vl.normalized())) / vl.norm_sq());

        image[(i*W + j) * 3 + 0] = 255;
        image[(i*W + j) * 3 + 1] = 255;
        image[(i*W + j) * 3 + 2] = 255;
      } else {
        image[(i*W + j) * 3 + 0] = 0;
        image[(i*W + j) * 3 + 1] = 0;
        image[(i*W + j) * 3 + 2] = 0;
      }
    }
  }
  stbi_write_png("image.png", W, H, 3, &image[0], 0);

  return 0;
}
