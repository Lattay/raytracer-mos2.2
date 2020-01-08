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

  double I0 = 500000;

  double minI = 1000000;
  double maxI = 0;


  std::vector<float> image(W*H * 3, 0);
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      double x, y, z;
      x = j - W/2.0 + 0.5;
      y = - i + H/2 + 0.5;
      z = -W/(2*tan(fov/2.0));

      Vec dir = Vec(x, y, z).normalized();
      Ray r(origin, dir);

      if(s.intersect(r)){
        Vec p = s.intersection(r);

        Vec n = (p - s.origin()).normalized();
        Vec vl = light - p;

        double I = I0 * std::max(0.0, n.dot(vl.normalized())) / vl.norm_sq();

        image[(i*W + j) * 3 + 0] = I;
        image[(i*W + j) * 3 + 1] = I;
        image[(i*W + j) * 3 + 2] = I;

        if(I > maxI){
          maxI = I;
        }

        if(I < minI){
          minI = I;
        }
      } else {
        image[(i*W + j) * 3 + 0] = 0.0;
        image[(i*W + j) * 3 + 1] = 0.0;
        image[(i*W + j) * 3 + 2] = 0.0;
      }
    }
  }

  std::cout << "I min max " << minI << " " << maxI << std::endl;

  // Convert float intensity to uchar pixel value
  std::vector<unsigned char> file_buffer(W*H * 3, 0);
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {
      file_buffer[(i * W + j) * 3 + 0] = std::max(0, std::min(255, (int) ceil(image[(i * W + j) * 3 + 0])));
      file_buffer[(i * W + j) * 3 + 1] = std::max(0, std::min(255, (int) ceil(image[(i * W + j) * 3 + 1])));
      file_buffer[(i * W + j) * 3 + 2] = std::max(0, std::min(255, (int) ceil(image[(i * W + j) * 3 + 2])));
    }
  }

  stbi_write_png("image.png", W, H, 3, &file_buffer[0], 0);

  return 0;
}
