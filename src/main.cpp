#include "main.hpp" 
 
const double pi = 3.141592654;

// Camera
const int W = 512;
const int H = 512;

const double fov = pi/3.0;

const Vec origin(0, 0, 55);

// Light source
const Light light(Vec(-10, 20, 40), 1e8);

// colors
const Material mirror(MIRROR);
const Material glass(TRANSPARENT);
const Material black(Vec(0, 0, 0));
const Material red(Vec(1, 0, 0));
const Material blue(Vec(0, 0, 1));
const Material green(Vec(0, 1, 0));
const Material purple(Vec(1, 0, 1));
const Material white(Vec(1, 1, 1));

// Scene
const Vec c(0, 0, 0);

int main() {
  Scene scene;

  // scene.add_new_sphere(Sphere(c, 10, white));
  scene.add_new_sphere(Sphere(c, 10, glass));
  scene.add_new_sphere(Sphere(c + Vec(-8, 8, 10), 3, mirror));

  scene.add_new_sphere(Sphere(Vec(0, 1000, 0), 940, red));
  scene.add_new_sphere(Sphere(Vec(0, 0, -1000), 940, green));
  scene.add_new_sphere(Sphere(Vec(0, -1000, 0), 990, blue));
  // This one should be invisible unless there is a bug (or a reflexion)
  scene.add_new_sphere(Sphere(Vec(0, 0, 1000), 940, purple));

  std::vector<float> image(W*H * 3, 0);
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      double x, y, z;
      x = j - W/2.0 + 0.5;
      y = - i + H/2 + 0.5;
      z = -W/(2*tan(fov/2.0));

      Vec dir = Vec(x, y, z).normalized();
      Ray r(origin, dir);

      Vec color = scene.get_color(r, light);

      image[(i * W + j) * 3 + 0] = color.r();
      image[(i * W + j) * 3 + 1] = color.g();
      image[(i * W + j) * 3 + 2] = color.b();
    }
  }

  // Convert float intensity to uchar pixel value
  std::vector<unsigned char> file_buffer(W*H * 3, 0);
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {
      file_buffer[(i * W + j) * 3 + 0] = std::max(0, std::min(255, (int) ceil(pow(image[(i * W + j) * 3 + 0], 0.45))));
      file_buffer[(i * W + j) * 3 + 1] = std::max(0, std::min(255, (int) ceil(pow(image[(i * W + j) * 3 + 1], 0.45))));
      file_buffer[(i * W + j) * 3 + 2] = std::max(0, std::min(255, (int) ceil(pow(image[(i * W + j) * 3 + 2], 0.45))));
    }
  }

  stbi_write_png("image.png", W, H, 3, &file_buffer[0], 0);

  return 0;
}
