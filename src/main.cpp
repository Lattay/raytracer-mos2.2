#include "main.hpp" 
 
const double pi = 3.141592654;

const int ray_number = 20;

// Camera
const int W = 512;
const int H = 512;

const double fov = pi/3.0;

const Vec origin(0, 0, 55);

// Light source
const Light light(Vec(-10, 20, 40), 5e8);

// colors
const Material mirror(MIRROR);
const Material glass(TRANSPARENT);
const Material water(TRANSPARENT, 1.2);
const Material black(Vec(0, 0, 0));
const Material red(Vec(1, 0, 0));
const Material blue(Vec(0, 0, 1));
const Material green(Vec(0, 1, 0));
const Material yellow(Vec(1, 1, 0));
const Material purple(Vec(1, 0, 1));
const Material white(Vec(1, 1, 1));

// Scene
const Vec c(0, 0, 0);

static int clamp(double d){
  return std::max(0, std::min(255, (int) ceil(d)));
}

static double _gamma(double d){
  return pow(d, 0.45);
}

static Vec vgamma(Vec const& v){
  return Vec(_gamma(v.x()), _gamma(v.y()), _gamma(v.z()));
}

int main() {
  Scene scene;

  // scene.add_new_sphere(Sphere(c, 10, white));
  scene.add_new_sphere(Sphere(c - Vec(15, 0, 0), 10, water));
  scene.add_new_sphere(Sphere(c + Vec(15, 0, 0), 10, glass));
  scene.add_new_sphere(Sphere(c + Vec(-8, 8, 10), 3, mirror));

  scene.add_new_sphere(Sphere(Vec(0, 1000, 0), 940, red));
  scene.add_new_sphere(Sphere(Vec(0, 0, -1000), 940, green));
  scene.add_new_sphere(Sphere(Vec(-1000, 0, 0), 940, yellow));
  scene.add_new_sphere(Sphere(Vec(1000, 0, 0), 940, yellow));
  scene.add_new_sphere(Sphere(Vec(0, -1000, 0), 985, blue));
  
  // This one should be invisible unless there is a bug (or a reflexion)
  // scene.add_new_sphere(Sphere(Vec(0, 0, 1000), 940, purple));

  std::vector<float> image(W*H * 3, 0);

#pragma omp parallel for schedule(dynamic, H/8)
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      double x, y, z;
      x = j - W/2.0 + 0.5;
      y = - i + H/2 + 0.5;
      z = -W/(2*tan(fov/2.0));

      Vec dir = Vec(x, y, z).normalized();
      Ray r(origin, dir);

      Vec color(0, 0, 0);
      for(int i = 0; i < ray_number; i++){
        color = color + vgamma(scene.get_color(r, light));
      }

      color = color / (double) ray_number;

      image[(i * W + j) * 3 + 0] = color.r();
      image[(i * W + j) * 3 + 1] = color.g();
      image[(i * W + j) * 3 + 2] = color.b();
    }
  }

  // Convert float intensity to uchar pixel value
  std::vector<unsigned char> file_buffer(W*H * 3, 0);
#pragma omp parallel for schedule(dynamic, H/8)
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {
      file_buffer[(i * W + j) * 3 + 0] = clamp(image[(i * W + j) * 3 + 0]);
      file_buffer[(i * W + j) * 3 + 1] = clamp(image[(i * W + j) * 3 + 1]);
      file_buffer[(i * W + j) * 3 + 2] = clamp(image[(i * W + j) * 3 + 2]);
    }
  }

  stbi_write_png("image.png", W, H, 3, &file_buffer[0], 0);

  return 0;
}
