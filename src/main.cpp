#include "main.hpp" 

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

  // Setup scene
  Config gconf; Scene scene;

  init_config("./raytracer.cfg", gconf, scene);

  Light light(Vec(-10, 20, 40), gconf.light_intensity, 10);

  Base cam_base = base_from(gconf.direction);

  std::cout << std::endl << "Current config:" << std::endl;
  show_config(gconf);

  std::cout << std::endl << "Scene ready." << std::endl;
  
  // Start rendering
  std::vector<float> image(gconf.W * gconf.H * 3, 0);

  std::cout << "Working on maximum of " << omp_get_max_threads() << " threads." << std::endl;

  for (int i = 0; i < gconf.H; i++) {
    #pragma omp parallel for schedule(dynamic, gconf.H/8)
    for (int j = 0; j < gconf.W; j++) {

      double x, y, z;
      x = gconf.W/2.0 - (j + 0.5);
      y = gconf.H/2.0 - (i + 0.5);
      z = gconf.W/(2*tan(gconf.fov/2.0));

      Vec color(0, 0, 0);
      for(int i = 0; i < gconf.ray_number; i++){
        // Randomized ray origin position for focal blur
        Vec2 v1 = box_muller(0.5);
        Vec delta_c = gconf.focal_opening * Vec(v1.x, v1.y, 0);
        Vec cam = gconf.camera + delta_c;

        // Randomized direction for anti-aliasing
        Vec2 v2 = box_muller(0.1 * gconf.antialiasing);
        Vec dir = gconf.field_depth * Vec(x + v2.x, y + v2.y, z).normalized() - delta_c;

        Ray r(cam, dir.rotate(cam_base));

        color = color + scene.get_color(r, light, gconf.recursive_depth);
      }

      color = vgamma(color / (double) gconf.ray_number);

      image[(i * gconf.W + j) * 3 + 0] = color.r();
      image[(i * gconf.W + j) * 3 + 1] = color.g();
      image[(i * gconf.W + j) * 3 + 2] = color.b();
    }
    if(i % (gconf.H/10) == 0 && i != 0){
      std::cout << "progress: " << 10 * i / (gconf.H / 10) << "%" << std::endl;
    }
  }

  // Convert float intensity to uchar pixel value
  std::vector<unsigned char> file_buffer(gconf.W*gconf.H * 3, 0);
  #pragma omp parallel for schedule(dynamic, gconf.H/8)
  for (int i = 0; i < gconf.H; i++) {
    for (int j = 0; j < gconf.W; j++) {
      file_buffer[(i * gconf.W + j) * 3 + 0] = clamp(image[(i * gconf.W + j) * 3 + 0]);
      file_buffer[(i * gconf.W + j) * 3 + 1] = clamp(image[(i * gconf.W + j) * 3 + 1]);
      file_buffer[(i * gconf.W + j) * 3 + 2] = clamp(image[(i * gconf.W + j) * 3 + 2]);
    }
  }

  stbi_write_png("image.png", gconf.W, gconf.H, 3, &file_buffer[0], 0);

  return 0;
}
