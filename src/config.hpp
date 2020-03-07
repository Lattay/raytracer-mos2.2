#ifndef CONFIG_HPP
#define CONFIG_HPP
#include "Scene.hpp"

// colors
const Transparent glass(1.5);
const Transparent water(1.2);
const Diffuse black(Vec(0, 0, 0));
const Diffuse red(Vec(1, 0, 0));
const Diffuse light_blue(Vec(0.5, 0.5, 1));
const Diffuse blue(Vec(0, 0, 1));
const Diffuse green(Vec(0, 1, 0));
const Diffuse yellow(Vec(1, 1, 0));
const Diffuse purple(Vec(1, 0, 1));
const Diffuse white(Vec(1, 1, 1));
const Reflective mirror(0.9);

typedef struct {
  int ray_number = 1;

  // Camera
  int W = 512;
  int H = 512;

  double field_depth = 55;
  double focal_opening = 0;

  double fov = pi/3.0;
  double antialiasing = 2;

  Vec camera = Vec(0, 0, 55);
  double light_intensity = 5e8;
} Config;

void init_config(const char* file_name, Config& conf, Scene& scene);

#endif
