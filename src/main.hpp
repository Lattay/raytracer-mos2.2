#ifndef MAIN_HPP
#define MAIN_HPP
#include <vector>
#include <iostream>
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../vendor/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../vendor/stb_image.h"

#include "Vec.hpp"
#include "Sphere.hpp"
#include "Ray.hpp"
#include "Scene.hpp"
#include "Light.hpp"
#include "random_tools.hpp"

const int ray_number = 2;

// Camera
const int W = 512;
const int H = 512;

const double field_depth = 55;
const double focal_opening = 0;

const double fov = pi/3.0;

const Vec origin(0, 0, 55);

// Light source
const Light light(Vec(-10, 20, 40), 5e8, 10);

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

// Scene
const Vec c(0, 0, 0);

#endif
