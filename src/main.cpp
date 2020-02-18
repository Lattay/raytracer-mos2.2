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

  std::cout << "Working on maximum of " << omp_get_max_threads() << " threads." << std::endl;

  Scene scene;

  // scene.add_new_sphere(Sphere(c, 10, white));
  scene.add_new_sphere(Sphere(c - Vec(15, 0, 0), 10, mirror));
  // scene.add_new_sphere(Sphere(c + Vec(15, 0, 40), 10, light_blue));
  scene.add_new_sphere(Sphere(c + Vec(-8, 8, 10), 3, purple));
  Mesh mesh = Mesh("./misc/suzanne.obj", 10, c);
  Vec box = mesh.box_size();
  std::cout << "Box size " << box.x() << ", " << box.y() << ", " << box.z() << std::endl;

  scene.add_mesh(mesh);

  scene.add_new_sphere(Sphere(Vec(0, 1000, 0), 940, white));
  scene.add_new_sphere(Sphere(Vec(0, 0, -1000), 940, green));
  scene.add_new_sphere(Sphere(Vec(-1000, 0, 0), 940, yellow));
  scene.add_new_sphere(Sphere(Vec(1000, 0, 0), 940, yellow));
  scene.add_new_sphere(Sphere(Vec(0, -1000, 0), 985, blue));

  std::cout << "Scene ready." << std::endl;
  
  // This one should be invisible unless there is a bug (or a reflexion)
  // scene.add_new_sphere(Sphere(Vec(0, 0, 1000), 940, purple));

  std::vector<float> image(W*H * 3, 0);

  for (int i = 0; i < H; i++) {
#pragma omp parallel for schedule(dynamic, H/8)
    for (int j = 0; j < W; j++) {

      double x, y, z;
      x = j - W/2.0 + 0.5;
      y = - i + H/2.0 + 0.5;
      z = -W/(2*tan(fov/2.0));

      Vec color(0, 0, 0);
      for(int i = 0; i < ray_number; i++){

        // Randomized direction to achieve anti-aliasing
        Vec2 v = box_muller(0.5);
        Vec2 v2 = box_muller(0.5);

        Vec delta_c = Vec(v2.x, v2.y, 0);
        Vec c_prim = origin + delta_c;

        Vec dir = Vec(x + v.x, y + v.y, z).normalized()*field_depth - delta_c;
        Ray r(c_prim, dir);

        color = color + scene.get_color(r, light);
      }

      color = vgamma(color / (double) ray_number);

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
