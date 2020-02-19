#include "main.hpp" 

static Config gconf; 
 
static int clamp(double d){
  return std::max(0, std::min(255, (int) ceil(d)));
}

static double _gamma(double d){
  return pow(d, 0.45);
}

static Vec vgamma(Vec const& v){
  return Vec(_gamma(v.x()), _gamma(v.y()), _gamma(v.z()));
}

static void init_config(const char* file_name, Config* conf){
  std::ifstream file;

  file.open(file_name);
  char buffer[255];

  while(!file.eof()){
    if(file.fail()){
      std::cerr << "Something wrong happened with config." << std::endl;
      exit(1);
    }
    file.getline(buffer, 255);
    if(strncmp("ray_number ", buffer, 11) == 0){
      sscanf(buffer, "ray_number %d", &conf->ray_number);
      std::cout << "ray_number " << conf->ray_number << std::endl;
    } else if(strncmp("W ", buffer, 2) == 0){
      sscanf(buffer, "W %d", &conf->W);
      std::cout << "W " << conf->W << std::endl;
    } else if(strncmp("H ", buffer, 2) == 0){
      sscanf(buffer, "H %d", &conf->H);
      std::cout << "H " << conf->H << std::endl;
    } else if(strncmp("field_depth ", buffer, 12) == 0){
      sscanf(buffer, "field_depth %f", &conf->field_depth);
      std::cout << "field_depth " << conf->field_depth << std::endl;
    } else if(strncmp("focal_opening ", buffer, 14) == 0){
      sscanf(buffer, "focal_opening %f", &conf->focal_opening);
      std::cout << "focal_opening " << conf->focal_opening << std::endl;
    } else if(strncmp("fov ", buffer, 4) == 0){
      sscanf(buffer, "fov %f", &conf->fov);
      std::cout << "fov " << conf->fov << std::endl;
    } else {
      std::cout << buffer << std::endl;
    }
  }
  file.close();
}

int main() {

  std::cout << "Working on maximum of " << omp_get_max_threads() << " threads." << std::endl;

  Scene scene;

  Config gconf;

  init_config("./raytracer.cfg", &gconf);

  // scene.add_new_sphere(Sphere(c, 10, white));
  scene.add_new_sphere(Sphere(c - Vec(15, 0, 0), 10, mirror));
  scene.add_new_sphere(Sphere(c + Vec(15, 0, 40), 10, light_blue));
  scene.add_new_sphere(Sphere(c + Vec(-8, 8, 10), 3, purple));
  // Mesh mesh("./misc/cube.obj", 1, c);
  // scene.add_mesh(mesh);
  // Vec box = mesh.box_size();
  // std::cout << "Box size " << box.x() << ", " << box.y() << ", " << box.z() << std::endl;


  scene.add_new_sphere(Sphere(Vec(0, 1000, 0), 940, white));
  scene.add_new_sphere(Sphere(Vec(0, 0, -1000), 940, green));
  scene.add_new_sphere(Sphere(Vec(-1000, 0, 0), 940, yellow));
  scene.add_new_sphere(Sphere(Vec(1000, 0, 0), 940, yellow));
  scene.add_new_sphere(Sphere(Vec(0, -1000, 0), 985, blue));

  std::cout << "Scene ready." << std::endl;
  
  // This one should be invisible unless there is a bug (or a reflexion)
  // scene.add_new_sphere(Sphere(Vec(0, 0, 1000), 940, purple));

  std::vector<float> image(gconf.W * gconf.H * 3, 0);

  for (int i = 0; i < gconf.H; i++) {
#pragma omp parallel for schedule(dynamic, gconf.H/8)
    for (int j = 0; j < gconf.W; j++) {

      double x, y, z;
      x = j - gconf.W/2.0 + 0.5;
      y = - i + gconf.H/2.0 + 0.5;
      z = -gconf.W/(2*tan(gconf.fov/2.0));

      Vec color(0, 0, 0);
      for(int i = 0; i < gconf.ray_number; i++){

        // Randomized direction to achieve anti-aliasing
        Vec2 v = box_muller(0.5);
        Vec2 v2 = box_muller(0.5);

        double dx = (roll() - 0.5) * gconf.focal_opening;
        double dy = (roll() - 0.5) * gconf.focal_opening;

        Vec delta_c = Vec(dx, dy, 0);
        Vec c_prim = origin + delta_c;

        Vec dir = Vec(x + v.x, y + v.y, z).normalized() * gconf.field_depth - delta_c;
        Ray r(c_prim, dir);

        color = color + scene.get_color(r, light);
      }

      color = vgamma(color / (double) gconf.ray_number);

      image[(i * gconf.W + j) * 3 + 0] = color.r();
      image[(i * gconf.W + j) * 3 + 1] = color.g();
      image[(i * gconf.W + j) * 3 + 2] = color.b();
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
