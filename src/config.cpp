#include "config.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
#include "Vec.hpp"
#include "Sphere.hpp"
#include "Mesh.hpp"

void init_config(const char* file_name, Config& conf, Scene& scene){
  std::ifstream file;

  file.open(file_name);
  char buffer[255];

  std::cout << "Loading config file " << file_name << "..." << std::endl;
  while(!file.eof()){
    if(file.fail()){
      std::cerr << "Something wrong happened with config." << std::endl;
      exit(1);
    }
    file.getline(buffer, 255);
    if(buffer[0] == '#'
        || buffer[0] == 0){
      // comment or empty line
    } else if(strncmp("ray_number ", buffer, 11) == 0){
      sscanf(buffer, "ray_number %d", &conf.ray_number);
    } else if(strncmp("W ", buffer, 2) == 0){
      sscanf(buffer, "W %d", &conf.W);
    } else if(strncmp("H ", buffer, 2) == 0){
      sscanf(buffer, "H %d", &conf.H);
    } else if(strncmp("field_depth ", buffer, 12) == 0){
      sscanf(buffer, "field_depth %lf", &conf.field_depth);
    } else if(strncmp("focal_opening ", buffer, 14) == 0){
      sscanf(buffer, "focal_opening %lf", &conf.focal_opening);
    } else if(strncmp("fov ", buffer, 4) == 0){
      sscanf(buffer, "fov %lf", &conf.fov);
    } else if(strncmp("light_intensity ", buffer, 16) == 0){
      sscanf(buffer, "light_intensity %lf", &conf.light_intensity);
    } else if(strncmp("recursive_depth ", buffer, 16) == 0){
      sscanf(buffer, "recursive_depth %d", &conf.recursive_depth);
    } else if(strncmp("antialiasing ", buffer, 13) == 0){
      sscanf(buffer, "antialiasing %lf", &conf.antialiasing);
    } else if(strncmp("camera ", buffer, 7) == 0){
      sscanf(buffer, "camera %lf %lf %lf %lf %lf %lf",
          &conf.camera[0], &conf.camera[1], &conf.camera[2],
          &conf.direction[0], &conf.direction[1], &conf.direction[2]);
    } else if(strncmp("sphere ", buffer, 7) == 0){
      Vec offset;
      double radius;
      int mat_num;
      sscanf(buffer, "sphere %d %lf %lf %lf %lf", &mat_num, &radius, &offset[0], &offset[1], &offset[2]);
      switch(mat_num){
        case 0:
          scene.add_new_sphere(Sphere(offset, radius, white));
          break;
        case 1:
          scene.add_new_sphere(Sphere(offset, radius, blue));
          break;
        case 2:
          scene.add_new_sphere(Sphere(offset, radius, light_blue));
          break;
        case 3:
          scene.add_new_sphere(Sphere(offset, radius, green));
          break;
        case 4:
          scene.add_new_sphere(Sphere(offset, radius, yellow));
          break;
        case 5:
          scene.add_new_sphere(Sphere(offset, radius, red));
          break;
        case 6:
          scene.add_new_sphere(Sphere(offset, radius, mirror));
          break;
        case 7:
          scene.add_new_sphere(Sphere(offset, radius, glass));
          break;
        default:
          scene.add_new_sphere(Sphere(offset, radius, purple));
      }
    } else if(strncmp("mesh ", buffer, 5) == 0){
      char obj_file[255];
      Vec offset;
      double scaling;
      sscanf(buffer, "mesh %s %lf %lf %lf %lf", obj_file, &scaling, &offset[0], &offset[1], &offset[2]);
      scene.add_mesh(*(new Mesh(obj_file, offset, scaling)));
    } else {
      std::cout << "unexpected line: " << buffer << std::endl;
    }
  }
  std::cout << "Closing config file...";
  file.close();
  std::cout << " done." << std::endl;
}

void show_config(Config& conf){
  std::cout << "ray_number " << conf.ray_number << std::endl;
  std::cout << "W " << conf.W << std::endl;
  std::cout << "H " << conf.H << std::endl;
  std::cout << "field_depth " << conf.field_depth << std::endl;
  std::cout << "focal_opening " << conf.focal_opening << std::endl;
  std::cout << "fov " << conf.fov << std::endl;
  std::cout << "antialiasing " << conf.antialiasing << std::endl;
  std::cout << "recursive_depth " << conf.recursive_depth << std::endl;
  std::cout << "camera (" << conf.camera.x() << ", " << conf.camera.y() << ", " << conf.camera.z() << ")" << std::endl;
  std::cout << "light_intensity " << conf.light_intensity << std::endl;
}
