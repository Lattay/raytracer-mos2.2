#include "../src/Mesh.hpp"


int main(){
  RawMesh mesh("./misc/cube.obj", Vec(0, 0, 0), 1);

  MeshBox box(mesh, 0, mesh.indices.size());

  BoxIntersection inter = box.intersection(Ray(Vec(0, 0, 20), Vec(-1, 0, -3)));

  std::cout << inter.valid() << "\n";

  return 0;
}
