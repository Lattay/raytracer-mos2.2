#include "../src/Mesh.hpp"


int main(){
  RawMesh mesh("./misc/hexa.obj", Vec(0, 0, 0), 1);

  MeshBox box(mesh, 0, mesh.indices.size());

  Ray r(Vec(0, 0, 1), Vec(0, 1, -3));
  BoxIntersection inter = box.intersection(r);

  if(inter.valid()){
    std::cout << "size " << inter.size() << "\n";
    for(size_t i = 0; i < inter.size(); i++){
      std::cout << inter.slice(i).first << " " << inter.slice(i).length << "\n";
    }

    std::cout << "\n";
    for(auto it = inter.begin(); it != inter.end(); ++it){
      std::cout << *it << " ";
      Triangle tri(mesh, mesh.indices[*it]);
      std::cout << tri.center().y() << "\n";
      std::cout << tri.intersection(r).valid() << "\n";
    }
  } else {
    std::cout << "nope\n";
  }

  return 0;
}
