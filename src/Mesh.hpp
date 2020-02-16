#ifndef MESH_HPP
#define MESH_HPP
#include <vector>
#include "Object.hpp"
#include "Vec.hpp"
#include "Ray.hpp"
#include "Intersection.hpp"
#include "Material.hpp"
#include <vector>

class TriangleIndices {
  public:
    int vtxi, vtxj, vtxk;
    int uvi, uvj, uvk;
    int ni, nj, nk;
    int faceGroup;

    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1):
      vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) { };
};

class Triangle {
  private:
    Vec m_i, m_j, m_k;
  public:
    Triangle(Vec i, Vec j, Vec k):
      m_i(i), m_j(j), m_k(k) {};

    Vec operator[](int i){
      switch(i){
        case 0: return m_i;
        case 1: return m_j;
        case 2: return m_k;
      }
    }
    Vec const& i(){ return m_i; };
    Vec const& j(){ return m_j; };
    Vec const& k(){ return m_k; };

    Intersection intersection(Ray const& r, Texture const& tex);
};

class RawMesh {

  public:
    std::vector<TriangleIndices> indices;
    std::vector<Vec> vertices;
    std::vector<Vec> normals;
    std::vector<Vec> uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vec> vertex_colors;
    std::vector<Texture> m_textures;

    RawMesh() {};
    RawMesh(const char* obj, double scaling, const Vec& offset);

    void read_OBJ(const char* obj);
    void add_texture(const char* filename); 
    Triangle get_triangle(int index) const{
      return Triangle(vertices[indices[index].ni],
          vertices[indices[index].nj],
          vertices[indices[index].nk]);
    };
};

class Mesh {

  private:
    std::vector<TriangleIndices> m_indices;
    std::vector<Vec> m_vertices;
    std::vector<Vec> m_normals;
    std::vector<Vec> m_uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vec> m_vertex_colors;

    std::vector<Texture> m_textures;

  public:
    Mesh() {};
    Mesh(RawMesh parent, std::vector<size_t> indices);

    Intersection intersection(Ray const& r) const;

    Texture const& get_texture(size_t i) const;
};

class MeshBox final: public Object {
  protected:
    std::vector<size_t> m_indices;
    bool m_terminal;
    Vec m_vmin, m_vmax; //, m_vmedian;
    MeshBox const* m_top;
    MeshBox const* m_bottom;
    Mesh const* m_mesh;
  public:
    MeshBox() {};
    ~MeshBox();
    MeshBox(RawMesh const& mesh, std::vector<size_t> indices);

    virtual Intersection intersection(Ray const& r) const;
    Vec box_size() const { return m_vmax - m_vmin; };
};


MeshBox* load_mesh(const char* obj, double scaling, const Vec& offset);
bool intersect_box(Ray r, Vec min, Vec max);
bool intersect_rectangle(Ray r, Vec corner, Vec n, Vec t_x, Vec t_y);

#endif
