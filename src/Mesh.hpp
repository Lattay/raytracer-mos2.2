#ifndef MESH_HPP
#define MESH_HPP
#include <vector>
#include "Object.hpp"
#include "Vec.hpp"
#include "Ray.hpp"
#include "Intersection.hpp"
#include "Material.hpp"

class TriangleIndices {
  public:
    int vtxi, vtxj, vtxk;
    int uvi, uvj, uvk;
    int ni, nj, nk;
    int faceGroup;

    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1):
      vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) { };
};

typedef std::vector<size_t> Indices;

class TriIntersection final: public Intersection {
  private:
    bool m_valid;
    Vec m_position;
    Indices const* m_indices;

  public:
    TriIntersection():
      m_valid(false) {};
    TriIntersection(Vec position, Indices const* indices):
      m_valid(true), m_position(position), m_indices(indices) {};

    bool valid() const {return m_valid;};
    Vec const& position() const{return m_position;};
    Indices const& indices() const{return *m_indices;};
};

class RawMesh {

  public:
    std::vector<TriangleIndices> indices;
    std::vector<Vec> vertices;
    std::vector<Vec> normals;
    std::vector<Vec> uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vec> vertex_colors;
    std::vector<Texture> textures;

    RawMesh() {};
    RawMesh(const char* obj, double scaling, const Vec& offset);

    void read_OBJ(const char* obj);
    void add_texture(const char* filename); 
};

class MeshBox {
  protected:
    Indices m_indices;
    bool m_terminal;
    Vec m_vmin, m_vmax, m_center; //, m_vmedian;
    MeshBox const* m_top;
    MeshBox const* m_bottom;
  public:
    MeshBox() {};
    ~MeshBox();
    MeshBox(RawMesh const& mesh, Indices indices);

    TriIntersection intersection(Ray const& r) const;
    Vec box_size() const { return m_vmax - m_vmin; };
};

class Mesh final: public Object{

  private:
    RawMesh* m_mesh;
    MeshBox* m_box;

  public:
    Mesh() {};
    Mesh(const char* obj, double scaling, const Vec& offset);
    ~Mesh();


    Texture const& get_texture(size_t i) const;
    virtual Intersection intersection(Ray const& r) const;
    Vec box_size() const{return m_box->box_size();};
};

#endif
