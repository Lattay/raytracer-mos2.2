#ifndef MESH_HPP
#define MESH_HPP
#include <vector>
#include "Object.hpp"
#include "Vec.hpp"
#include "Ray.hpp"
#include "Intersection.hpp"
#include "Material.hpp"
#include "Texture.hpp"

class TriangleIndices {
  public:
    int vtxi, vtxj, vtxk;
    int uvi, uvj, uvk;
    int ni, nj, nk;
    int faceGroup = -1;

    TriangleIndices(int vtxi, int vtxj, int vtxk, int ni, int nj, int nk, int uvi, int uvj, int uvk):
      vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) { };
};

class Triangle;

class RawMesh {

  public:
    std::vector<TriangleIndices> indices;
    std::vector<Vec> vertices;
    std::vector<Vec> normals;
    std::vector<Vec> uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vec> vertex_colors;
    std::vector<Texture*> textures;

    RawMesh() {};
    RawMesh(const char* obj, const Vec& offset, double scaling);

    void read_OBJ(const char* obj);
    void add_texture(Texture* text);
    Texture const* get_texture(size_t i) const;
};

class Triangle final: public Object {
  private:
    Vec m_i, m_j, m_k, m_ni, m_nj, m_nk, m_uvi, m_uvj, m_uvk;
    Texture const* m_tex = nullptr;

  public:
    Triangle() {};
    Triangle(RawMesh const& mesh, TriangleIndices idc):
      m_i(mesh.vertices.at(idc.vtxi)), m_j(mesh.vertices.at(idc.vtxj)), m_k(mesh.vertices.at(idc.vtxk)),
      m_ni(mesh.normals.at(idc.ni)), m_nj(mesh.normals.at(idc.nj)), m_nk(mesh.normals.at(idc.nk)),
      m_uvi(mesh.uvs.at(idc.uvi)), m_uvj(mesh.uvs.at(idc.uvj)), m_uvk(mesh.uvs.at(idc.uvk)),
      m_tex(mesh.get_texture(idc.faceGroup)) {};
    virtual Intersection intersection(Ray const& r) const;
    Vec center() const{ return (m_i + m_j + m_k)/3; };
    Vec i() const{ return m_i; };
    Vec j() const{ return m_j; };
    Vec k() const{ return m_k; };
    Diffuse* get_color(Vec pos) const;
};

typedef struct {size_t first, length;} Slice;

class BoxIntersection {
  private:
    std::vector<size_t> m_indices;

  public:
    BoxIntersection():
      m_indices(0) {};
    BoxIntersection(std::vector<size_t> const& indices):
      m_indices(indices.begin(), indices.end()) {};

    bool empty() const{ return m_indices.empty(); }
    size_t size() const { return m_indices.size(); };
    void append(size_t i) { m_indices.push_back(i); };
    size_t get(size_t i) const{ return m_indices[i]; }
    std::vector<size_t>::const_iterator begin() const{ return m_indices.begin(); };
    std::vector<size_t>::const_iterator end() const{ return m_indices.end(); };
};

BoxIntersection operator&&(BoxIntersection const& a, BoxIntersection const& b);

class BoundingBox {
  private:
    Vec m_vmin, m_vmax;
  public:
    BoundingBox():
      m_vmin(), m_vmax() {};
    BoundingBox(Vec min, Vec max):
      m_vmin(min), m_vmax(max) {};
    bool intersect(Ray const& r) const;
    Vec min() const { return m_vmin; };
    Vec max() const { return m_vmax; };
    double vol() const {
      Vec d = m_vmax - m_vmin; return d.x()*d.y()*d.z();
    };
    BoundingBox& operator+=(BoundingBox const& b);
};

BoundingBox operator+(BoundingBox const& a, BoundingBox const& b);
BoundingBox operator*(BoundingBox const& a, BoundingBox const& b);

class MeshBox {
  protected:
    bool m_terminal;
    bool m_empty;
    BoundingBox m_bounding_box;
    std::vector<size_t> m_indices;
    MeshBox const* m_top = nullptr;
    MeshBox const* m_bottom = nullptr;
  public:
    MeshBox(): m_terminal(true), m_empty(true), m_bounding_box(), m_indices(0) {};
    ~MeshBox();
    MeshBox(RawMesh const& mesh, std::vector<size_t> const& indices);

    BoxIntersection intersection(Ray const& r) const;
    BoundingBox const& box() const{ return m_bounding_box; };
};

class Mesh final: public Object{

  private:
    RawMesh* m_mesh = nullptr;
    MeshBox* m_box = nullptr;

  public:
    // Mesh() {};
    Mesh(const char* obj, const Vec& offset, double scaling);
    ~Mesh();


    void load_texture(const char* texture);
    Texture const* get_texture(size_t i) const;
    virtual Intersection intersection(Ray const& r) const;
    BoundingBox const& box() const{ return m_box->box(); };
};

#endif
