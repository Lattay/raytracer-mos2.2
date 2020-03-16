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

    TriangleIndices(int vtxi, int vtxj, int vtxk, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1):
      vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) { };
};

typedef std::vector<size_t> Indices;

class Triangle;

class RawMesh {

  public:
    std::vector<TriangleIndices> indices;
    std::vector<Vec> vertices;
    std::vector<Vec> normals;
    std::vector<Vec> uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vec> vertex_colors;
    std::vector<Texture> textures;

    RawMesh() {};
    RawMesh(const char* obj, const Vec& offset, double scaling);

    void read_OBJ(const char* obj);
    void add_texture(const char* filename); 
};

class Triangle final: public Object {
  private:
    Vec m_i, m_j, m_k, m_ni, m_nj, m_nk, m_uvi, m_uvj, m_uvk;
    // Texture m_tex;

  public:
    Triangle(RawMesh const& mesh, TriangleIndices idc):
      m_i(mesh.vertices[idc.vtxi]), m_j(mesh.vertices[idc.vtxj]), m_k(mesh.vertices[idc.vtxk]),
      m_ni(mesh.normals[idc.ni]), m_nj(mesh.normals[idc.nj]), m_nk(mesh.normals[idc.nk]),
      m_uvi(mesh.uvs[idc.uvi]), m_uvj(mesh.uvs[idc.uvj]), m_uvk(mesh.uvs[idc.uvk])
      /*, m_tex(mesh.textures[0]) */ {};
    virtual Intersection intersection(Ray const& r) const;
    Vec center() const{ return (m_i + m_j + m_k)/3; };
    Vec i() const{ return m_i; };
    Vec j() const{ return m_j; };
    Vec k() const{ return m_k; };
};

typedef struct {size_t first, length;} Slice;

class BoxIntersection {
  private:
    bool m_valid;
    std::vector<Slice> m_slices;

  public:
    class iterator: public std::iterator<std::input_iterator_tag, size_t, size_t, const size_t*, size_t >{
      size_t m_k;
      size_t m_slice_k;
      BoxIntersection const* m_inter;
      Slice m_slice;

      public:
        iterator(BoxIntersection const* inter):
          m_k(0), m_slice_k(0), m_inter(inter), m_slice(inter->slice(0)) {}
        iterator(BoxIntersection const* inter, size_t k, size_t slice_k):
          m_k(k), m_slice_k(slice_k), m_inter(inter), m_slice(inter->slice(0)) {}
        iterator& operator++() {
          if(m_k < m_slice.length - 1){
            m_k++;
          } else {
            m_k = 0;
            m_slice_k++;
            if(m_slice_k < m_inter->size()){
              m_slice = m_inter->slice(m_slice_k);
            }
          }
          return *this;
        }
        bool operator==(iterator other) const {
          return m_k == other.m_k && m_slice_k == other.m_slice_k;
        }
        bool operator!=(iterator other) const { return !(*this == other); }
        reference operator*() const {return m_k + m_slice.first;}
    };
    BoxIntersection():
      m_valid(false), m_slices(0) {};
    BoxIntersection(size_t start, size_t length):
      m_valid(true), m_slices(1, {start, length}) {};
    bool valid() const{ return m_valid; }
    size_t size() const { return m_slices.size(); };
    void append(Slice slice) { m_slices.push_back(slice); };
    Slice slice(size_t i) const{ return m_slices[i]; }
    iterator begin() const { return iterator(this); };
    iterator end() const { return iterator(this, 0, m_slices.size()); };
    std::vector<Slice> const& slices() const{ return m_slices; };
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
    BoundingBox& operator+=(BoundingBox const& b);
};

BoundingBox operator+(BoundingBox const& a, BoundingBox const& b);

class MeshBox {
  protected:
    bool m_terminal;
    BoundingBox m_bounding_box;
    int m_first = 0, m_length = 0;
    MeshBox const* m_top = NULL;
    MeshBox const* m_bottom = NULL;
  public:
    MeshBox(): m_terminal(true), m_bounding_box() {};
    ~MeshBox();
    MeshBox(RawMesh const& mesh, int first, int length);

    BoxIntersection intersection(Ray const& r) const;
    BoundingBox const& box() const{ return m_bounding_box; };
};

class Mesh final: public Object{

  private:
    RawMesh* m_mesh;
    MeshBox m_box;

  public:
    Mesh() {};
    Mesh(const char* obj, const Vec& offset, double scaling);
    ~Mesh();


    Texture const& get_texture(size_t i) const;
    virtual Intersection intersection(Ray const& r) const;
    BoundingBox const& box() const{ return m_box.box(); };
};

#endif
