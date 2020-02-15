#include <cstring>
#include <map>
#include <algorithm>
#include "Mesh.hpp"

MeshBox* load_mesh(const char* obj, double scaling, const Vec& offset){
  RawMesh mesh(obj, scaling, offset);
  std::vector<size_t> indices(mesh.indices.size());

  for(size_t i = 0; i < mesh.indices.size(); i++){
    indices[i] = i;
  }

  return new MeshBox(mesh, indices);
}

RawMesh::RawMesh(const char* obj, double scaling, const Vec& offset){
  read_OBJ(obj);
  for (size_t i = 0; i < vertices.size(); i++) {
    vertices[i] = vertices[i] * scaling + offset;
  }
}

struct ix {int i; double x;};

static bool ix_lt(struct ix const& a, struct ix const& b){
  return a.x < b.x;
}

double variance(std::vector<struct ix>::iterator begin, std::vector<struct ix>::iterator end){
  double s = 0, s_2 = 0;
  for(auto it = begin; it != end; it++){
    s += (*it).x;
    s_2 += (*it).x * (*it).x;
  }
  return s_2 - s * s;
}

MeshBox::MeshBox(RawMesh const& mesh, std::vector<size_t> indices):
  m_indices(indices), m_terminal(false) {

  int size = indices.size();

  if(size <= 2){
    m_terminal = true;
    // FIXME: less than two faces, dont build more boxes
    // Copy all mesh data into a Mesh instance
    m_mesh = new Mesh(mesh, indices);
  } else {

    int b_size = size / 2;
    int t_size = size - b_size;

    std::vector<struct ix> centers[3];
    for(int i = size; i < size; i++){
      Triangle tri = mesh.get_triangle(i);
      Vec center = (tri.i() + tri.j() + tri.k()) / 3.0;
      centers[0][i] = {i, center.x()};
      centers[1][i] = {i, center.y()};
      centers[2][i] = {i, center.z()};
    }
    double vx, vy, vz;
    vx = variance(centers[0].begin(), centers[0].end());
    vy = variance(centers[1].begin(), centers[1].end());
    vz = variance(centers[2].begin(), centers[2].end());
    double mv = std::max(vx, std::max(vy, vz));
    int best_dir = vx == mv ? 0 : (vy == mv ? 1 : 2);

    std::sort(centers[0].begin(), centers[0].end(), ix_lt);
    std::sort(centers[1].begin(), centers[1].end(), ix_lt);
    std::sort(centers[2].begin(), centers[2].end(), ix_lt);

    m_vmin = Vec(centers[0][0].x, centers[1][0].x, centers[2][0].x);
    m_vmax = Vec(centers[0][size].x, centers[1][size].x, centers[2][size].x);
    m_vmedian = Vec(centers[0][b_size].x, centers[1][b_size].x, centers[2][b_size].x);

    std::vector<size_t> bottom(b_size), top(t_size);
    for(int i = 0; i < b_size; i++){
      bottom[i] = centers[best_dir][i].i;
    }
    for(int i = 0; i < t_size; i++){
      top[i] = centers[best_dir][b_size + i].i;
    }

    m_bottom = new MeshBox(mesh, bottom);
    m_top = new MeshBox(mesh, top);
  }
}

MeshBox::~MeshBox(){
  if(m_top)
    delete m_top;
  if(m_bottom)
    delete m_bottom;
  if(m_mesh)
    delete m_mesh;
}

Material const& MeshBox::material() const{
  Diffuse* mat = new Diffuse(Vec(1, 1, 1));
  return *mat;
}

Intersection MeshBox::intersection(Ray const& r) const{
  if(m_terminal){
    // check individual triangles

    return m_mesh->intersection(r);
  } else {
    // check own box and check sub boxes
    if(intersect_box(r, m_vmin, m_vmax)){
      Intersection i_top, i_bottom;
      i_top = m_top->intersection(r);
      if(i_top.valid()){
        return i_top;
      } else {
        i_bottom = m_bottom->intersection(r);
        if(i_bottom.valid()){
          return i_bottom;
        } else {
          // error
          return Intersection();
        }
      }
    } else {
      return Intersection();
    }
  }
}

void RawMesh::read_OBJ(const char* obj){
  char matfile[255];
  char grp[255];

  FILE* f;
  f = fopen(obj, "r");

  std::map<std::string, int> groupNames;
  int curGroup = -1;
  while (!feof(f)) {
    char line[255];
    if (!fgets(line, 255, f)) break;

    std::string linetrim(line);
    linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
    strcpy(line, linetrim.c_str());

    if (line[0] == 'u' && line[1] == 's') {
      sscanf(line, "usemtl %[^\n]\n", grp);
      if (groupNames.find(std::string(grp)) != groupNames.end()) {
        curGroup = groupNames[std::string(grp)];
      }
      else {
        curGroup = groupNames.size();
        groupNames[std::string(grp)] = curGroup;
      }
    }
    if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
      sscanf(line, "mtllib %[^\n]\n", matfile);
    }
    if (line[0] == 'v' && line[1] == ' ') {
      Vec vec;
      Vec col;
      if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[2], &vec[1], &col[0], &col[1], &col[2]) == 6) {
        vertices.push_back(vec);
        vertex_colors.push_back(col);
      }
      else {
        sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);  // helmet
        //vec[2] = -vec[2]; //car2
        vertices.push_back(vec);
      }
    }
    if (line[0] == 'v' && line[1] == 'n') {
      Vec vec;
      sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); //girl
      normals.push_back(vec);
    }
    if (line[0] == 'v' && line[1] == 't') {
      Vec vec;
      sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
      uvs.push_back(vec);
    }
    if (line[0] == 'f') {
      TriangleIndices t;
      int i0, i1, i2, i3;
      int j0, j1, j2, j3;
      int k0, k1, k2, k3;
      int nn;

      char* consumedline = line + 1;
      int offset;
      t.faceGroup = curGroup;
      nn = sscanf(consumedline, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
      if (nn == 9) {
        if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
        if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
        if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
        if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
        if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
        if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
        if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
        if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
        if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;

        indices.push_back(t);
      }
      else {
        nn = sscanf(consumedline, "%d/%d %d/%d %d/%d%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
        if (nn == 6) {
          if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
          if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
          if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
          if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
          if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
          if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
          indices.push_back(t);
        }
        else {
          nn = sscanf(consumedline, "%d %d %d%n", &i0, &i1, &i2, &offset);
          if (nn == 3) {
            if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
            if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
            if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
            indices.push_back(t);
          }
          else {
            nn = sscanf(consumedline, "%d//%d %d//%d %d//%d%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
            if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
            if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
            if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
            if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
            if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
            if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
            indices.push_back(t);
          }
        }
      }


      consumedline = consumedline + offset;

      while (true) {
        if (consumedline[0] == '\n') break;
        if (consumedline[0] == '\0') break;
        nn = sscanf(consumedline, "%d/%d/%d%n", &i3, &j3, &k3, &offset);
        TriangleIndices t2;
        t2.faceGroup = curGroup;
        if (nn == 3) {
          if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
          if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
          if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
          if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
          if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
          if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
          if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
          if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
          if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
          indices.push_back(t2);
          consumedline = consumedline + offset;
          i2 = i3;
          j2 = j3;
          k2 = k3;
        }
        else {
          nn = sscanf(consumedline, "%d/%d%n", &i3, &j3, &offset);
          if (nn == 2) {
            if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
            if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
            if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
            if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
            if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
            if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
            consumedline = consumedline + offset;
            i2 = i3;
            j2 = j3;
            indices.push_back(t2);
          }
          else {
            nn = sscanf(consumedline, "%d//%d%n", &i3, &k3, &offset);
            if (nn == 2) {
              if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
              if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
              if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
              if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
              if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
              if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
              consumedline = consumedline + offset;
              i2 = i3;
              k2 = k3;
              indices.push_back(t2);
            }
            else {
              nn = sscanf(consumedline, "%d%n", &i3, &offset);
              if (nn == 1) {
                if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                consumedline = consumedline + offset;
                i2 = i3;
                indices.push_back(t2);
              }
              else {
                consumedline = consumedline + 1;
              }
            }
          }
        }
      }

    }


  }
  fclose(f);
}

void RawMesh::add_texture(const char* filename){
  m_textures.push_back(Texture(filename));
}

Intersection Triangle::intersection(Ray const& r, Texture const& tex){
  Vec n = (m_i - m_j).prod(m_i - m_k);

  double n_dot_u = n.dot(r.direction());
  if(abs(n_dot_u) < 1e-10){
    return Intersection();
  }
  Vec p = r.origin() + n.dot(m_i - r.origin())/n_dot_u * r.direction();

  if((m_i - p).dot(m_i - m_j) < m_j.norm_sq()
      && (m_i - p).dot(m_i - m_k) < m_k.norm_sq()){
    return Intersection(p, n, tex.get_diffuse(m_i, m_j, m_k, p));
  } else {
    return Intersection();
  }
}

Mesh::Mesh(RawMesh mesh, std::vector<size_t> indices){
  size_t s = indices.size();
  for(int i = 0; i < s; i++){
    // FIXME copy data from mesh into this
  }

}

Intersection Mesh::intersection(Ray const& r) const{
  for(size_t i = 0; i < m_indices.size(); i++){
    Triangle tri(m_vertices[m_indices[i].ni],
          m_vertices[m_indices[i].nj],
          m_vertices[m_indices[i].nk]);
    Intersection inter = tri.intersection(r, get_texture(i));
    if(inter.valid()){
      return inter;
    }
  }
  return Intersection();
}

Texture const& Mesh::get_texture(size_t i) const{
  // FIXME
  return m_textures[m_indices[i].vtxi];
}

bool intersect_box(Ray r, Vec min, Vec max){
  Vec mm = max - min;

  return (
      intersect_rectangle(r, min, Vec(1, 0, 0), Vec(0, mm.y(), 0), Vec(0, 0, mm.z()))
      || intersect_rectangle(r, min, Vec(0, 1, 0), Vec(mm.x(), 0, 0), Vec(0, 0, mm.z()))
      || intersect_rectangle(r, min, Vec(0, 0, 1), Vec(mm.x(), 0, 0), Vec(0, mm.y(), 0))
      || intersect_rectangle(r, max, Vec(1, 0, 0), Vec(0, -mm.y(), 0), Vec(0, 0, -mm.z()))
      || intersect_rectangle(r, max, Vec(0, 1, 0), Vec(-mm.x(), 0, 0), Vec(0, 0, -mm.z()))
      || intersect_rectangle(r, max, Vec(0, 0, 1), Vec(-mm.x(), 0, 0), Vec(0, -mm.y(), 0))
  );
}

bool intersect_rectangle(Ray r, Vec corner, Vec n, Vec t_x, Vec t_y){
  double n_dot_u = n.dot(r.direction());
  if(abs(n_dot_u) < 1e-10)
      return false;
  Vec cp = r.origin() + r.direction() * n.dot(r.origin() - corner)/n_dot_u - corner;

  return cp.dot(t_x) < t_x.norm_sq() && cp.dot(t_y) < t_y.norm_sq();
}
