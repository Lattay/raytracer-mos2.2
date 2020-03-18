#include <cstring>
#include <map>
#include <algorithm>
#include <iostream>
#include "Mesh.hpp"
#include "Sphere.hpp"

static RawMesh* current_mesh;  // Thread unsafe !
static bool lt_tri(TriangleIndices a, TriangleIndices b){
  Triangle tri1(*current_mesh, a);
  Triangle tri2(*current_mesh, b);

  return tri1.center().y() < tri2.center().y();
}

Mesh::Mesh(const char* obj, const Vec& offset, double scaling){
  std::cout << "Loading mesh from " << obj
    << " at (" << offset.x() << ", " << offset.y() << ", " << offset.z() << ")... "  << std::endl;
  m_mesh = new RawMesh(obj, offset, scaling);

  std::cout << "Sort tris..." << std::endl;
  current_mesh = m_mesh;
  std::sort(m_mesh->indices.begin(), m_mesh->indices.end(), lt_tri);

  std::cout << "Building mesh box tree..." << std::endl;
  m_box = MeshBox(*m_mesh, 0, m_mesh->indices.size());

  std::cout << "Box: ("
    << m_box.box().min().x() << ", "
    << m_box.box().min().y() << ", "
    << m_box.box().min().z() << ") -> ("
    << m_box.box().max().x() << ", "
    << m_box.box().max().y() << ", "
    << m_box.box().max().z() << ")" << std::endl;


  std::cout << "Loading mesh done." << std::endl;
}

Mesh::~Mesh(){
  delete m_mesh;
}

RawMesh::RawMesh(const char* obj, const Vec& offset, double scaling){
  read_OBJ(obj);
  for(size_t i = 0; i < vertices.size(); i++) {
    vertices[i] = (scaling * vertices[i]) + offset;
  }
}

/*
typedef struct {int i; double x;} LabeledX;

static bool lt(LabeledX const& a, LabeledX const& b){
  return a.x < b.x;
}
*/

BoxIntersection operator&&(BoxIntersection const& a, BoxIntersection const& b){
  if(a.empty()) {
    return b;
  } else if(b.empty()){
    return a;
  } else {
    BoxIntersection sum;
    for(auto it = a.begin(); it != a.end(); ++it){
      sum.append(*it);
    }
    for(auto it = b.begin(); it != b.end(); ++it){
      sum.append(*it);
    }
    return sum;
  }
}

static const Vec veps(0.1, 0.1, 0.1);

MeshBox::MeshBox(RawMesh const& mesh, int first, int length):
  m_terminal(true), m_bounding_box(), m_first(first), m_length(length) {

  if(length <= 0){
    throw "Length should not be null.";
  } else if(length <= 2){
    Triangle tri(mesh, mesh.indices[m_first]);
    m_bounding_box = BoundingBox(tri.i() - veps, tri.i() + veps)
      + BoundingBox(tri.j() - veps, tri.j() + veps)
      + BoundingBox(tri.k() - veps, tri.k() + veps);
    if(length == 2){
      Triangle tri2(mesh, mesh.indices[m_first + 1]);
      m_bounding_box  += BoundingBox(tri2.i() - veps, tri2.i() + veps)
        + BoundingBox(tri2.j() - veps, tri2.j() + veps)
        + BoundingBox(tri2.k() - veps, tri2.k() + veps);
    }
  } else {
    m_terminal = false;
    size_t first_half = length/2;
    size_t second_half = length - first_half;
    m_bottom = new MeshBox(mesh, first, first_half);
    m_top = new MeshBox(mesh, first + first_half, second_half);

    m_bounding_box = m_top->m_bounding_box + m_bottom->m_bounding_box;

    /*
    std::vector<LabeledX> verts[3] = {
      std::vector<LabeledX>(length),
      std::vector<LabeledX>(length),
      std::vector<LabeledX>(length),
    };

    for(int i = 0; i < length; i++){
      Vec v = Triangle(mesh, i+first);
      verts[0][i] = {i, v.x()};
      verts[1][i] = {i, v.y()};
      verts[2][i] = {i, v.z()};
    }

    std::sort(verts[0].begin(), verts[0].end(), lt);
    std::sort(verts[1].begin(), verts[1].end(), lt);
    std::sort(verts[2].begin(), verts[2].end(), lt);

    m_vmin = Vec(verts[0][0].x, verts[1][0].x, verts[2][0].x);
    m_vmax = Vec(verts[0][length-1].x, verts[1][length-1].x, verts[2][length-1].x);
    */

    /*
    double vx, vy, vz;
    vx = (*verts[0].end()).x - (*verts[0].begin()).x;
    vy = (*verts[1].end()).x - (*verts[1].begin()).x;
    vz = (*verts[2].end()).x - (*verts[2].begin()).x;
    double mv = std::max(vx, std::max(vy, vz));

    int best_dir = vx == mv ? 0 : (vy == mv ? 1 : 2);

    Indices bottom(b_size), top(t_size);
    for(int i = 0; i < b_size; i++){
      bottom[i] = verts[best_dir][i].i;
    }
    for(int i = 0; i < t_size; i++){
      top[i] = verts[best_dir][b_size + i].i;
    }

    m_bottom = new MeshBox(mesh, bottom);
    m_top = new MeshBox(mesh, top);
    */
  }
}

MeshBox::~MeshBox(){
  if(!m_terminal){
    delete m_top;
    delete m_bottom;
  }
}

BoxIntersection MeshBox::intersection(Ray const& r) const{
  if(m_terminal){  // checking two triangles is faster than checking six rectangles
    return BoxIntersection(m_first, m_length);
  } else {
    if(m_bounding_box.intersect(r)){
      return m_bottom->intersection(r) && m_top->intersection(r);
    } else {
      return BoxIntersection();
    }
  }
}

static inline double abs(double x){ return (x > 0.0) ? x : -x; }

/* Return true if ray r intersect the rectangle defined in the plane (corner, n)
 * by the point corner and the side vectors t_x and t_y.
 */
static bool intersect_rectangle(Ray const& r, Vec corner, Vec n, Vec t_x, Vec t_y){
  double n_dot_u = n.dot(r.direction());
  if(abs(n_dot_u) < 1e-30)  // u is tangeant to the plane
      return false;
  double alpha = n.dot(corner - r.origin())/n_dot_u;
  if(alpha < 0)
    return false;
  // cp is the corner to intersection vector
  Vec cp = r.origin() + r.direction() * alpha - corner;

  // u and v are the coordinates in the local referential (times the respective magnitude of t_x and t_y)
  double u = cp.dot(t_x), v = cp.dot(t_y);
  return (
      (u >= 0) && (u <= t_x.norm_sq())
      && (v >= 0) && (v <= t_y.norm_sq())
  );
}

bool BoundingBox::intersect(Ray const& r) const{
  Vec mm = m_vmax - m_vmin;

  return (
      intersect_rectangle(r, m_vmin, Vec(1, 0, 0), Vec(0, mm.y(), 0), Vec(0, 0, mm.z()))
      || intersect_rectangle(r, m_vmin, Vec(0, 1, 0), Vec(mm.x(), 0, 0), Vec(0, 0, mm.z()))
      || intersect_rectangle(r, m_vmin, Vec(0, 0, 1), Vec(mm.x(), 0, 0), Vec(0, mm.y(), 0))
      || intersect_rectangle(r, m_vmax, Vec(1, 0, 0), Vec(0, -mm.y(), 0), Vec(0, 0, -mm.z()))
      || intersect_rectangle(r, m_vmax, Vec(0, 1, 0), Vec(-mm.x(), 0, 0), Vec(0, 0, -mm.z()))
      || intersect_rectangle(r, m_vmax, Vec(0, 0, 1), Vec(-mm.x(), 0, 0), Vec(0, -mm.y(), 0))
  );
}

BoundingBox operator+(BoundingBox const& a, BoundingBox const& b){
  Vec min, max;
  for(int i = 0; i < 3; i++){
    min[i] = std::min(a.min()[i], b.min()[i]);
    max[i] = std::max(a.max()[i], b.max()[i]);
  }
  return BoundingBox(min, max);
}

BoundingBox& BoundingBox::operator+=(BoundingBox const& b){
  for(int i = 0; i < 3; i++){
    m_vmin[i] = std::min(m_vmin[i], b.min()[i]);
    m_vmax[i] = std::max(m_vmax[i], b.max()[i]);
  }
  return *this;
}

void RawMesh::read_OBJ(const char* obj){
  char matfile[255];
  char grp[255];

  FILE* f;
  f = fopen(obj, "r");
  if(!f){
    std::cerr << "Cannot open " << obj << std::endl;
    exit(1);
  }

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
      } else {
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
      if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
        vertices.push_back(vec);
        vertex_colors.push_back(col);
      }
      else {
        sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);  // helmet
        //vec[2] = -vec[2]; //car2
        vertices.push_back(vec);
      }
    }
    if (line[0] == 'v' && line[1] == 'n') {
      Vec vec;
      sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]); //girl
      normals.push_back(vec);
    }
    if (line[0] == 'v' && line[1] == 't') {
      Vec vec;
      sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
      uvs.push_back(vec);
    }
    if (line[0] == 'f') {
      TriangleIndices t(0, 0, 0, 0, 0, 0, 0, 0, 0);
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
        TriangleIndices t2(0, 0, 0, 0, 0, 0, 0, 0, 0);
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
  textures.push_back(Texture(filename));
}

Intersection Mesh::intersection(Ray const& r) const{

  double d_2, min_d_2 = 1e30;
  
  BoxIntersection box_inter(m_box.intersection(r));

  Intersection closest;  // invalid by default

  if(!box_inter.empty()){
    for(auto it = box_inter.begin(); it != box_inter.end(); ++it){
      Triangle tri(*m_mesh, m_mesh->indices[*it]);

      Intersection inter = tri.intersection(r);
      if(inter.valid()){
        d_2 = (inter.position() - r.origin()).norm_sq();
        if(d_2 < min_d_2){
          closest = inter;
          min_d_2 = d_2;
        }
      }
    }
  }

  return closest;
}

Texture const& Mesh::get_texture(size_t i) const{
  return m_mesh->textures[m_mesh->indices[i].vtxi];
}

Intersection Triangle::intersection(Ray const& r) const{
  Vec vk = m_k - m_i;
  Vec vj = m_j - m_i;
  Vec n = m_ni + m_nk + m_nj;

  double n_dot_u = n.dot(r.direction());
  if(n_dot_u < -1e-30){
    double n_dot_oi = n.dot(m_i - r.origin());
    if(n_dot_oi < 0.0){
      Vec vp = r.origin() + n_dot_oi/n_dot_u * r.direction() - m_i;

      double k_dot_k = vk.norm_sq();
      double k_dot_j = vk.dot(vj);
      double k_dot_p = vk.dot(vp);
      double j_dot_j = vj.norm_sq();
      double j_dot_p = vj.dot(vp);

      double inv_denom = 1/(k_dot_k * j_dot_j - k_dot_j * k_dot_j);

      // u, v and w (= 1 - u - v) are the barycentric coordinates
      double u = (j_dot_j * k_dot_p - k_dot_j * j_dot_p) * inv_denom;
      double v = (k_dot_k * j_dot_p - k_dot_j * k_dot_p) * inv_denom;

      if(u >= 0 && v >= 0 && (u + v) < 1){
        Vec n_rend(u * m_ni + v * m_nj + (1 - u - v) * m_nk);
        return Intersection(vp + m_i, n_rend, white_emit /*m_tex.get_diffuse(m_i, m_j, m_k, vp + m_i)*/);
      } else {  // intersection is out of the triangle
        return Intersection();
      }
    } else {  // intersection is behind the camera
      return Intersection();
    }
  } else {  // ray intersect the face from the back or is almost tangeant to it
    return Intersection();
  }
}

