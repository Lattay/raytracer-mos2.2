#include <cstring>
#include <map>
#include <algorithm>
#include <iostream>
#include "../vendor/stb_image.h"
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
  size_t length = m_mesh->indices.size();
  std::vector<size_t> v(length);
  for(size_t i = 0; i < length; i++)
    v[i] = i;
  m_box = new MeshBox(*m_mesh, v);

  std::cout << "Box: ("
    << m_box->box().min().x() << ", "
    << m_box->box().min().y() << ", "
    << m_box->box().min().z() << ") -> ("
    << m_box->box().max().x() << ", "
    << m_box->box().max().y() << ", "
    << m_box->box().max().z() << ")" << std::endl;


  std::cout << "Loading mesh done." << std::endl;
}

Mesh::~Mesh(){
  delete m_mesh;
  delete m_box;
}

RawMesh::RawMesh(const char* obj, const Vec& offset, double scaling){
  read_OBJ(obj);
  for(size_t i = 0; i < vertices.size(); i++) {
    vertices[i] = (scaling * vertices[i]) + offset;
  }
}

BoxIntersection operator&&(BoxIntersection const& a, BoxIntersection const& b){
  BoxIntersection sum;
  if(!a.empty()){
    for(auto it = a.begin(); it != a.end(); ++it){
      sum.append(*it);
    }
  }
  if(!b.empty()){
    for(auto it = b.begin(); it != b.end(); ++it){
      sum.append(*it);
    }
  }
  return sum;
}

static const Vec veps(1e-5, 1e-5, 1e-5);

MeshBox::MeshBox(RawMesh const& mesh, std::vector<size_t> const& indices):
  m_terminal(true), m_bounding_box(), m_indices(indices) {

  size_t length = m_indices.size();

  if(length <= 0){
    // empty indices
  } else if(length <= 2){
    Triangle tri(mesh, mesh.indices[m_indices[0]]);
    m_bounding_box = BoundingBox(tri.i() - veps, tri.i() + veps)
      + BoundingBox(tri.j() - veps, tri.j() + veps)
      + BoundingBox(tri.k() - veps, tri.k() + veps);
    if(length == 2){
      Triangle tri2(mesh, mesh.indices[m_indices[1]]);
      m_bounding_box  += BoundingBox(tri2.i() - veps, tri2.i() + veps)
        + BoundingBox(tri2.j() - veps, tri2.j() + veps)
        + BoundingBox(tri2.k() - veps, tri2.k() + veps);
    }
  } else {
    m_terminal = false;

    std::vector<Triangle> tris(length);
    Vec center(0, 0, 0);
    for(size_t i = 0; i < length; i++){
      tris[i] = Triangle(mesh, mesh.indices[m_indices[i]]);
      center = center + tris[i].center()/length;
    }

    Vec deviation;

    for(size_t i = 0; i < length; i++){
      Vec dev = tris[i].center()/length - center;
      deviation[0] += dev[0]*dev[0]/length;
      deviation[1] += dev[1]*dev[1]/length;
      deviation[2] += dev[2]*dev[2]/length;
    }

    double mv = std::max(deviation.x(), std::max(deviation.y(), deviation.z()));
    int best_dir;
    Vec median_normal;
    if(mv == deviation.x()){
      best_dir = 0;
      median_normal = Vec(1, 0, 0);
    } else if(mv == deviation.y()){
      best_dir = 1;
      median_normal = Vec(0, 1, 0);
    } else {
      best_dir = 2;
      median_normal = Vec(0, 0, 1);
    }

    size_t m = length/2;

    // nth_element partitionate m_indices such that the mth element is greater
    // than any before and lower than any after in term of best_dir coordinate.
    nth_element(m_indices.begin(), m_indices.begin()+m, m_indices.end(),
        [mesh, best_dir](size_t i, size_t j)->bool{
          return Triangle(mesh, mesh.indices[i]).center()[best_dir]
            < Triangle(mesh, mesh.indices[j]).center()[best_dir];
        }
    );

    // then 
    std::vector<size_t> bottom(m_indices.begin(), m_indices.begin()+m);
    std::vector<size_t> top(m_indices.begin()+m, m_indices.end());

    m_bottom = new MeshBox(mesh, bottom);
    m_top = new MeshBox(mesh, top);

    m_bounding_box = m_top->m_bounding_box + m_bottom->m_bounding_box;

    // double superpos = (m_top->m_bounding_box * m_bottom->m_bounding_box).vol();
    // std::cout << "Superposition: " << 2*superpos/(m_top->m_bounding_box.vol() + m_bottom->m_bounding_box.vol())*100 << "% "
    //   << 0.5*(m_top->m_bounding_box.vol() + m_bottom->m_bounding_box.vol()) << "\n";
  }
}

MeshBox::~MeshBox(){
  if(!m_terminal){
    delete m_top;
    delete m_bottom;
  }
}

BoxIntersection MeshBox::intersection(Ray const& r) const{
  if(m_indices.empty()){
    return BoxIntersection();
  } else if(m_terminal){  // checking two triangles is faster than checking six rectangles
      return BoxIntersection(m_indices);
  } else if(m_bounding_box.intersect(r)){
    return m_bottom->intersection(r) && m_top->intersection(r);
  } else {
    return BoxIntersection();
  }
}

static inline double dabs(double x){ return (x > 0.0) ? x : -x; }

/* Return true if ray r intersect the rectangle defined in the plane (corner, n)
 * by the point corner and the side vectors t_x and t_y.
 */
static bool intersect_rectangle(Ray const& r, Vec corner, Vec n, Vec t_x, Vec t_y){
  double n_dot_u = n.dot(r.direction());
  if(dabs(n_dot_u) < 1e-30)  // u is tangeant to the plane
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

BoundingBox operator*(BoundingBox const& a, BoundingBox const& b){
  Vec min, max;
  for(int i = 0; i < 3; i++){
    min[i] = std::max(a.min()[i], b.min()[i]);
    max[i] = std::min(a.max()[i], b.max()[i]);

    if(min[i] > max[i]){
      min[i] = max[i];
    }
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

void RawMesh::add_texture(Texture* tex){
  textures.push_back(tex);
}

Texture const* RawMesh::get_texture(size_t i) const{
  if(i < textures.size())
    return textures[i];
  else
    return nullptr;
}

Intersection Mesh::intersection(Ray const& r) const{

  double d_2, min_d_2 = 1e30;
  
  BoxIntersection box_inter(m_box->intersection(r));

  Intersection closest;  // invalid by default

  if(!box_inter.empty()){
    int i = 0;
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
      i++;
    }
    // if(i > m){
    //   m = i;
    //   std::cout << "up " << m << "\n";
    // }

  }

  return closest;
}

void Mesh::load_texture(const char* texture){
  m_mesh->add_texture(new Texture(texture));
}

Texture const* Mesh::get_texture(size_t i) const{
  return m_mesh->get_texture(i);
}

Intersection Triangle::intersection(Ray const& r) const{
  Vec vk = m_k - m_i;
  Vec vj = m_j - m_i;
  Vec n = vj.prod(vk);

  if(n.dot(m_ni + m_nj + m_nk) < 0){
    n = -n;
  }

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

      // u, v and t (= 1 - u - v) are the barycentric coordinates
      // t + u + v = 1
      // and (all equivalent):
      // 0 = t PI + u PJ + v PK
      // P = t I + u J + v K
      // P = I + u IJ + v IK
      double u = (k_dot_k * j_dot_p - k_dot_j * k_dot_p) * inv_denom;
      double v = (j_dot_j * k_dot_p - k_dot_j * j_dot_p) * inv_denom;

      if(u >= 0 && v >= 0 && (u + v) < 1){
        Vec n_rend((1 - u - v) * m_ni + u * m_nj + v * m_nk);
        Vec vtex((1 - u - v) * m_uvi + u * m_uvj + v * m_uvk);
        return Intersection(vp + m_i, n_rend, get_color(vtex));
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

Diffuse* Triangle::get_color(Vec pos) const{
  if(m_tex != nullptr){
    return new Diffuse(m_tex->get_color(pos));
  } else {
    return new Diffuse(white_emit);
  }
}
