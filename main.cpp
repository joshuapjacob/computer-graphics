#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <cstring>
#include <stdio.h>
#include <algorithm>
#include <list>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "vector.h"

// INTERSECTION ----------------------------------------------------------------

struct Intersection {
  bool intersected = false;
  bool reflective = false;
  bool is_light = false;
  double refractive_index = 1.;
  double t;
  Vector sphere_center;
  Vector P;
  Vector N;
  Vector albedo;
};

// BOUNDING BOX ----------------------------------------------------------------

struct BoundingBox {
  Vector B_min;
  Vector B_max;
};

// RAY -------------------------------------------------------------------------

class Ray {
  public:
    Ray(const Vector& O, const Vector& u) {
        this->O = O;
        this->u = u;
    }
    Vector O;
    Vector u;
};

// GEOMETRY --------------------------------------------------------------------

class Geometry {
  public:
    virtual Intersection intersect(const Ray &ray) = 0;
    bool reflective;
};

// SPHERE ----------------------------------------------------------------------

class Sphere : public Geometry {
  public:
    Sphere(
      const Vector& C,
      const double& R,
      const Vector& albedo,
      bool reflective = false,
      double refractive_index = 1.,
      bool invert_normal = false,
      bool is_light_source = false
    ) {
        this->C = C;
        this->R = R;
        this->albedo = albedo;
        this->reflective = reflective;
        this->refractive_index = refractive_index;
        this->invert_normal = invert_normal;
        this->is_light_source = is_light_source;
    }
    Intersection intersect(const Ray &ray) override {
      Intersection intersection;
      Vector O_C = ray.O - C;
      double discriminant = pow(dot(ray.u, O_C), 2) - (pow(norm(O_C), 2) - pow(R, 2));
      intersection.intersected = discriminant >= 0;
      intersection.t = 0;
      if (intersection.intersected) {
          double t_1 = dot(ray.u, O_C*(-1.)) - sqrt(discriminant);
          double t_2 = dot(ray.u, O_C*(-1.)) + sqrt(discriminant);
          if (t_2 < 0) intersection.intersected = false;
          else intersection.t = (t_1 >= 0) ? t_1 : t_2;
      }
      intersection.P = ray.O + ray.u*intersection.t;
      intersection.N = normalize(intersection.P - C);
      intersection.albedo = albedo;
      intersection.sphere_center = this->C;
      intersection.refractive_index = refractive_index;
      if (this->reflective) intersection.reflective = true;
      if (this->invert_normal) intersection.N = (-1.)*intersection.N;
      if (this->is_light_source){
        intersection.is_light = true;
      }
      return intersection;
    }
    Vector albedo;
  private:
    Vector C;
    double R;
    double refractive_index;
    bool invert_normal;
    bool is_light_source;
};

// UTILITIES -------------------------------------------------------------------

void boxMuller(double stdev, double &x, double &y) {
  double r1 = ((double) rand() / (RAND_MAX));
  double r2 = ((double) rand() / (RAND_MAX));
  x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
  y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

Vector random_cos(const Vector &N) {
  double r1 = ((double) rand() / (RAND_MAX));
  double r2 = ((double) rand() / (RAND_MAX));
  double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
  double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
  double z = sqrt(r2);

  int min_i = 0;
  double min = abs(N[0]);
  for (int i = 1; i < 3; i++) {
    if (abs(N[i]) < min) {
      min = abs(N[i]);
      min_i = i;
    }
  }

  Vector T1;
  if (min_i == 0) T1 = Vector(0, N[2], -N[1]);
  else if (min_i == 1) T1 = Vector(N[2], 0, -N[0]);
  else if (min_i == 2) T1 = Vector(N[1], -N[0], 0);
  T1 = normalize(T1);

  Vector T2 = cross(T1, N);

  return T1*x + T2*y + N*z;
}

// MESH ------------------------------------------------------------------------

class TriangleIndices {
public:
	TriangleIndices(
	  	int vtxi = -1, int vtxj = -1, int vtxk = -1,
	  	int ni = -1, int nj = -1, int nk = -1,
	  	int uvi = -1, int uvj = -1, int uvk = -1,
	  	int group = -1
	) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk),
		uvi(uvi), uvj(uvj), uvk(uvk),
		ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

struct Node {
  Node *child_left;
  Node *child_right;
  BoundingBox bounding_box;
  int starting_triangle;
  int ending_triangle;
};

class TriangleMesh : public Geometry {
  public:
    ~TriangleMesh() {}
    TriangleMesh(
      double scaling_factor,
      const Vector &translation,
      const Vector &albedo,
      bool reflective = false
    ) {
      this->reflective = reflective;
      this->scaling_factor = scaling_factor;
      this->translation = translation;
      this->albedo = albedo;
      this->bvh_root = new Node;
      // TODO: Rotataion
    };

    void readOBJ(const char* obj) {

      char matfile[255];
      char grp[255];

      FILE* f;
      f = fopen(obj, "r");
      int curGroup = -1;
      while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
          sscanf(line, "usemtl %[^\n]\n", grp);
          curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ') {
          Vector vec;

          Vector col;
          if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
            col[0] = std::min(1., std::max(0., col[0]));
            col[1] = std::min(1., std::max(0., col[1]));
            col[2] = std::min(1., std::max(0., col[2]));

            vertices.push_back(vec);
            vertexcolors.push_back(col);

          } else {
            sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            vertices.push_back(vec);
          }
        }
        if (line[0] == 'v' && line[1] == 'n') {
          Vector vec;
          sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
          normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
          Vector vec;
          sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
          uvs.push_back(vec);
        }
        if (line[0] == 'f') {
          TriangleIndices t;
          int i0, i1, i2, i3;
          int j0, j1, j2, j3;
          int k0, k1, k2, k3;
          int nn;
          t.group = curGroup;

          char* consumedline = line + 1;
          int offset;

          nn = sscanf(consumedline, "%i/%i/%i %i/%i/%i %i/%i/%i%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
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
          } else {
            nn = sscanf(consumedline, "%i/%i %i/%i %i/%i%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
            if (nn == 6) {
              if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
              if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
              if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
              if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
              if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
              if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
              indices.push_back(t);
            } else {
              nn = sscanf(consumedline, "%i %i %i%n", &i0, &i1, &i2, &offset);
              if (nn == 3) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                indices.push_back(t);
              } else {
                nn = sscanf(consumedline, "%i//%i %i//%i %i//%i%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
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
            nn = sscanf(consumedline, "%i/%i/%i%n", &i3, &j3, &k3, &offset);
            TriangleIndices t2;
            t2.group = curGroup;
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
            } else {
              nn = sscanf(consumedline, "%i/%i%n", &i3, &j3, &offset);
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
              } else {
                nn = sscanf(consumedline, "%i//%i%n", &i3, &k3, &offset);
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
                } else {
                  nn = sscanf(consumedline, "%i%n", &i3, &offset);
                  if (nn == 1) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                    consumedline = consumedline + offset;
                    i2 = i3;
                    indices.push_back(t2);
                  } else {
                    consumedline = consumedline + 1;
                  }
                }
              }
            }
          }

        }

      }
      fclose(f);
      build_bvh(this->bvh_root, 0, indices.size());
    }

    Intersection intersect(const Ray &ray) override {
      Intersection intersection;
      intersection.intersected = false;
      double t;
      double min_t = MAXFLOAT;

      if (not intersect_bounding_box(ray, this->bvh_root->bounding_box, t))
        return intersection;

      // NO BVH
      // Vector A, B, C, N, e1, e2;
      // for (int i = 0; i < indices.size(); i++) {
      //   TriangleIndices triangle = this->indices[i];
      //   A = scaling_factor*vertices[triangle.vtxi] + translation;
      //   B = scaling_factor*vertices[triangle.vtxj] + translation;
      //   C = scaling_factor*vertices[triangle.vtxk] + translation;
      //   e1 = B - A;
      //   e2 = C - A;
      //   N = cross(e1, e2);

      //   double beta  =  dot(e2, cross(A - ray.O, ray.u)) / dot(ray.u, N);
      //   double gamma = -dot(e1, cross(A - ray.O, ray.u)) / dot(ray.u, N);
      //   double alpha = 1. - beta - gamma;

      //   if (alpha > 0. && beta > 0. && gamma > 0.) {
      //     double t = dot(A - ray.O, N) / dot(ray.u, N);
      //     if (0 < t && t < min_t) {
      //       min_t = t;
      //       intersection.intersected = true;
      //       intersection.t = t;
      //       intersection.P = A + beta*e1 + gamma*e2;
      //       intersection.N = N;
      //       if (this->reflective) intersection.reflective = true;
      //       intersection.albedo = this->albedo;
      //     }
      //   }
      // }
      // return intersection;

      // BVH
      std::list<Node*> nodes_to_visit;
      nodes_to_visit.push_front(this->bvh_root);
      while (not nodes_to_visit.empty()) {
        Node* curNode = nodes_to_visit.back();
        nodes_to_visit.pop_back();
        if (curNode->child_left) {
          if (intersect_bounding_box(ray, curNode->child_left->bounding_box, t)) {
            if (t < min_t) nodes_to_visit.push_back(curNode->child_left);
          }
          if (intersect_bounding_box(ray, curNode->child_right->bounding_box, t)) {
            if (t < min_t) nodes_to_visit.push_back(curNode->child_right);
          }
        } else {
          Vector A, B, C, N, e1, e2;
          for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++) {
            TriangleIndices triangle = this->indices[i];
            A = scaling_factor*vertices[triangle.vtxi] + translation;
            B = scaling_factor*vertices[triangle.vtxj] + translation;
            C = scaling_factor*vertices[triangle.vtxk] + translation;
            e1 = B - A;
            e2 = C - A;
            N = cross(e1, e2);

            double beta  =  dot(e2, cross(A - ray.O, ray.u)) / dot(ray.u, N);
            double gamma = -dot(e1, cross(A - ray.O, ray.u)) / dot(ray.u, N);
            double alpha = 1. - beta - gamma;

            if (alpha > 0. && beta > 0. && gamma > 0.) {
              double t = dot(A - ray.O, N) / dot(ray.u, N);
              if (0 < t && t < min_t) {
                min_t = t;
                intersection.intersected = true;
                intersection.t = t;
                intersection.P = A + beta*e1 + gamma*e2;
                intersection.N = N;
                if (this->reflective) intersection.reflective = true;
                intersection.albedo = this->albedo;
              }
            }
          }
        }
      }
      return intersection;
    }

  private:

    void build_bvh(Node *node, int starting_triangle, int ending_triangle) {
      node->bounding_box =  compute_bounding_box(starting_triangle, ending_triangle);
      node->starting_triangle = starting_triangle;
      node->ending_triangle = ending_triangle;

      Vector diag = node->bounding_box.B_max - node->bounding_box.B_min;
      Vector middle_diag = node->bounding_box.B_min + diag*0.5;

      int longest_axis = 0;
      double max = abs(diag[0]);
      for (int i = 1; i < 3; i++) {
        if (abs(diag[i]) > max) {
          max = abs(diag[i]);
          longest_axis = i;
        }
      }

      int pivot_index = starting_triangle;
      for (int i = starting_triangle; i < ending_triangle; i++) {
        Vector vertex1 = scaling_factor*this->vertices[this->indices[i].vtxi] + translation;
        Vector vertex2 = scaling_factor*this->vertices[this->indices[i].vtxj] + translation;
        Vector vertex3 = scaling_factor*this->vertices[this->indices[i].vtxk] + translation;
        Vector barycenter = (vertex1 + vertex2 + vertex3)/3.;
        if (barycenter[longest_axis] < middle_diag[longest_axis]) {
          std::swap(indices[i], indices[pivot_index]);
          pivot_index++;
        }
      }

      if ( // Stopping Criteria
        pivot_index <= starting_triangle ||
        pivot_index >= ending_triangle - 5 ||
        ending_triangle - starting_triangle < 5
      ) return;

      node->child_left = new Node();
      node->child_right = new Node();
      build_bvh(node->child_left, starting_triangle, pivot_index);
      build_bvh(node->child_right, pivot_index, ending_triangle);
    }

    BoundingBox compute_bounding_box(int starting_triangle, int ending_triangle) {
      double min_x = MAXFLOAT, min_y = MAXFLOAT, min_z = MAXFLOAT;
      double max_x = - MAXFLOAT, max_y = - MAXFLOAT, max_z = - MAXFLOAT;
      for (int i = starting_triangle; i < ending_triangle; i++) {
        auto original_triangle_vertices = {
          this->vertices[this->indices[i].vtxi],
          this->vertices[this->indices[i].vtxj],
          this->vertices[this->indices[i].vtxk]
        };
        for (auto const& vertex: original_triangle_vertices) {
          Vector V = scaling_factor*vertex + translation;
          if (V[0] < min_x) min_x = V[0];
          else if (V[0] > max_x) max_x = V[0];
          if (V[1] < min_y) min_y = V[1];
          else if (V[1] > max_y) max_y = V[1];
          if (V[2] < min_z) min_z = V[2];
          else if (V[2] > max_z) max_z = V[2];
        }
      }
      BoundingBox bounding_box;
      bounding_box.B_min = Vector(min_x, min_y, min_z);
      bounding_box.B_max = Vector(max_x, max_y, max_z);
      return bounding_box;
    }

    bool intersect_bounding_box(const Ray &ray, BoundingBox bounding_box, double &t) {
      double tx1, ty1, tz1;
      double tx0, ty0, tz0;
      double t_B_min, t_B_max;
      Vector N;

      N = Vector(1,0,0);
      t_B_min = dot(bounding_box.B_min - ray.O, N) / dot(ray.u, N);
      t_B_max = dot(bounding_box.B_max - ray.O, N) / dot(ray.u, N);
      tx0 = std::min(t_B_min, t_B_max);
      tx1 = std::max(t_B_min, t_B_max);

      N = Vector(0,1,0);
      t_B_min = dot(bounding_box.B_min - ray.O, N) / dot(ray.u, N);
      t_B_max = dot(bounding_box.B_max - ray.O, N) / dot(ray.u, N);
      ty0 = std::min(t_B_min, t_B_max);
      ty1 = std::max(t_B_min, t_B_max);

      N = Vector(0,0,1);
      t_B_min = dot(bounding_box.B_min - ray.O, N) / dot(ray.u, N);
      t_B_max = dot(bounding_box.B_max - ray.O, N) / dot(ray.u, N);
      tz0 = std::min(t_B_min, t_B_max);
      tz1 = std::max(t_B_min, t_B_max);

      double first_intersection_t = std::max({tx0,ty0,tz0});
      if (std::min({tx1,ty1,tz1}) > first_intersection_t > 0) {
        t = first_intersection_t;
        return true;
      }
      return false;
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    Vector albedo;
    Vector translation;
    double scaling_factor;
    BoundingBox full_bounding_box;
    Node* bvh_root;
};

// SCENE  ----------------------------------------------------------------------

class BasicScene {
  public:
    explicit BasicScene(
      const Vector& light_position,
      double light_radius,
      double light_intensity
    ){
      Sphere* left_wall = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
      Sphere* front_wall = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
      Sphere* right_wall = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
      Sphere* back_wall = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
      Sphere* ceiling = new Sphere(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
      Sphere* ground = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
      geometries.push_back(left_wall);
      geometries.push_back(front_wall);
      geometries.push_back(right_wall);
      geometries.push_back(back_wall);
      geometries.push_back(ceiling);
      geometries.push_back(ground);
      this->light_position = light_position;
      this->light_radius = light_radius;
      this->light_intensity = light_intensity;
    }

    void addGeometry(Geometry* geometry) {
      geometries.push_back(geometry);
    }

    Intersection intersect(const Ray& ray) {
      Intersection best_intersection;
      best_intersection.intersected = false;
      double min_t = MAXFLOAT;
      for (auto &geometry : geometries) {
        Intersection intersection = geometry->intersect(ray);
        if (intersection.intersected && intersection.t < min_t) {
          min_t = intersection.t;
          best_intersection = intersection;
        }
      }
      return best_intersection;
    }

    Vector getColor(const Ray& ray, int ray_depth, bool last_bounce_diffuse = false) {
      // TODO: Refactor
      if (ray_depth < 0) return Vector(0., 0., 0.);
      Intersection intersection = intersect(ray);
      Vector Lo(0., 0., 0.);

      if (intersection.intersected) {
        double eps = 1e-10;
        Vector N = intersection.N;
        Vector P = intersection.P + N*eps;

        if (intersection.is_light) {
          double R = light_radius;
          if (last_bounce_diffuse) return Vector(0., 0. ,0.);
          else return Vector(1., 1. ,1.)*light_intensity/(4*M_PI*M_PI*R*R);
        }

        if (intersection.reflective) {
          // Reflection
          Ray reflected_ray = Ray(P, ray.u - (2*dot(N,ray.u)) * N);
          return getColor(reflected_ray, ray_depth - 1);
        }
        else if (intersection.refractive_index != 1.) {
          // Refraction
          double n1, n2;
          if (dot(N,ray.u) > 0) {
            N = (-1.)*N;
            n1 = intersection.refractive_index;
            n2 = 1.;
          } else {
            n1 = 1.;
            n2 = intersection.refractive_index;
          }
          double k0 = pow(n1 - n2, 2.)/pow(n1 + n2, 2.);
          P = intersection.P - N*eps;
          double dot_N_u = dot(N,ray.u);
          if (1. - pow((n1/n2), 2.) * (1 - pow(dot_N_u, 2.)) > 0) {
            // Standard Refraction
            Vector w_T = (n1/n2) * (ray.u - dot_N_u*N);
            Vector w_N = (-1.)*N * sqrt(1 - pow((n1/n2),2.)*(1 - pow(dot_N_u,2.)));
            Vector w = w_T + w_N;
            double x = ((double) rand() / (RAND_MAX));
            if (x < k0 + (1-k0)*pow(1 - abs(dot(N,w)),5.)) {
              Ray reflected_ray = Ray(P, ray.u - (2*dot(intersection.N,ray.u)) * intersection.N);
              return getColor(reflected_ray, ray_depth - 1);
            } else {
              Ray refracted_ray = Ray(P, w);
              return getColor(refracted_ray, ray_depth - 1);
            }
          } else {
            // Total Internal Relection
            Ray internal_reflected_ray = Ray(P, ray.u - (2*dot(intersection.N,ray.u)) * intersection.N);
            return getColor(internal_reflected_ray, ray_depth - 1);
          }
        }
        else {
          double I = light_intensity;
          double R = light_radius;
          Vector x = intersection.sphere_center;

          // Direct Lighting
          Vector xprime = R*random_cos(normalize(x - light_position)) + light_position;
          Vector Nprime = normalize(xprime - light_position);
          double d = norm(xprime - P);
          Vector omega = normalize(xprime - P);
          Ray lightRay = Ray(P, omega);
          Intersection lightIntersection = intersect(lightRay);
          double visibility = (lightIntersection.intersected && lightIntersection.t <= d) ? 0. : 1.;
          double pdf = dot(Nprime, normalize(x-light_position))/(M_PI*R*R);
          Vector rho = intersection.albedo;
          Lo = I/(4*M_PI*M_PI*R*R)*rho/M_PI*visibility*std::max(dot(N,omega),0.)*std::max(dot(Nprime,(-1.)*omega),0.)/(pow(norm(xprime-P),2.)*pdf);

          // Indirect Lighting
          Ray random_ray = Ray(P, random_cos(N));
          Lo += rho * getColor(random_ray, ray_depth - 1, true);
        }
      }

      return Lo;
    }

  private:
    std::vector<Geometry*> geometries;
    Vector light_position;
    double light_radius;
    double light_intensity;
};

// MAIN ------------------------------------------------------------------------

int main() {

  auto start = std::chrono::high_resolution_clock::now();

  Vector light_position(0, 20, 5);
  double light_radius = 3; // > 0
  double light_intensity = 1e5;

  BasicScene scene = BasicScene(light_position, light_radius, light_intensity);

  // Spherical Light Source
  Sphere* light = new Sphere(light_position, light_radius, Vector(1., 0., 0.), false, 1., false, true);
  scene.addGeometry(light);

  // Cat
  TriangleMesh* cat = new TriangleMesh(0.6, Vector(-10, -10, 0), Vector(0.718, 0.514, 0.263), false);
  cat->readOBJ("objs/cat.obj");
  scene.addGeometry(cat);

  // Gun
  TriangleMesh* gun = new TriangleMesh(15, Vector(20, 5, 0), Vector(0.1, 0.1, 0.1), false);
  gun->readOBJ("objs/gun.obj");
  scene.addGeometry(gun);

  // White Sphere
  Sphere* white_sphere = new Sphere(Vector(20, 0, -20), 10, Vector(1., 1., 1.));
  scene.addGeometry(white_sphere);

  // Reflective Sphere
  Sphere* reflective_sphere = new Sphere(Vector(-5, 15, -15), 12, Vector(1., 1., 1.), true);
  scene.addGeometry(reflective_sphere);

  // Solid Refractive Sphere
  Sphere* solid_refractive_sphere = new Sphere(Vector(-17, -7, 15), 3, Vector(1., 1., 1.), false, 1.5);
  scene.addGeometry(solid_refractive_sphere);

  // Hollow Refractive Sphere
  Sphere* hollow_refractive_sphere_outer = new Sphere(Vector(4, -6, 30), 4, Vector(1., 1., 1.), false, 1.5);
  Sphere* hollow_refractive_sphere_inner = new Sphere(Vector(4, -6, 30), 3.75, Vector(1., 1., 1.), false, 1.5, true);
  scene.addGeometry(hollow_refractive_sphere_outer);
  scene.addGeometry(hollow_refractive_sphere_inner);

  int W = 512;
  int H = 512;
  int rays_per_pixel = 10;
  std::vector<unsigned char> image(W*H * 3, 0);
  double fov = 1.0472; // 60 deg
  Vector camera_position = Vector(0, 0, 55);
  double max_ray_depth = 5;
  double gamma = 2.2;

  #pragma omp parallel for schedule(dynamic,1)
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      Vector color;
      double x,y;

      for (int k = 0; k < rays_per_pixel; k++) {
        boxMuller(0.5,x,y);
        Vector V;
        V[0] = (camera_position[0] + (j + x) + 0.5 - W / 2);
        V[1] = (camera_position[1] - (i + y) - 0.5 + H / 2);
        V[2] = camera_position[2] - W / (2 * tan(fov / 2));

        Ray ray = Ray(camera_position, normalize(V - camera_position));
        color += scene.getColor(ray, max_ray_depth);
      }

      color /= rays_per_pixel;

      image[(i*W + j) * 3 + 0] = std::max(std::min(255., pow(color[0], 1./gamma)*255),0.);
      image[(i*W + j) * 3 + 1] = std::max(std::min(255., pow(color[1], 1./gamma)*255),0.);
      image[(i*W + j) * 3 + 2] = std::max(std::min(255., pow(color[2], 1./gamma)*255),0.);
    }
  }

  stbi_write_png("render.png", W, H, 3, &image[0], 0);

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
  std::cout << "Duration: " << duration.count() / 1000. << "s" << std::endl;

  return 0;
}