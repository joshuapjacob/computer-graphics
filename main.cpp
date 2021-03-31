#include <iostream>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "vector.h"

struct Intersection {
  bool intersected;
  double t;
  Vector P;
  Vector N;
  Vector albedo;
};

class Ray {
  public:
    Ray(const Vector& O, const Vector& u) {
        this->O = O;
        this->u = u;
    }
    Vector O;
    Vector u;
};

class Sphere {
  public:
    Sphere(const Vector& C, const double& R, const Vector& albedo) {
        this->C = C;
        this->R = R;
        this->albedo = albedo;
    }
    Intersection intersect(const Ray &ray) {
      Intersection intersection;
      Vector O_C = ray.O - C;
      double discriminant = pow(dot(ray.u, O_C), 2) - (dot(O_C, O_C) - pow(R, 2));
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
      return intersection;
    }
    Vector albedo;
  private:
    Vector C;
    double R;
};

class BasicScene {
  public:
    BasicScene(){
      Sphere* left_wall = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
      Sphere* front_wall = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
      Sphere* right_wall = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
      Sphere* back_wall = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
      Sphere* ceiling = new Sphere(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
      Sphere* ground = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
      spheres.push_back(left_wall);
      spheres.push_back(front_wall);
      spheres.push_back(right_wall);
      spheres.push_back(back_wall);
      spheres.push_back(ceiling);
      spheres.push_back(ground);
    }
    void addSphere(Sphere* sphere) {
      spheres.push_back(sphere);
    }

    Vector getColor(const Ray& ray, int ray_depth) {
      if (ray_depth < 0) return Vector(0., 0., 0.);
      double min_t = 1e5;
      Sphere* closest_sphere;
      for (auto &sphere : spheres) {
        Intersection intersection = sphere->intersect(ray);
        if (intersection.intersected && intersection.t < min_t) {
          min_t = intersection.t;
          closest_sphere = sphere;
        }
      }
      return closest_sphere->albedo;
    }
  private:
    std::vector<Sphere*> spheres;
};

int main() {
  
  BasicScene scene = BasicScene();

  Sphere* sphere = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.));
  scene.addSphere(sphere);

  int W = 1024;
  int H = 1024;
  std::vector<unsigned char> image(W*H * 3, 0);
  double fov = 1.0472; // 60 deg
  Vector Q = Vector(0, 0, 55);
  double max_path_length = 1e5;
  double gamma = 2.2;

  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      Vector V;
      V[0] = Q[0] + j + 0.5 - W / 2;
      V[1] = Q[1] - i - 0.5 + H / 2;
      V[2] = Q[2] - W / (2 * tan(fov / 2));

      Ray ray = Ray(Q, normalize(V - Q));
      Vector color = scene.getColor(ray, max_path_length);
            
      image[(i*W + j) * 3 + 0] = std::min(255., pow(color[0], 1./gamma)*255);
      image[(i*W + j) * 3 + 1] = std::min(255., pow(color[1], 1./gamma)*255);
      image[(i*W + j) * 3 + 2] = std::min(255., pow(color[2], 1./gamma)*255);
    }
  }
  
  stbi_write_png("render.png", W, H, 3, &image[0], 0);

  return 0;
}