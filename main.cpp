#include <iostream>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "vector.h"

struct Intersection {
  bool intersected = false;
  bool reflective = false;
  double refractive_index = 1.;
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
    Sphere(const Vector& C,
           const double& R,
           const Vector& albedo,
           bool reflective = false,
           double refractive_index = 1.,
           bool hollow_inner_sphere = false) {
        this->C = C;
        this->R = R;
        this->albedo = albedo;
        this->reflective = reflective;
        this->refractive_index = refractive_index;
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
      intersection.albedo = albedo;
      intersection.refractive_index = refractive_index;
      if (this->reflective) intersection.reflective = true;
      return intersection;
    }
    Vector albedo;
  private:
    Vector C;
    double R;
    bool reflective;
    double refractive_index;
};

class BasicScene {
  public:
    explicit BasicScene(const Vector& light_position){
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
      this->S = light_position;
    }

    void addSphere(Sphere* sphere) {
      spheres.push_back(sphere);
    }

    Intersection intersect(const Ray& ray) {
      Intersection best_intersection;
      best_intersection.intersected = false;
      double min_t = 1e5;
      for (auto &sphere : spheres) {
        Intersection intersection = sphere->intersect(ray);
        if (intersection.intersected && intersection.t < min_t) {
          min_t = intersection.t;
          best_intersection = intersection;
        }
      }
      return best_intersection;
    }

    Vector getColor(const Ray& ray, int ray_depth) {
      if (ray_depth < 0) return Vector(0., 0., 0.);
      Intersection intersection = intersect(ray);
      Vector L;

      if (intersection.intersected) {
        double eps = 1e-10;
        Vector N = intersection.N;
        Vector P = intersection.P + N*eps;

        if (intersection.reflective) {
          Ray reflected_ray = Ray(P, ray.u - (2*dot(N,ray.u)) * N);
          return getColor(reflected_ray, ray_depth - 1);

        } else if (intersection.refractive_index != 1.) {
          double n1, n2;
          if (dot(N,ray.u) > 0) {
            N = (-1.)*N;
            n1 = intersection.refractive_index;
            n2 = 1.;
          } else {
            n1 = 1.;
            n2 = intersection.refractive_index;
          }
          P = intersection.P - N*eps;
          double dot_N_u = dot(N,ray.u);
          Vector u_T = (n1/n2) * (ray.u - dot_N_u*N);
          Vector u_N = (-1.)*N * sqrt(1 - pow((n1/n2),2.)*(1 - pow(dot_N_u,2.)));
          Vector u = u_T + u_N;
          Ray refracted_ray = Ray(P, u);
          return getColor(refracted_ray, ray_depth - 1);

        } else {
          double d = norm(S - P);
          Vector omega = normalize(S - P);
          Ray lightRay = Ray(S, omega*(-1.));
          Intersection lightIntersection = intersect(lightRay);
          bool V_p = !(lightIntersection.intersected && lightIntersection.t <= d);
          Vector rho = lightIntersection.albedo;
          L = rho / M_PI * I / (4 * M_PI * pow(d, 2)) * V_p * std::max(dot(N, omega), 0.);
        }
      }

      return L;
    }

  private:
    std::vector<Sphere*> spheres;
    Vector S;
    double I = 1e5;
};

int main() {
  
  BasicScene scene = BasicScene(Vector(-10, 20, 40));

  // Sphere* sphere1 = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5, false);
  // scene.addSphere(sphere1);
  Sphere* sphere2 = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
  scene.addSphere(sphere2);

  int W = 1024;
  int H = 1024;
  std::vector<unsigned char> image(W*H * 3, 0);
  double fov = 1.0472; // 60 deg
  Vector Q = Vector(0, 0, 55);
  double max_ray_depth = 10;
  double gamma = 2.2;

  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      Vector V;
      V[0] = Q[0] + j + 0.5 - W / 2;
      V[1] = Q[1] - i - 0.5 + H / 2;
      V[2] = Q[2] - W / (2 * tan(fov / 2));

      Ray ray = Ray(Q, normalize(V - Q));
      Vector color = scene.getColor(ray, max_ray_depth);
            
      image[(i*W + j) * 3 + 0] = std::min(255., pow(color[0], 1./gamma)*255);
      image[(i*W + j) * 3 + 1] = std::min(255., pow(color[1], 1./gamma)*255);
      image[(i*W + j) * 3 + 2] = std::min(255., pow(color[2], 1./gamma)*255);
    }
  }
  
  stbi_write_png("render.png", W, H, 3, &image[0], 0);

  return 0;
}