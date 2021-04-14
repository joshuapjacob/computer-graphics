#include <iostream>
#include <vector>
#include <chrono>

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
           bool invert_normal = false) {
        this->C = C;
        this->R = R;
        this->albedo = albedo;
        this->reflective = reflective;
        this->refractive_index = refractive_index;
        this->invert_normal = invert_normal;
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
      if (this->invert_normal) intersection.N = (-1.)*intersection.N;
      return intersection;
    }
    Vector albedo;
  private:
    Vector C;
    double R;
    bool reflective;
    double refractive_index;
    bool invert_normal;
};

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

  return x * T1 + y * T2 + z * N;
}

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
      // TODO: Refactor this function.
      if (ray_depth < 0) return Vector(0., 0., 0.);
      Intersection intersection = intersect(ray);
      Vector Lo;

      if (intersection.intersected) {
        double eps = 1e-10;
        Vector N = intersection.N;
        Vector P = intersection.P + N*eps;

        if (intersection.reflective) { 
          // Reflection
          Ray reflected_ray = Ray(P, ray.u - (2*dot(N,ray.u)) * N);
          return getColor(reflected_ray, ray_depth - 1);

        } else if (intersection.refractive_index != 1.) { 
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

        } else { 
          // Diffuse
          double d = norm(S - P);
          Vector omega = normalize(S - P);
          Ray lightRay = Ray(S, omega*(-1.));
          Intersection lightIntersection = intersect(lightRay);

          // Direct Lighting
          double visibility = (lightIntersection.intersected && lightIntersection.t <= d) ? 0 : 1;
          Vector rho = lightIntersection.albedo;
          Lo = I/(4*M_PI*pow(d, 2)) * rho/M_PI * visibility * std::max(dot(N,omega), 0.);
          
          // TODO: Fix Indirect Lighting
          Ray random_ray = Ray(P, random_cos(N));
          Lo += rho * getColor(random_ray, ray_depth - 1);
        }
      }

      return Lo;
    }
  
  private:
    std::vector<Sphere*> spheres;
    Vector S;
    double I = 1e5;
};

int main() {

  auto start = std::chrono::high_resolution_clock::now();
  
  BasicScene scene = BasicScene(Vector(-10, 20, 40));

  // White Sphere
  Sphere* white_sphere = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.));
  scene.addSphere(white_sphere);

  // // Reflective Sphere
  // Sphere* reflective_sphere = new Sphere(Vector(-20, 0, 0), 10, Vector(1., 1., 1.), true);
  // scene.addSphere(reflective_sphere);

  // // Solid Refractive Sphere
  // Sphere* solid_refractive_sphere = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
  // scene.addSphere(solid_refractive_sphere);

  // // Hollow Refractive Sphere
  // Sphere* hollow_refractive_sphere_outer = new Sphere(Vector(20, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
  // Sphere* hollow_refractive_sphere_inner = new Sphere(Vector(20, 0, 0), 9.5, Vector(1., 1., 1.), false, 1.5, true);
  // scene.addSphere(hollow_refractive_sphere_outer);
  // scene.addSphere(hollow_refractive_sphere_inner);

  int W = 512;
  int H = 512;
  int rays_per_pixel = 10;
  std::vector<unsigned char> image(W*H * 3, 0);
  double fov = 1.0472; // 60 deg
  Vector Q = Vector(0, 0, 55);
  double max_ray_depth = 5;
  double gamma = 2.2;

  #pragma omp parallel for
  for (int i = 0; i < H; i++) {
    for (int j = 0; j < W; j++) {

      Vector color;
      double x,y;
      
      for (int k = 0; k < rays_per_pixel; k++) {
        boxMuller(0.5,x,y);
        Vector V;
        V[0] = (Q[0] + (j + x) + 0.5 - W / 2);
        V[1] = (Q[1] - (i + y) - 0.5 + H / 2);
        V[2] = Q[2] - W / (2 * tan(fov / 2));

        Ray ray = Ray(Q, normalize(V - Q));
        color += scene.getColor(ray, max_ray_depth);
      }

      color /= rays_per_pixel;
            
      image[(i*W + j) * 3 + 0] = std::min(255., pow(color[0], 1./gamma)*255);
      image[(i*W + j) * 3 + 1] = std::min(255., pow(color[1], 1./gamma)*255);
      image[(i*W + j) * 3 + 2] = std::min(255., pow(color[2], 1./gamma)*255);
    }
  }
  
  stbi_write_png("render.png", W, H, 3, &image[0], 0);

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
  std::cout << "Duration: " << duration.count() / 1000. << "s" << std::endl;

  return 0;
}