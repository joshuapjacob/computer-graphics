#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <cstring>
#include <stdio.h>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "vector.h"

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

class Ray {
  public:
    Ray(const Vector& O, const Vector& u) {
        this->O = O;
        this->u = u;
    }
    Vector O;
    Vector u;
};

class Geometry {
  public:
    virtual Intersection intersect(const Ray &ray) = 0;
    bool reflective;
};

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

class TriangleMesh : public Geometry {
public:
	~TriangleMesh() {}
	TriangleMesh(bool reflective = false) {
    this->reflective = reflective;
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

	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	Intersection intersect(const Ray &ray) override {
		Intersection intersection;
		intersection.intersected = false;

		Vector A, B, C, N;
		double min_t = MAXFLOAT;
		TriangleIndices closest_triangle;
		for(auto const& triangle: indices) {
			A = vertices[triangle.vtxi];
			B = vertices[triangle.vtxj];
			C = vertices[triangle.vtxk];
			N = cross(B - A, C - A);
			double t = dot(A - ray.O, N) / dot(ray.u, N);
			if (0 < t && t < min_t){
				min_t = t;
				closest_triangle = triangle;
			}
		}

		if (min_t == MAXFLOAT) return intersection;

		A = vertices[closest_triangle.vtxi];
		B = vertices[closest_triangle.vtxj];
		C = vertices[closest_triangle.vtxk];
		Vector e1 = B - A;
		Vector e2 = C - A;
		N = cross(e1, e2);

		double beta  =  dot(e2, cross(A - ray.O, ray.u)) / dot(ray.u, N);
		double gamma = -dot(e1, cross(A - ray.O, ray.u)) / dot(ray.u, N);
		double alpha = 1 - beta - gamma;

		if ((alpha > 0) && (beta > 0) && (gamma > 0)) {
			intersection.intersected = true;
			intersection.t = min_t;
			intersection.P = A + beta*e1 + gamma*e2;
			intersection.N = N;
      if (this->reflective) intersection.reflective = true;
			if (this->vertexcolors.size() > 0){
				intersection.albedo = vertexcolors[closest_triangle.group];
			}
      intersection.albedo = Vector(1.,1.,1.);
		}

		return intersection;
	}
	
};

class BasicScene {
  public:
    explicit BasicScene(
      const Vector& light_position,
      double light_radius,
      double light_intensity
    ){
      // Sphere* left_wall = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
      // Sphere* front_wall = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
      // Sphere* right_wall = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
      // Sphere* back_wall = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
      // Sphere* ceiling = new Sphere(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
      // Sphere* ground = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
      // geometries.push_back(left_wall);
      // geometries.push_back(front_wall);
      // geometries.push_back(right_wall);
      // geometries.push_back(back_wall);
      // geometries.push_back(ceiling);
      // geometries.push_back(ground);
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
          Lo = intersection.albedo;

          
          // double R = light_radius;
          // Vector x = intersection.sphere_center;

          // // Direct Lighting
          // Vector xprime = R*random_cos(normalize(x - light_position)) + light_position;
          // Vector Nprime = normalize(xprime - light_position);
          // double d = norm(xprime - P);
          // Vector omega = normalize(xprime - P);
          // Ray lightRay = Ray(P, omega);
          // Intersection lightIntersection = intersect(lightRay);
          // double visibility = (lightIntersection.intersected && lightIntersection.t <= d) ? 0. : 1.;
          // double pdf = dot(Nprime, normalize(x-light_position))/(M_PI*R*R);
          // Vector rho = intersection.albedo;
          // Lo = I/(4*M_PI*M_PI*R*R)*rho/M_PI*visibility*std::max(dot(N,omega),0.)*std::max(dot(Nprime,(-1.)*omega),0.)/(pow(norm(xprime-P),2.)*pdf);
          
          // // Indirect Lighting
          // Ray random_ray = Ray(P, random_cos(N));
          // Lo += rho * getColor(random_ray, ray_depth - 1, true);
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

int main() {

  auto start = std::chrono::high_resolution_clock::now();
  
  Vector light_position(0, 15, 0);
  double light_radius = 2;
  double light_intensity = 1e5;

  BasicScene scene = BasicScene(light_position, light_radius, light_intensity);

  // Spherical Light Source
  Sphere* light = new Sphere(light_position, light_radius, Vector(1., 1., 1.), false, 1., false, true);
  scene.addGeometry(light);

  // Pyramid
  TriangleMesh* pyramid = new TriangleMesh(false);
  pyramid->readOBJ("objs/triangle.obj");
  scene.addGeometry(pyramid);

  // Cat
  TriangleMesh* cat = new TriangleMesh(true);
  cat->readOBJ("objs/simple_cat.obj");
  scene.addGeometry(cat);

  // White Sphere
  // Sphere* white_sphere = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.));
  // scene.addGeometry(white_sphere);

  // Reflective Sphere
  // Sphere* reflective_sphere = new Sphere(Vector(-20, 0, 0), 10, Vector(1., 1., 1.), true);
  // scene.addSphere(reflective_sphere);

  // Solid Refractive Sphere
  // Sphere* solid_refractive_sphere = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
  // scene.addSphere(solid_refractive_sphere);

  // Hollow Refractive Sphere
  // Sphere* hollow_refractive_sphere_outer = new Sphere(Vector(20, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
  // Sphere* hollow_refractive_sphere_inner = new Sphere(Vector(20, 0, 0), 9.5, Vector(1., 1., 1.), false, 1.5, true);
  // scene.addSphere(hollow_refractive_sphere_outer);
  // scene.addSphere(hollow_refractive_sphere_inner);

  int W = 512;
  int H = 512;
  int rays_per_pixel = 1;
  std::vector<unsigned char> image(W*H * 3, 0);
  double fov = 1.0472; // 60 deg
  Vector Q = Vector(0, 0, 100);
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