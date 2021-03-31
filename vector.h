#pragma once

class Vector {
  public:
    explicit Vector(double x = 0., double y = 0., double z = 0.);
    Vector& operator+=(const Vector& b);
    Vector &operator-=(const Vector& b);
    Vector& operator*=(double t);
    Vector& operator/=(double t);
    const double& operator[](int i) const;
    double& operator[](int i);

  private:
    double coords[3];
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator/(const Vector& a, double t);
Vector operator*(const Vector& a, double t);
Vector operator*(double t, const Vector& a);
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector &a, const Vector &b);
double norm(const Vector& a);
Vector normalize(const Vector& a);