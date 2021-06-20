#pragma once

#include <string>
#include <vector>
#pragma once
#include "vector.h"

class Polygon {
public:
    Polygon() {};
    explicit Polygon(const std::vector<Vector> &vertices) {
        this->vertices = vertices;
    };
    void add(const Vector& a) {
        this->vertices.push_back(a);
    };
    double area();
    Vector centroid();
	std::vector<Vector> vertices;
};

void save_svg(
    const std::vector<Polygon> &polygons,
    std::string filename,
    std::string fillcol = "none"
);

// void save_svg_animated(
//     const std::vector<Polygon> &polygons,
//     std::string filename,
//     int frameid,
//     int nbframes
// );

void save_frame(
    const std::vector<Polygon> &cells,
    std::string filename, int frameid,
    uint N
);
