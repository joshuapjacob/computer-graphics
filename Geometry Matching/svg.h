#pragma once

#include <string>
#include <vector>

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
	std::vector<Vector> vertices;
};

void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none");

void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes);
