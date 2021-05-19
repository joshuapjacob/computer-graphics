#include <vector>
#include <stdexcept>

#include "vector.h"
#include "svg.h"

Vector intersect(const Vector &A, const Vector &B, const std::vector<Vector> &edge) {
    Vector u = edge[0];
    Vector v = edge[1];
    Vector N = normalize(Vector(v[1] - u[1], u[0] - v[0]));
    double t = dot(u - A, N) / dot(B - A, N);
    if (0 <= t && t <= 1) return A + t*(B - A);
    return Vector();
}

bool inside(const Vector &P, const std::vector<Vector> &edge) {
    Vector u = edge[0];
    Vector v = edge[1];
    Vector N = Vector(v[1] - u[1], u[0] - v[0]);
    if (dot(P - u, N) <= 0) return true;
    return false;
}

// Sutherland-Hodgman Polygon Clipping
Polygon clip(Polygon subjectPolygon, Polygon clipPolygon) {
    Polygon outPolygon;
    for (size_t i = 0; i < clipPolygon.vertices.size(); i++) {
        Vector u = clipPolygon.vertices[i];
        Vector v = clipPolygon.vertices[
            (i > 0) ? (i - 1) : (clipPolygon.vertices.size() - 1)
        ];
        std::vector<Vector> edge = {u, v};
        outPolygon = Polygon();
        for (size_t j = 0; j < subjectPolygon.vertices.size(); j++) {
            Vector curVertex = subjectPolygon.vertices[j];
            Vector prevVertex = subjectPolygon.vertices[
                (j > 0) ? (j - 1) : (subjectPolygon.vertices.size() - 1)
            ];
            Vector intersection = intersect(prevVertex, curVertex, edge);
            if (inside(curVertex, edge)) {
                if (not inside(prevVertex, edge)) {
                    outPolygon.add(intersection);
                }
                outPolygon.add(curVertex);
            } else if (inside(prevVertex, edge)) {
                outPolygon.add(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}

// MAIN ------------------------------------------------------------------------

int main() {

    // TEST: Sutherland-Hodgman Polygon Clipping
    Polygon subjectPolygon({
        Vector(0.1, 0.2), Vector(0.3, 0.8), Vector(0.6, 0.5), 
        Vector(0.8, 0.9), Vector(0.9, 0.1), Vector(0.5, 0.4)
    });
    Polygon clipPolygon({
        Vector(0.3, 0.3), Vector(0.3, 0.7), 
        Vector(0.7, 0.7), Vector(0.7, 0.3)
    });
    save_svg({subjectPolygon, clipPolygon}, "imgs/before_clipping.svg");
    save_svg({clip(subjectPolygon, clipPolygon)}, "imgs/after_clipping.svg");
}
