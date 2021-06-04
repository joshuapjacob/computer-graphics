#include <vector>
#include <cmath>
#include <iostream>

#include <lbfgs.h>

#include "vector.h"
#include "svg.h"

// POLYGON CLIPPING ------------------------------------------------------------

Vector intersect(
    const Vector &A, const Vector &B,
    const std::vector<Vector> &edge
) {
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
    for (uint i = 0; i < clipPolygon.vertices.size(); i++) {
        Vector u = clipPolygon.vertices[i];
        Vector v = clipPolygon.vertices[
            (i > 0) ? (i - 1) : (clipPolygon.vertices.size() - 1)
        ];
        std::vector<Vector> edge = {u, v};
        outPolygon = Polygon();
        for (uint j = 0; j < subjectPolygon.vertices.size(); j++) {
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

// VORONOI ---------------------------------------------------------------------

Vector voronoi_intersect(
    const Vector &A, const Vector &B,
    const std::vector<Vector> &edge,
    const std::vector<double> &weight
) {
    Vector u = edge[0];
    Vector v = edge[1];
    Vector N = (u + v)/2 
        + (weight[0] - weight[1])/(2*std::pow(norm(u - v), 2))*(v - u);
    double t = dot(N - A, u - v)/dot(B - A, u - v);
    if (0 <= t && t <= 1) return A + t*(B - A);
    return Vector();
}

bool voronoi_inside(
    const Vector &P,
    const std::vector<Vector> &edge,
    const std::vector<double> &weight
) {
    Vector u = edge[0];
    Vector v = edge[1];
    Vector N = (u + v)/2 
        + (weight[0] - weight[1])/(2*std::pow(norm(u - v), 2))*(v - u);
    if (dot(P - N, v - u) < 0) return true;
    return false;
}

 // Voronoi Parallel Linear Enumeration
 std::vector<Polygon> voronoi_diagram(
    const std::vector<Vector> &points,
    const Polygon &bounds,
    const std::vector<double> &weights = {}
) {
    std::vector<Polygon> result(points.size());
    #pragma omp parallel for
    for(uint i = 0; i < points.size(); i++) {
        Vector Pi = points[i];
        Polygon cell = bounds;
        for(uint j = 0; j < points.size(); j++) {
            if(i == j) continue;
            Vector Pj = points[j];
            std::vector<Vector> edge = {Pi, Pj};
            std::vector<double> weight = {1, 1};
            if (weights.size()) weight = {weights[i], weights[j]};
            Polygon outPolygon;
            for(uint k = 0; k < cell.vertices.size(); k++) {
                Vector curVertex = cell.vertices[k];
                Vector prevVertex = cell.vertices[
                    (k > 0)? (k - 1) : (cell.vertices.size() - 1)
                ];
                Vector intersection = voronoi_intersect(
                    prevVertex, curVertex, edge, weight
                );
                if (voronoi_inside(curVertex, edge, weight)) {
                    if (not voronoi_inside(prevVertex, edge, weight)) {
                        outPolygon.add(intersection);
                    }
                    outPolygon.add(curVertex);
                } else if (voronoi_inside(prevVertex, edge, weight)) {
                    outPolygon.add(intersection);
                }
            }
            cell = outPolygon;
        }
        result[i] = cell;
    }
    return result;
}

// L-BFGS WEIGHT OPTIMIZATION --------------------------------------------------

void weight_optimization(
    std::vector<double> &weights
) {
    // lbfgsfloatval_t fx;
    uint N = weights.size();
    lbfgsfloatval_t *m_weights = lbfgs_malloc(N);
    for (uint i = 0; i < N; i++) {
        m_weights[i] = weights[i];
    }

    lbfgs(N, m_weights, &fx, evaluate, _progress, instance, NULL);

    for (uint i = 0; i < N; i++) {
        weights[i] = m_weights[i];
        std::cout << m_weights[i] << std::endl;
    }
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

    // TEST: Voronoi Diagram
    uint n = 10;
    Polygon bounds({
        Vector(0, 0), Vector(0, 1), 
        Vector(1, 1), Vector(1, 0)
    });
    std::vector<Vector> points(n);
    srand(time(NULL));
    for (uint i = 0; i < n; i++) {
        double x = (double) rand() / RAND_MAX;
        double y = (double) rand() / RAND_MAX;
        points[i] = Vector(x, y);
    }
    save_svg(voronoi_diagram(points, bounds), "imgs/voronoi_diagram.svg");

    // TEST: Power Diagram
    std::vector<double> weights(n);
    for (uint i = 0; i < n; i++) {
        if (points[i][0] < 0.1 || points[i][0] > 0.9 || points[i][1] < 0.1 || points[i][1] > 0.9) {
            weights[i] = 0.9;
        } else {
            weights[i] = 1;
        }
    }
    save_svg(voronoi_diagram(points, bounds, weights), "imgs/power_diagram.svg");

    // TEST: Power Diagram Weight Optimization
    weight_optimization(weights);
    save_svg(voronoi_diagram(points, bounds, weights), "imgs/power_diagram_optimized.svg");


}
