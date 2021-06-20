#pragma once
#include <lbfgs.h>
#include <stdio.h>
#include <functional>
#include <vector>
#include "vector.h"
#include "svg.h"

class objective_function {
protected:
    lbfgsfloatval_t *m_x;
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    );
    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    );
    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    );
    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    );
    std::vector<double> lambdas;
    int iterations;
    std::function<std::vector<Polygon>(
        const std::vector<Vector>&,
        const Polygon&,
        const std::vector<double>&
    )> voronoi_function;
    std::vector<Vector> points;
    Polygon bounds;
    std::vector<Polygon> polygons;
public:
    objective_function(
        const std::vector<double> &lambdas,
        std::function<std::vector<Polygon>(
            const std::vector<Vector>&,
            const Polygon&,
            const std::vector<double>&
        )> voronoi_function,
        const std::vector<Vector> &points,
        const Polygon &bounds
    );
    virtual ~objective_function();
    std::vector<Polygon> run(int N);
};