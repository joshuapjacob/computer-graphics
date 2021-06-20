#include "lbfgs.h"
#include <vector>
#include "vector.h"
#include "svg.h"

objective_function::objective_function(
        const std::vector<double> &lambdas,
        std::function<std::vector<Polygon>(
            const std::vector<Vector>&,
            const Polygon&,
            const std::vector<double>&
        )> voronoi_function,
        const std::vector<Vector> &points,
        const Polygon &bounds
    ) {
    this->m_x = NULL;
    this->lambdas = lambdas;
    this->voronoi_function = voronoi_function;
    this->points = points;
    this->bounds = bounds;
    this->iterations = 0;
}

objective_function::~objective_function() {
    if (m_x != NULL) {
        lbfgs_free(m_x);
        m_x = NULL;
    }
}

std::vector<Polygon> objective_function::run(int N)
{
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(N);

    if (m_x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        exit(1);
    }

    /* Initialize the variables. */
    for (int i = 0; i < N; i ++) {
        m_x[i] = -1;
    }

    /*
        Start the L-BFGS optimization; this will invoke the callback functions
        evaluate() and progress() when necessary.
        */
    int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, NULL);

    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);

    return this->polygons;
}

lbfgsfloatval_t objective_function::_evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ) {
    return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
}

lbfgsfloatval_t objective_function::evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ) {
    lbfgsfloatval_t fx = 0.0;
    std::vector<double> weights(x, x + n);
    this->polygons = this->voronoi_function(points, bounds, weights);
    if (iterations % 100 == 0) {
        std::string filename = "imgs/optimized_" + std::to_string(iterations) + ".svg";
        // save_svg(polygons, filename);
    }
    for (uint i = 0; i < n; i++) {
        std::vector<Vector> vertices = this->polygons[i].vertices;
        size_t n = vertices.size();
        Vector point = this->points[i];
        double area = this->polygons[i].area();
        double tmp = 0;
        if (n > 0) {
            Vector c1 = vertices[0];
            for (uint i = 0; i < n - 2; i++) {
                Vector c2 = vertices[i + 1];
                Vector c3 = vertices[i + 2];
                double T = Polygon({c1,c2,c3}).area();
                tmp += (T/6.) * (
                    dot(c1 - point, c1 - point) +
                    dot(c1 - point, c2 - point) +
                    dot(c1 - point, c3 - point) +
                    dot(c2 - point, c2 - point) +
                    dot(c2 - point, c3 - point) +
                    dot(c3 - point, c3 - point)
                );
            }
        }
        fx += tmp - x[i]*area + this->lambdas[i]*x[i];
        g[i] = area - this->lambdas[i];
    }
    iterations += 1;
    return -fx;
}

int objective_function::_progress(
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
    ) {
    return reinterpret_cast<objective_function*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

int objective_function::progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    ) {
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}