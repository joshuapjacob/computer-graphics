#include "lbfgs.h"

objective_function::~objective_function() {
    if (m_x != NULL) {
        lbfgs_free(m_x);
        m_x = NULL;
    }
}

int objective_function::run(int N)
{
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(N);

    if (m_x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }

    /* Initialize the variables. */
    for (int i = 0;i < N;i += 2) {
        m_x[i] = -1.2;
        m_x[i+1] = 1.0;
    }

    /*
        Start the L-BFGS optimization; this will invoke the callback functions
        evaluate() and progress() when necessary.
        */
    int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, NULL);

    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);
    
    return ret;
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

    return fx;
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
    // printf("Iteration %d:\n", k);
    // printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    // printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    // printf("\n");
    return 0;
}