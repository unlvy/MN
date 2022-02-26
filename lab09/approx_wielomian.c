#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdio.h>

#define x0 2
#define sigma 4
#define a0 -0.125
#define a1 0.125
#define a2 -0.03125

double delta(double alpha);
double g(double x);
double g2(double x, double alpha);
double uppercaseG(double x, gsl_vector *b);
double f(double x);
void approximate(int N, int m, double alpha, FILE* fPtr);

/* ----------------------------------------- */

int main() {

    FILE *fPtr1 = fopen("N_11_alpha_0.0.txt", "w");
    FILE *fPtr2 = fopen("N_11_alpha_0.5.txt", "w");
    FILE *fPtr3 = fopen("N_101_alpha_0.0.txt", "w");
    FILE *fPtr4 = fopen("N_101_alpha_0.5.txt", "w");

    approximate(11, 4, 0.0, fPtr1);
    approximate(11, 4, 0.5, fPtr2);
    approximate(101, 4, 0.0, fPtr3);
    approximate(101, 4, 0.5, fPtr4);

    fclose(fPtr1);
    fclose(fPtr2);
    fclose(fPtr3);
    fclose(fPtr4);
}

/* ----------------------------------------- */

double delta(double alpha) {
    return alpha * (rand() / ((double)RAND_MAX + 1.0) - 0.5);
}

double g(double x) {
    return exp(a0 + a1 * x + a2 * pow(x, 2));
}

double g2(double x, double alpha) {
    return g(x) * (1.0 + delta(alpha));
}

double uppercaseG(double x, gsl_vector *b) {
    return exp(
                gsl_vector_get(b, 0) + 
                gsl_vector_get(b, 1) * x + 
                gsl_vector_get(b, 2) * pow(x, 2) +
                gsl_vector_get(b, 3) * pow(x, 3)
               );
}

double f(double x) {
    return a0 + a1 * x + a2 * pow(x, 2);
}

void approximate(int N, int m, double alpha, FILE* fPtr) {

    gsl_vector *xVec = gsl_vector_calloc(N);
    gsl_vector *gVec = gsl_vector_calloc(N);
    gsl_vector *fVec = gsl_vector_calloc(N);

    double xmin = -3.0 * sigma + x0;
    double xmax = 3.0 * sigma + x0;
    double step = (xmax - xmin) / (N - 1);
    double value;

    if (fPtr) {
        fprintf(fPtr, "gj = [");
        for (int i = 0; i < N; i++) {
            value = xmin + i * step;
            gsl_vector_set(xVec, i, value);
            gsl_vector_set(gVec, i, g2(value, alpha));
            gsl_vector_set(fVec, i, log(gsl_vector_get(gVec, i)));
            fprintf(fPtr, "%.6lf ", gsl_vector_get(gVec, i));
        }
        fprintf(fPtr, "]\n");
    }

    gsl_matrix *G = gsl_matrix_calloc(m, m);
    gsl_vector *rVec = gsl_vector_calloc(m);

    double sum;
    for (int k = 0; k < m; k++) {
        sum = 0.0;
        for (int j = 0; j < N; j++) {
            sum += gsl_vector_get(fVec, j) * pow(gsl_vector_get(xVec, j), k);
        }
        printf("%.6lf ", sum);
        gsl_vector_set(rVec, k, sum);
        for (int i = 0; i < m; i++) {
            sum = 0.0;
            for (int j = 0; j < N; j++) {
                sum += pow(gsl_vector_get(xVec, j), i + k);
            }
            gsl_matrix_set(G, i, k, sum);
        }
        printf("\n");
    }

    gsl_linalg_HH_svx(G, rVec);

    if (fPtr) {
        fprintf(fPtr, "bi = [");
        for (int i = 0; i < m; i++) {
            fprintf(fPtr, "%g ", gsl_vector_get(rVec, i));
        }
        fprintf(fPtr, "]\n");

        fprintf(fPtr, "G(x) = ");
        for (double x = xmin; x <= xmax; x += 0.01) {
            fprintf(fPtr, "%.6lf ", uppercaseG(x, rVec));
        }
        fprintf(fPtr, "]\n");
    }

    gsl_vector_free(xVec);
    gsl_vector_free(gVec);
    gsl_vector_free(fVec);
    gsl_matrix_free(G);
}
