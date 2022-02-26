#include <gsl/gsl_linalg.h>
#include <math.h>

#define IT_MAX 300
#define N 7

void matrixXvector(gsl_vector *dest, gsl_matrix *matrix, gsl_vector *vector, int n);
double vectorXvector(gsl_vector *vector1, gsl_vector *vector2, int n);
void normalizeVector(gsl_vector *dest, gsl_vector *src, int n);
void insertColumn(gsl_vector *vector, gsl_matrix *matrix, int k, int n);
void matrixTXmatrixXmatrix(gsl_matrix *dest, gsl_matrix *matrix1, gsl_matrix *matrix2, int n);
void orthogonalization(gsl_vector *vector, gsl_matrix *matrix, int k, int n);

/* ------------------------- ---------------------------------------*/

int main(){

    /* wszystkie potrzebne zmienne */
    gsl_matrix* A = gsl_matrix_calloc(N, N);
    gsl_matrix* X = gsl_matrix_calloc(N, N);
    gsl_matrix* D = gsl_matrix_calloc(N, N);
    gsl_vector* x = gsl_vector_calloc(N);
    gsl_vector* x1 = gsl_vector_calloc(N);

    double lambda;

    /* wypelnienie macierzy A */
    for( int i = 0; i < N; i++ ){
            for( int j = 0; j < N; j++ ) {
                gsl_matrix_set(A, i, j, 1.0 / (sqrt(2.0 + fabs((double)i - (double)j))));
        }
    }

    /* glowna pelta iteracyjna */
    FILE *f1 = fopen("lambdas.txt", "w"); 
    FILE *f2 = fopen("D.txt", "w");
    for (int k = 0; k < N; k++) {
        /* inicjalizacja wektora startowego*/
        gsl_vector_set_all(x, 1.0);
        /* petla wewnetrzna */
        for (int i = 1; i <= IT_MAX; i++) {
            /* xk^i+1 = Akxk^i */
            matrixXvector(x1, A, x, N);
            /* ortogonalizacja G-S */
            for (int j = 0; j < k; j++) {
            	orthogonalization(x1, X, j, N);
            }
            /* przypisanie do lambdy */
            lambda = vectorXvector(x1, x, N) / vectorXvector(x, x, N);
            /* normalizacja wektora x */
            normalizeVector(x, x1, N);
            /* zapis lambda do pliku */
        	if (f1) {
            	fprintf(f1, "%lf ", lambda);
        	}
        }
        if (f1) {fprintf(f1, "\n");}
        /* zapisanie x jako k-tej kolumny X */
        insertColumn(x, X, k, N);
     
    }
    // puts("");

    /* wypelnienie macierzy D, zapis do pliku */
    matrixTXmatrixXmatrix(D, X, A, N);
    if (f2) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(f2, "%g ", gsl_matrix_get(D, i, j));
            }
            fprintf(f2, "\n");
        }
    }   
    
    /* sprzatanie pamieci */
    fclose(f1);
    fclose(f2);
    gsl_matrix_free(A);
    gsl_matrix_free(X);
    gsl_matrix_free(D);
    gsl_vector_free(x);
    gsl_vector_free(x1);

    return 0;
}

/* ------------------------- ---------------------------------------*/

void matrixXvector(gsl_vector *dest, gsl_matrix *matrix, gsl_vector *vector, int n) {
    double value;
    for (int i = 0; i < n; i++) {
        value = 0.0;
        for (int j = 0; j < n; j++) {
            value += gsl_vector_get(vector, j) * gsl_matrix_get(matrix, i, j);
        }
        gsl_vector_set(dest, i, value);
    }
}

double vectorXvector(gsl_vector *vector1, gsl_vector *vector2, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += gsl_vector_get(vector1, i) * gsl_vector_get(vector2, i);
    }
    return result;
}

void normalizeVector(gsl_vector *dest, gsl_vector *src, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += pow(gsl_vector_get(src, i), 2);
    }
    norm = sqrt(norm);
    for (int i = 0; i < n; i++) {
        gsl_vector_set(dest, i, gsl_vector_get(src, i) / norm);
    }
}

void insertColumn(gsl_vector *vector, gsl_matrix *matrix, int k, int n) {
    for (int i = 0; i < n; i++) {
        gsl_matrix_set(matrix, i, k, gsl_vector_get(vector, i));
    }
}

void matrixTXmatrixXmatrix(gsl_matrix *dest, gsl_matrix *matrix1, gsl_matrix *matrix2, int n) {
    gsl_matrix *temp = gsl_matrix_calloc(N, N);
    /* 1. temp = X^t * A */
    double value;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            value = 0.0;
            for (int k = 0; k < n; k++) {
                value += gsl_matrix_get(matrix1, k, i) * gsl_matrix_get(matrix2, k, j);
            }
            gsl_matrix_set(temp, i, j, value);
        }
    }
    /* 2. D = temp * X */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            value = 0.0;
            for (int k = 0; k <n; k++) {
                value += gsl_matrix_get(temp, i, k) * gsl_matrix_get(matrix1, k, j);
            }
            gsl_matrix_set(dest, i, j, value);
        }
    }
    gsl_matrix_free(temp);
}

void orthogonalization(gsl_vector *vector, gsl_matrix *matrix, int k, int n) {
	double value;
	double alpha = 0.0;
	for (int i = 0; i < n; i++) {
		alpha += gsl_vector_get(vector, i) * gsl_matrix_get(matrix, i, k);
	}
	for (int i = 0; i < n; i++) {
		value = gsl_vector_get(vector, i) - alpha * gsl_matrix_get(matrix, i, k);
		gsl_vector_set(vector, i, value);
	}
}
