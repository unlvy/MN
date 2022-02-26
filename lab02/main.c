#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#define N 4

int main () {

	/* utworzenie i inicjalizacja macierzy A*/
	gsl_matrix *A = gsl_matrix_calloc(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			gsl_matrix_set(A, i, j, 1 / (i + j + 2.0));
		}
	}

	/* utworzenie wektora p i inta signum */
	gsl_permutation *p = gsl_permutation_calloc(N);
	int signum;	

	/* rozklad LU macierzy */
	gsl_linalg_LU_decomp(A, p, &signum);

	/* zapis do pliku */
	FILE *f = fopen("out_1.txt", "w");
	    if (f) {
		fprintf(f, "Elementy diagonalne macierzy U:\n");
		double determinant = 1.0;
		for (int i = 0; i < N; i ++) {
			fprintf(f, "%lf\t", gsl_matrix_get(A, i, i));
			determinant *= gsl_matrix_get(A, i, i);
		}
		fprintf(f, "\nWyznacznik = %.12lf\n", determinant);
	}
	fclose(f);

	/* utworzenie wektorow b i x */
	gsl_vector *b = gsl_vector_calloc(N);
	gsl_vector *x = gsl_vector_calloc(N);

	/* utworzenie macierzy B (A^(-1)) */
	gsl_matrix *B = gsl_matrix_calloc(N, N);

	/* wyznaczenie macierzy odwrotnej*/
	for (int i = 0; i < N; i++) {

		/* wypelnienie b aktualnymi wartosciami */
		for (int j = 0; j < N; j++) {
			gsl_vector_set(b, j, 0.0);
		}
		gsl_vector_set(b, i, 1.0);

		/* wywolanie procedury */
		gsl_linalg_LU_solve(A, p, b, x);

		/* wypelnienie kolumny macierzy B aktualnymi wartosciami wektora x */
		for (int k = 0; k < N; k++) {
			gsl_matrix_set(B, k, i, gsl_vector_get(x, k));
		}
	}

	/* zapis do pliku */
	f = fopen("out_2.txt", "w");
    if (f) {
		fprintf(f, "Macierz odwrotna:\n");
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				fprintf(f, "%lf\t", gsl_matrix_get(B, i, j));
			}
		fprintf(f, "\n");
		}
	}
	fclose(f);
	

	/* utworzenie macierzy C (do iloczynu A * A^(-1)) */
	gsl_matrix *C = gsl_matrix_calloc(N, N);

	/* zapisanie do A pierwotnych wartosci */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			gsl_matrix_set(A, i, j, 1 / (i + j + 2.0));
		}
	}

	/* obliczenie ilorazu A * A^(-1) */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double value = 0;

			for (int k = 0; k < N; k++) {
				value += gsl_matrix_get(A, i, k) * gsl_matrix_get(B, k, j);
			}

			gsl_matrix_set(C, i, j, value);
		}
	}

	/* zapis do pliku */
	f = fopen("out_3.txt", "w");
    if (f) {
		fprintf(f, "Iloczyn A * A^(-1)\n");
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				fprintf(f, "%1.5g\t", gsl_matrix_get(C, i, j));
			}
		fprintf(f, "\n");
		}
	}
	fclose(f);

	/* obliczanie wskaznika uwarunkowania */
	double maxA;
	double maxB;
	for (int i = 0; i < N; i++)	{
		for (int j = 0; j < N; j++) {

			if (maxA < fabs(gsl_matrix_get(A, i, j))) {
				maxA = fabs(gsl_matrix_get(A, i, j));
			}

			if (maxB < fabs(gsl_matrix_get(B, i, j))) {
				maxB = fabs(gsl_matrix_get(B, i, j));
			}

		}
	}
	double cond = maxA * maxB;

	/* zapis do pliku */
	f = fopen("out_4.txt", "w");
    if (f) {
		fprintf(f, "cond =  %lf\n", cond);
	}
	fclose(f);
	
	/* zwalnianie pamieci */
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(C);
	gsl_permutation_free(p);
	gsl_vector_free(x);
	gsl_vector_free(b);

	return 0;
}


