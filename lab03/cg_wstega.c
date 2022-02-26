#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#define N 10000
#define M 5

#define max(X, Y) ((X) > (Y) ? (X) : (Y))
#define min(X, Y) ((X) < (Y) ? (X) : (Y))
#define abs(X) ((X) > 0 ? (X) : -(X))

/* do iloczynu skalarnego wektorow */
double vecTxvec(gsl_vector* v1, gsl_vector* v2);
/* do mnozenia macierz wektor */
void matrixxvec(gsl_matrix* A, gsl_vector* v, gsl_vector* _Avk);
/* do mnozenia wektora przez liczbe */
void doublexvec(double multipler, gsl_vector* source, gsl_vector* result);
/* do dodawania wektorow */
void vecplusvec(gsl_vector* s1, gsl_vector* s2, gsl_vector* result);

int main() {

	clock_t start = clock();
	/* utworzenie macierzy A */
	gsl_matrix *A = gsl_matrix_calloc(N, N);

	/* inicjalizacja macierzy A */
	double value;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			value = (abs(i - j) > M) ? 0.0 : 1.0 / (1.0 + (double)abs(i - j));
			gsl_matrix_set(A, i, j, value);	
		}
	}

	/* utworzenie wektora wyrazow wolnych b */
	gsl_vector *b = gsl_vector_calloc(N);

	/* inicjalizacja wektora b */
	for (int i = 0; i < N; i++) {
		gsl_vector_set(b, i, (double)i + 1.0);
	}

	/* utworzenie wektora startowego x */
	gsl_vector *x = gsl_vector_calloc(N);

	/* inicjalizacja wektora x */
	for (int i = 0; i < N; i++) {
		gsl_vector_set(x, i, 0.0);
	}

	/* metoda CG */
	{
		/* utworzenie wektorow r, v, Avk, support i wartosci alfa, beta */
		gsl_vector *r = gsl_vector_calloc(N);
		gsl_vector *v = gsl_vector_calloc(N);
		gsl_vector *_Avk = gsl_vector_calloc(N);
		gsl_vector *support = gsl_vector_calloc(N);

		double alfa;
		double beta;

		/* inicjalizacja (jako ze x = 0, wypelniam tak jak b) */
		for (int i = 0; i < N; i++) {
			gsl_vector_set(r, i, (double)i + 1.0);
			gsl_vector_set(v, i, (double)i + 1.0);
		}

		/* petla iteracyjna CG*/
		int k = 0;
		double vkAvk;
		double _rtr = vecTxvec(r, r);
		// FILE* f = fopen("outxd.txt", "w");
		// if (f) {
  //   		fprintf(f, "iteracja\t ||rk||2\t alfak\t\t betak\t\t ||xk||2\n");
  //   		fprintf(f, "%8d\t %7lf\t\t\t\t\t %7lf\t\n",
  //   					k, sqrt(_rtr), sqrt(vecTxvec(x, x)));
  //   	}
		while (sqrt(_rtr) > 1e-6) {
			matrixxvec(A, v, _Avk);

			k++;

			/* 1. nadpisanie wartosci alfa*/
			alfa = _rtr / vecTxvec(v, _Avk);

			/* 2. nadpisanie wektora x */
			doublexvec(alfa, v, support);
			vecplusvec(x, support, x);

			/* 3. nadpisanie wektora r */
			doublexvec(-alfa, _Avk, support);
			vecplusvec(r, support, r);

			/* 4. nadpisanie wartosci beta */
			beta = vecTxvec(r, r) / _rtr;

			/* 5. nadpisanie wektora v */
			doublexvec(beta, v, support);
			vecplusvec(r, support, v);

			/* 6. aktualizacja wartosci _rtr */
		    _rtr = vecTxvec(r, r);

			// /* 7. zapis do pliku */
			// if (f) {
   //  			fprintf(f, "%lf ", sqrt(vecTxvec(x, x)));
		 //    }

		}
		clock_t end = clock();
		float seconds = (float)(end - start) / CLOCKS_PER_SEC;
		printf("%f", seconds);
		// fclose(f);

		/* zwalnianie pamieci */
		gsl_vector_free(r);
		gsl_vector_free(v);
		gsl_vector_free(_Avk);
		gsl_vector_free(support);
	}

	/* zwalnianie pamieci */
	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_vector_free(x);
	

	return 0;
}

double vecTxvec(gsl_vector* v1, gsl_vector* v2) {
	double result = 0.0;
	for(int i = 0; i < N; i++) {
		result += gsl_vector_get(v1, i) * gsl_vector_get(v2, i);
	}
	return result;
}

void matrixxvec(gsl_matrix* A, gsl_vector* v, gsl_vector* result) {
	double value, jmin, jmax;
	for (int i = 0; i < N; i++) {
		jmin=max(0, i - M);
		jmax=min(i + M, N - 1);
		value = 0.0;
		for (int j = jmin; j <= jmax; j++) {
			value += gsl_matrix_get(A, i, j) * gsl_vector_get(v, j);
		}
		gsl_vector_set(result, i, value);
	}
}

void doublexvec(double multipler, gsl_vector* source, gsl_vector* result) {
	double value;
	for (int i = 0; i < N; i++) {
		value = multipler * gsl_vector_get(source, i);
		gsl_vector_set(result, i, value);
	}
}

void vecplusvec(gsl_vector* s1, gsl_vector* s2, gsl_vector* result) {
	double value;
	for (int i = 0; i < N; i++) {
		value = gsl_vector_get(s1, i) + gsl_vector_get(s2, i);
		gsl_vector_set(result, i, value);
	}
}
