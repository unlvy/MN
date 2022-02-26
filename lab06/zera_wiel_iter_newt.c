#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IT_MAX 30
#define N 5

double licz_r(double* a, double* b, int n, double x0);

int main() {

	/* utworzenie wektorow a, b, c */
	double* a = (double*)malloc((N + 1) * sizeof(double));
	double* b = (double*)malloc((N + 1) * sizeof(double));
	double* c = (double*)malloc(N * sizeof(double));

	/* proces wyznaczania zer */
	{
		/* inicjalizacja wektora danych */
		a[0] = 240.0;
		a[1] = -196.0;
		a[2] = -92.0;
		a[3] = 33.0;
		a[4] = 14.0;
		a[5] = 1.0;

		/* petla po kolejnych zerach wielomianu */
		int n;
		double Rj, Rj_prim, x0, x1;
		FILE *f = fopen("results.txt", "w");
		if (f) {
				fprintf(f, " L     it         xit         Rit           R`it\n");
			}
		for (int L = 1; L <= N; L++) {

			/* ustalenie aktualnego stopnia wielomianu */
			n = N - L + 1;

			/* inicjalizacja wzoru iteracyjnego */
			x0 = 0.0;

			/* petla wewnetrzna */
			for (int it = 1; it < IT_MAX; it++) {

				/* wyznaczenie Rj */
				Rj = licz_r(a, b, n, x0);

				/* wyznaczenie Rj_prim */
				Rj_prim = licz_r(b, c, n - 1, x0);

				x1 = x0 - Rj / Rj_prim;

				/* zapis do pliku */
				if (f) {
					fprintf(f, "%2d    %2d    %10lf    %10lf    %10lf\n", L, it, x1, Rj, Rj_prim);
				}

				/* warunek wczesniejszego zakonczenia petli */
				if(fabs(x1-x0) < 1.0e-7) break; 

				/* zachowanie nowego przyblizenia */
				x0 = x1;

			}

			/* redukcja stopnia wielomianu o 1 */
			for (int i = 0; i <= n - 1; i++) a[i] = b[i];

			if (f) {
				fprintf(f,"\n");
			}
		}
		fclose(f);
	}

	/* zwalnianie pamieci */
	free(a);
	free(b);
	free(c);

	return 0;
}

double licz_r(double* a, double* b, int n, double x0) {
	b[n] = 0.0;
	for (int k = n - 1; k >= 0; k--) {
		b[k] = a[k + 1] + x0 * b[k + 1];
	}
	return a[0] + x0 * b[0];
}
