#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <math.h>

#define CAPITAL_N 1.0
#define L 10.0
#define N 200

double ro(double x, int alpha);
double delta(int i, int j);

int main() {

	double value, xi;
	const double deltaX = L / (N + 1.0);

	/* utworzenie macierzy A */
	gsl_matrix *A = gsl_matrix_calloc(N, N);
	/* utworzenie macierzy B */
	gsl_matrix *B = gsl_matrix_calloc(N, N);
	/* utworzenie wektora wlasnosci wlasnych eval */
	gsl_vector *eval = gsl_vector_calloc(N);
	/* utowrzenie macierzy evec */
	gsl_matrix *evec = gsl_matrix_alloc(N, N);
	/* utworzenie wektora pomocniczego w */
	gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(N);


	/* przygotowanie do petli */
	FILE* f1 = fopen("out1.txt", "w");
	FILE* f2 = fopen("out2.txt", "w");

	/* petla iteracyjna z alfa */
	for (int alpha = 0; alpha <= 100; alpha += 2) {

		/* wypelnienie macierzy A */
		for (int i = 0; i < N; i++) {
			xi = (-L / 2.0) + deltaX * (i + 1.0);
			for (int j = 0; j < N; j++) {
				value = (-delta(i, j + 1) + 2 * delta(i, j)
						 - delta(i, j - 1)) / pow((double)(L / (N + 1.0)), 2);
				gsl_matrix_set(A, i, j, value);
			}
		}


		/* wypelnienie macierzy B */
		for (int i = 0; i < N; i++) {
			xi = (-L / 2.0) + (deltaX * (i + 1.0));
			value = ro(xi, alpha) / CAPITAL_N;
			gsl_matrix_set(B, i, i, value);
		}	

		/* rozwiazanie rownania */
		gsl_eigen_gensymmv(A, B, eval, evec, w);

		/* sortowanie wartosci i wektorow */
		gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		/* zapisanie do pliku wartosci pierwiastkow z 6 kolejnych najmniejszych wartosci wlasnych */
		// if (f1) {
		// 	fprintf(f1, "%1.5lf\t", sqrt(gsl_vector_get(eval, 5)));
		// 	// fprintf(f1, "%d ", alpha);
		// 	// for (int i = 0; i < 6; i++) {
		//  //    		fprintf(f1, "%1.10lf\t", sqrt(gsl_vector_get(eval, i)));
	 //  //   		}
	 //  //   		fprintf(f1, "\n");
  //   	}

  //   	 zapisanie do pliku wektorow wlasnych odpowiadajacych 6 najnizszym wartosciom 
    	if (alpha == 100 && f2) {
    		for (int i = 0; i < N; i++) {
	    		fprintf(f2, "%lf ", gsl_matrix_get(evec, i, 5));
	    	}
    		// for (int j = 0; j < N; j++) {
    		// 	for (int i = 0; i < 6; i++) {
    		// 		fprintf(f2, "%1.5lf  ", gsl_matrix_get(evec, j, i));
    		// 	}
    		// 	fprintf(f2, "\n");
    		// }
    	// 	// fprintf(f2, "\n\n\n");
    	// }

	}
	// if (f2) {
	// for (int i = 0; i < N; i++) {
	// 	xi = (-L / 2.0) + deltaX * (i + 1.0);
		
	// 		fprintf(f2, "%lf ", xi);
		
	// }
}
	fclose(f1);
	fclose(f2);	


	/* zwalnianie pamieci */
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_eigen_gensymmv_free(w);

	return 0;
}


double ro(double x, int alpha) {
	return (double)(1.0 + 4.0 * alpha * pow(x, 2));
}

double delta(int i, int j) {
	return i == j ? 1.0 : 0.0;
}
