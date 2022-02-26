#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"
#include "gauleg.c"
#include "gaulag.c"
#include "gauher.c"
#include "gammln.c"

#define INTEGRAL_1 ((float)(M_PI / 3.0f))
#define INTEGRAL_2 -0.8700577
#define INTEGRAL_3 (2.0f / 13.0f)

float f1(float x);
float f2a(float x);
float f2b(float x);
float f3(float x);

/* --------------------- */

int main() {

	FILE* file = fopen("out.dat", "w");

	float* x;
	float* w;

	float sum, error;

	/* --------------- 2.1 --------------- */
	for (int n = 2; n <= 100; n++) {

		sum = 0.0f;
		x = vector(1, n);
		w = vector(1, n);

		gauleg(1.0f, 2.0f, x, w, n);

		for (int i = 1; i <= n; i++) {
			sum += w[i] * f1(x[i]);
		}

		error = fabs(INTEGRAL_1 - sum);

		if(file) {
			fprintf(file, "%d %f\n", n, error);
		}
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	if(file) { fprintf(file, "\n\n"); }

	/* --------------- 2.2a --------------- */
	for (int n = 100; n <= 100; n+=2) {

		sum = 0.0f;
		x = vector(1, n);
		w = vector(1, n);

		gauher(x, w, n);

		for (int i = 1; i <= n; i++) {
			sum += w[i] * f2a(x[i]);
			printf("%f ", x[i]);
		}
		printf("\n");


		error = fabs(INTEGRAL_2 - sum);

		if(file) {
			fprintf(file, "%d %f\n", n, error);
		}
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	if(file) { fprintf(file, "\n\n"); }

	/* --------------- 2.2b --------------- */
	for (int n = 100; n <= 100; n++) {

		sum = 0.0f;
		x = vector(1, n);
		w = vector(1, n);

		gauleg(0.0f, 5.0f, x, w, n);

		for (int i = 1; i <= n; i++) {
			sum += w[i] * f2b(x[i]);
			printf("%f ", x[i]);
		}
		printf("\n");
		
		error = fabs(INTEGRAL_2 - sum);

		if(file) {
			fprintf(file, "%d %f\n\n", n, error);
		}
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	if(file) { fprintf(file, "\n"); }

	/* --------------- 2.3 --------------- */
	for (int n = 2; n <= 20; n++) {

		sum = 0.0f;
		x = vector(1, n);
		w = vector(1, n);

		gaulag(x, w, n, 0.0f);

		for (int i = 1; i <= n; i++) {
			sum += w[i] * f3(x[i]);
		}

		error = fabs(INTEGRAL_3 - sum);

		if(file) {
			fprintf(file, "%d %f\n", n, error);
		}
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	if(file) { fprintf(file, "\n\n"); }


	fclose(file);
	return 0;

}

/* --------------------- */

float f1(float x) {
	return 1.0f / (x * ( sqrt( pow( x, 2 ) - 1) ) );
}

float f2a(float x) {
	return 0.5f * log( fabs( x ) );
}

float f2b(float x) {
	return log( x ) * exp( - pow( x, 2 ));
}

float f3(float x) {
	return sin( 2 * x ) * exp( -2 * x);
}