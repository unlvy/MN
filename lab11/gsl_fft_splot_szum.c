#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define T 1.0
#define OMEGA (2.0 * M_PI / T)
#define SIGMA 0.02 
#define TMAX (3.0 * T)

/** oblicza wartosc f0 */
double f0(double x);
/** oblicza wartosc g */
double g(double x);
/** oblicza losowa delte */
double delta();
/** przeprowadza odszumianie */
void solve(int k, FILE* f);

/* --------------------- */

int main() {

	FILE* f1 = fopen("k8.dat", "w");
	FILE* f2 = fopen("k10.dat", "w");
	FILE* f3 = fopen("k12.dat", "w");

	if (f1) {
		solve(8, f1);	
	}
	if (f2) {
		solve(10, f2);	
	}
	if (f3) {
		solve(12, f3);	
	}

	fclose(f1);
	fclose(f2);
	fclose(f3);

	return 0;

}

/* --------------------- */

double f0(double x) {
	return sin(OMEGA * x) + sin(2.0 * OMEGA * x) + sin(3.0 * OMEGA * x);
}

double g(double x) {
	return 1.0 / (SIGMA * sqrt(2.0 * M_PI)) * exp(-(pow(x, 2) / (2.0 * pow(SIGMA, 2))));
}

double delta() {
	return rand() / (RAND_MAX + 1.0) - 0.5;
}

void solve(int k, FILE* f) {

	int n = pow(2, k);
	double dt = TMAX / (double)n;
	double ti, a1, a2, b1, b2;

	/* wszystkie niezbedne tablice */
	double f_data[2 * n];
	double g1_data[2 * n];
	double g2_data[2 * n];
	double g_data[2 * n];

	/* wypelnianie tablic, zapis do pliku zaszumionych wartosci */
	for (int i = 0; i < n; i++) {
		ti = dt * i;
		f_data[2 * i] = f0(ti) + delta();
		f_data[2 * i + 1] = 0.0;
		g1_data[2 * i] = g(ti);
		g1_data[2 * i + 1] = 0.0;
		g2_data[2 * i] = g(ti);
		g2_data[2 * i + 1] = 0.0;
		fprintf(f, "%lf %lf\n", ti, f_data[2 * i]);
	}

	fprintf(f, "\n\n");

	/* liczenie FFT */
	gsl_fft_complex_radix2_forward(f_data, 1, n);
	gsl_fft_complex_radix2_forward(g1_data, 1, n);
	gsl_fft_complex_radix2_backward(g2_data, 1, n);

	/* g(k) = FFT{g(t)} + FFT^(-1){g(t)} */
	for (int i = 0; i < 2 * n; i++) {
		g_data[i] = g1_data[i] + g2_data[i];
	}

	/* obliczenie splotu */
	for (int i = 0; i < n; i++) {
		a1 = f_data[2 * i];
		b1 = f_data[2 * i + 1];
		a2 = g_data[2 * i];
		b2 = g_data[2 * i + 1];
		f_data[2 * i] = a1 * a2 - b1 * b2;
		f_data[2 * i + 1] = a1 * b2 + b1 * a2;
	}

	/* transformacja odwrotna */
	gsl_fft_complex_radix2_backward(f_data, 1, n);

    /* poszukiwania maximum */
	double max = f_data[0];
	for (int i = 0; i < n; i++) {
		if (f_data[2 * i] > max) {
			max = f_data[2 * i];
		}
	}

	/* zapis do pliku unormalizowanych wartosci */
	for (int i = 0; i < n; i++) {
		ti = dt * i;
		fprintf(f, "%lf %lf\n", ti, f_data[2 * i] * 2.5 / max);
	}

}
