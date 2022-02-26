#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MU 4.0
#define DELTA 3.0

/** funkcja generatora 1 */
double gen1();

/** funkcja generatora 2 */
double gen2();

/** funkcja generatora 3 */
double gen3();

/** zwraca liczbe rozkladu trojkatnego */
double triangularDistribution();

/** oblicza srednia */
double calculateMu(double* results, int n);

/** oblicza odchylenie standardowe */
double calculateSigma(double* results, int n);

/** oblicza ilosci wystapien w danym przedziale */
void calculateHValues(double* results, int* h, int n, int range, double a, double b);

/** oblicza dystrybuante dla rozkladu trojkatnego */
double F(double x);

/* ------------------------------------------------------------ */

int main() {

	/* 2.1 */
	int n = 10000;
	double gen1Results[n];
	double gen2Results[n];

	for (int i = 0; i < n; i++) {
		gen1Results[i] = gen1();
		gen2Results[i] = gen2();
	}

	/* zapis do plikow kolumn xi i xi+1 */
	FILE* f1 = fopen("U.dat", "w");
	if (f1) {
		for (int i = 0; i < n - 1; i++) {
			fprintf(f1, "%lf %lf\n", gen1Results[i], gen1Results[i+1]);
		}
		fprintf(f1, "\n\n");
		for (int i = 0; i < n - 1; i++) {
			fprintf(f1, "%lf %lf\n", gen2Results[i], gen2Results[i+1]);
		}
	}
	fclose(f1);

	/* obliczenie sredniej i odchylenia std dla gen1 i gen2 */
	printf("gen1: mu = %lf, sigma = %lf\n", calculateMu(gen1Results, n), calculateSigma(gen1Results, n));
	printf("gen2: mu = %lf, sigma = %lf\n", calculateMu(gen2Results, n), calculateSigma(gen2Results, n));

	/* histogramy dla gen1 i gen2 */
	int gen1h[12] = {0};
	int gen2h[12] = {0};
	double a = 0.0;
	double b = 1.0;
	double step = (b - a) / 12;
	calculateHValues(gen1Results, gen1h, n, 12, a, b);
	calculateHValues(gen2Results, gen2h, n, 12, a, b);
	/* zapis do pliku */
	FILE* f2 = fopen("U_hist.dat", "w");
	if (f2) {
		for (int i = 0; i < 12; i++) {
			fprintf(f1, "%lf %lf\n",  (a + i * step) + step / 2.0, (double)gen1h[i] / n);
		}
		fprintf(f2, "\n\n");
		for (int i = 0; i < 12; i++) {
			fprintf(f2, "%lf %lf\n", (a + i * step) + step / 2.0, (double)gen2h[i] / n);
		}
	}
	fclose(f2);

	/* 2.2 */
	n = 1000;
	double triangularDistributionResults[n];
	int tdh[10] = {0};

	for (int i = 0; i < n; i++) {
		/* 2.2.1 */
		triangularDistributionResults[i] = triangularDistribution();
	}

	/* 2.2.2 */
	a = MU - DELTA;
	b = MU + DELTA;
	step =  (b - a) / 10.0;
	calculateHValues(triangularDistributionResults, tdh, n, 10, a, b);


	/* 2.2.3 */
	FILE* f3 = fopen("T_hist.dat", "w");
	double X2 = 0.0;
	double pi;
	for (int i = 0; i < 10; i++) {
		pi = F(a + (i + 1) * step) - F(a + i * step);
		X2 += pow(tdh[i] - n * pi, 2) / (pi * n);

		/* zapis do pliku na potrzeby histogramu*/
		if (f3) {
			fprintf(f3, "%lf %lf %lf\n",  (a + i * step) + step / 2.0, (double)tdh[i] / n, pi);
		}		
	}
	printf("Test dla gen trojkatnego: X2 = %lf\n", X2);
	fclose(f3);

	return 0;
}

/* ------------------------------------------------------------ */

double gen1() {
	// static long long int x = 10;
	// int a = 123;
	// int c = 1;
	// long int m = pow(2, 15);
	// x = (a * x + c) % m;
	// return x / (m + 1.0);
	return rand() / (RAND_MAX + 1.0);
}

double gen2() {
	static long long int x = 10;
	int a = 69069;
	int c = 1;
	long int m = pow(2, 32);
	x = (a * x + c) % m;
	return x / (m + 1.0);
}

double gen3() {
	static long long int x = 10;
	int a = 69069;
	int c = 1;
	long int m = pow(2, 32);
	x = (a * x + c) % m;
	return x / (m + 1.0);
}

double calculateMu(double* results, int n) {
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += results[i];
	}
	return sum / n;
}

double calculateSigma(double* results, int n) {
	double mu = calculateMu(results, n);
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += pow(results[i] - mu, 2);
	}
	return sqrt(sum / n);
}

void calculateHValues(double* results, int* h, int n, int range, double a, double b) {
	double step =  (b - a) / range;
	int j; 
	for (int i = 0; i < n; i++) {
		j = (results[i] - a) / step;
		h[j]++;
	}

}

double triangularDistribution() {
	return MU + ( gen3() + gen3() - 1.0 ) * DELTA;
}

double F(double x) {
	if (x <= MU) {
		return  ( -(  1.0 / pow(DELTA, 2) ) * ( ( -pow(x, 2) / 2.0 + MU * x) ) ) + ( x / DELTA );
	} else {
		return -( 1.0 / pow (DELTA, 2) ) * ( ( pow(x, 2) / 2.0 ) - MU * x + pow(MU, 2) ) + x / DELTA;
	}	
}
