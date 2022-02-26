#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define H 0.1
#define DELTA 0.0001

double f(double x, double y);
double derivative(double x, double y, double ex, double ey);
double norm(double x0, double y0, double x1, double y1);
void solve(double epsilon, FILE* f);

/* ----------------------------------------- */

int main() {

	FILE* f1 = fopen("eps1.dat", "w");
	FILE* f2 = fopen("eps2.dat", "w");

	solve(0.01, f1);
	solve(0.005, f2);

	fclose(f1);
	fclose(f2);

	return 0;
}

/* ----------------------------------------- */

double f(double x, double y) {
	return (5.0 / 2.0) * pow(pow(x, 2) - y, 2) + pow(1.0 - x, 2);
}

double derivative(double x, double y, double ex, double ey) {
	return (f(x + ex * DELTA, y + ey * DELTA) - f(x - ex * DELTA, y - ey * DELTA)) / (2.0 * DELTA);
}

double norm(double x0, double y0, double x1, double y1) {
	return sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
}

void solve(double epsilon, FILE* f) {
	double x1, y1;
	double x0 = -0.75;
	double y0 = 1.75;

	int i;
	for (i = 1; i <= 1000; i++) {
		
		x1 = x0 - H * derivative(x0, y0, 1.0, 0.0);
		y1 = y0 - H * derivative(x0, y0, 0.0, 1.0);

		if (f) {
			fprintf(f, "%lf %lf\n", x1, y1);
		}

		if (norm(x0, y0, x1, y1) < epsilon) {
			break;
		}

		x0 = x1;
		y0 = y1;
	}
	printf("%lf %lf %d\n", x1, y1, i);
}
