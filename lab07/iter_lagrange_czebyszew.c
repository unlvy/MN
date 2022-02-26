#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define n 20

double functionValue(double x);
double lagrangeInterpolation(double* vecX, double* vecY, double x);
double Czebyszew(int m);

int main() {

	/* utworzenie wektorow argumentow i wartosci */
	double vecX[n + 1];
	double vecY[n + 1];
	double vecXCzebyszew[n + 1];
	double vecYCzebyszew[n + 1];

	/* krok, o jaki beda zwiekszac sie wartosci dla vecX */
	double step = (10.0 / n);

	/* wypelnianie wektorow */
	for (int i = 0; i <= n; i++) {
		vecX[i] = -5.0 + step * i;
		vecY[i] = functionValue(vecX[i]);
		vecXCzebyszew[i] = Czebyszew(i);
		vecYCzebyszew[i] = functionValue(vecXCzebyszew[i]);
	}

	/* glowna petla iteracyjna z zapisami wartosci do plikow */
	FILE *f1 = fopen("resultsLagrange.txt", "w");
	FILE *f2 = fopen("resultsCzebyszew.txt", "w");
	if(f1 && f2) {
		/* wewnetrzna petla iteracyjna po wartosciach */
		for (double x = -5.0; x <= 5.0; x += 0.01) {
			/* zapois do pliku wartosci interpolacji Lagrangea dla wezlow rownoodleglych */
			fprintf(f1, "%lf ", lagrangeInterpolation(vecX, vecY, x));
			/* zapis do pliku wartosci interpolacji Lagrangea dla wezlow Czebyszewa */
			fprintf(f2, "%lf ", lagrangeInterpolation(vecXCzebyszew, vecYCzebyszew, x));
		}
	}
	fclose(f1);
	fclose(f2);

	return 0;
}

double functionValue(double x) {
	return pow(M_E, -pow(x, 2));
}

double lagrangeInterpolation(double* vecX, double* vecY, double x) {

	double polynomialSum = 0.0;
	    for (int i = 0; i <= n; i++) {
	        double multipler = 1.0;
	        for (int j = 0; j <= n; j++) {
	            if (j != i) {
	                multipler *= ((x - vecX[j]) / (vecX[i] - vecX[j]));
	            }
	        }

	        polynomialSum += vecY[i] * multipler;
	    }
	    return polynomialSum;
}

double Czebyszew(int m) {
	return (10.0 * cos(M_PI * (2.0 * (double)m + 1.0) / (2.0 * (double)n + 2.0))) / 2.0;
}
