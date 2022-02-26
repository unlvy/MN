#include <gsl/gsl_linalg.h>
#include <math.h>
#include <stdio.h>

#define XMIN -5.0
#define XMAX 5.0
#define ALPHA 0.0
#define BETA 0.0

void wyznacz_M(double *xw, double *yw, double *m, int n, double alfa, double beta);
double wyznacz_Sx(double *xw, double *yw, double *m, int n, double x);
double function1(double x);
double function2(double x);
void solve(int n, FILE* f1, FILE* f2);

/* ------------------------------------------------------------ */

int main() {

	/* otwieranie plikow */
	FILE* f1 = fopen("fun1_n5.txt", "w");
	FILE* f2 = fopen("fun2_n5.txt", "w");
	FILE* f3 = fopen("fun1_n8.txt", "w");
	FILE* f4 = fopen("fun2_n8.txt", "w");
	FILE* f5 = fopen("fun1_n21.txt", "w");
	FILE* f6 = fopen("fun2_n21.txt", "w");
	FILE* f7 = fopen("fun1_n10_pochodne.txt", "w");

	/* rozwiazanie dla n = 5 */
	solve(5, f1, f2);
	/* rozwiazanie dla n = 8 */
	solve(8, f3, f4);
	/* rozwiazanie dla n = 21 */
	solve(21, f5, f6);

	/* pochodne dla n = 10 */
	int n = 10;
	/* utworzenie i inicjalizacja wektora xw */
	double xw[n];
	double dx = (XMAX - XMIN) / (n - 1);
	for (int i = 0; i < n; i++) {
		xw[i] = XMIN + dx * (double)i;
	}

	/* utworzenie i inicjalizacja wektora yw */
	double yw[n];
	for (int i = 0; i < n; i++) {
		yw[i] = function1(xw[i]);
	}

	/* utworzenie i inicjalizacja wektora m */
	double m[n];
	wyznacz_M(xw, yw, m, n, ALPHA, BETA);

	/* utworzenie i inicjalizacja drugiego wektora pochodnych */
	double mw[n];
	dx = 0.01;
	for (int i = 0; i < n; i++) {
		mw[i] = (function1(xw[i] - dx) - 2 * function1(xw[i]) + function1(xw[i] + dx))
				/ pow(dx, 2);
	}

	/* zapis do pliku */
	if(f7) {
		fprintf(f7, "   x      m   pochodna\n");
		for (int i = 0; i < n; i++) {
	        fprintf(f7, "%.3lf %lf %lf\n", xw[i], m[i], mw[i]);
	    }
	}

	/* zamykanie plikow */
	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);
	fclose(f5);
	fclose(f6);
	fclose(f7);

    return 0;
}

/* ------------------------------------------------------------ */

void wyznacz_M(double *xw, double *yw, double *m, int n, double alfa, double beta) {
    
	/* 1. incjalizacja i utworzenie macierzy A */
    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_matrix_set(A, 0, 0, 1.0);
    gsl_matrix_set(A, n - 1, n - 1, 1.0);

    {
	    double lambda, mi;
	    for (int i = 1; i < n - 1; i++) {

	    	lambda = (xw[i + 1] - xw[i]) / ((xw[i] - xw[i - 1]) + (xw[i + 1] - xw[i]));
	    	mi = 1.0 - lambda;

	        gsl_matrix_set(A, i, i, 2.0);
	        gsl_matrix_set(A,  i, i - 1, mi);
	        gsl_matrix_set(A,  i, i + 1, lambda);
	    }
	}

    /* 2. utworzenie i inicjalizacja wektora d */
    gsl_vector *d = gsl_vector_calloc(n);
    gsl_vector_set(d, 0, alfa);
    gsl_vector_set(d, n - 1, beta);

    {
	    double value;
	    for (int i = 1; i < n - 1; i++) {
	    	value = 6.0 / ((xw[i] - xw[i - 1]) + (xw[i + 1] - xw[i])) *  
	    			((yw[i + 1] - yw[i]) / (xw[i + 1] - xw[i]) - (yw[i] - yw[i - 1]) / (xw[i] - xw[i - 1]));
	    	gsl_vector_set(d, i, value);
	    }
	}

	/* 3. rozwiazanie ukladu */
	gsl_linalg_HH_svx(A, d);

	/* 4. przepisanie wyniku do m */
	for (int i = 0; i < n; i++) {
        m[i] = gsl_vector_get(d, i);
    }
    
    /* 5. zwolnienie pamieci */
    gsl_vector_free(d);
    gsl_matrix_free(A);
}

double wyznacz_Sx(double *xw, double *yw, double *m, int n, double x) {

	/* 1. potrzebne zmienne */
    double result;
    int i;

    /* 2. znalezienie pierwszego podprzedzialu */
    for (i = 1; i < n; i++) {
        if (xw[i - 1] <= x && x <= xw[i]) {
            break;
        }
    }
    
    /* 3. wyznaczenie Sx */
    double h = xw[i] - xw[i - 1];
    double A = ((yw[i] - yw[i - 1]) / h) - (h / 6) * (m[i] - m[i - 1]);
    double B = yw[i - 1] - (m[i - 1] * h * h) / 6;

    result = m[i - 1] * (pow((xw[i] - x), 3) / (6.0 * h)) +
             m[i] * (pow((x - xw[i - 1]), 3) / (6.0 * h)) +
             A * (x - xw[i - 1]) + B;

    return result;
}

double function1(double x) {
	return 1.0 / (1.0 + pow(x, 2));
}

double function2(double x) {
	return cos(2 * x); 
}

void solve(int n, FILE* f1, FILE* f2) {
	
	/* 1. interpolacja f1 */
	/* 1.1 utworzenie i inicjalizacja wektora xw */
	double xw[n];
	double dx = (XMAX - XMIN) / (n - 1);
	for (int i = 0; i < n; i++) {
		xw[i] = XMIN + dx * (double)i;
	}

	/* 1.2. utworzenie i inicjalizacja wektora yw */
	double yw[n];
	for (int i = 0; i < n; i++) {
		yw[i] = function1(xw[i]);
	}

	/* 1.3. utworzenie i inicjalizacja wektora m */
	double m[n];
	wyznacz_M(xw, yw, m, n, ALPHA, BETA);
	
	/* 1.4. zapis wynikow do pliku */
	if (f1) {
		for (int i = 0; i <= 1000; i++) {
	        fprintf(f1, "%lf ", wyznacz_Sx(xw, yw, m, n, XMIN + i * 0.01));
	    }
	}

	/* 2. interpolacja f2 */
	/* 2.1. inicjalizacja wektora yw */
	for (int i = 0; i < n; i++) {
		yw[i] = function2(xw[i]);
	}

	/* 2.2. inicjalizacja wektora m */
	wyznacz_M(xw, yw, m, n, ALPHA, BETA);

	/* 2.3. zapis wynikow do pliku */
	if (f2) {
		for (int i = 0; i <= 1000; i++) {
	        fprintf(f2, "%lf ", wyznacz_Sx(xw, yw, m, n, XMIN + i * 0.01));
	    }
	}
}
