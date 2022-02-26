#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define START 0.0
#define STOP 1.0
#define LOWERCASE_N 8

double f(double x);
void createMatrix(double*** matrix, int n);
void deleteMatrix(double*** matrix, int n);
double Simpson(double step, int N);
double Milne(double step, int N);
void solve(FILE* file, bool version); // true - simpson, false - milne

/* --------------------- */

int main() {

	FILE* fileSimpson = fopen("SimpsonResults.txt", "w");
	FILE* fileMilne = fopen("MilneResults.txt", "w");

	if (fileSimpson) {
		solve(fileSimpson, true);
	}
	if (fileMilne) {
		solve(fileMilne, false);
	}

	fclose(fileSimpson);
	fclose(fileMilne);

}

/* --------------------- */

double f(double x) {
	return log(pow(x, 3) + 3 * pow(x, 2) + x + 0.1) * sin(18 * x);
}

void createMatrix(double*** matrix, int n) {
	*matrix = (double**)malloc(n * sizeof(double*));
	for (int i = 0; i < n; i++) {
		(*matrix)[i] = (double*)malloc(n * sizeof(double));
	}
}

void deleteMatrix(double*** matrix, int n) {
	for (int i = 0; i < n; i++) {
		free((*matrix)[i]);
	}
	free(*matrix);
}

double Simpson(double step, int N) {
	double result = 0.0;
	for (int i = 0; i <= (N / 2) - 1; i++) {
		result += (step / 3.0) * ( f( (2 * i) * step) + 4 * f( (2 * i + 1) * step) 
			    + f( (2 * i + 2) * step) );
	} 
	return result;
}

double Milne(double step, int N) {
	double result = 0.0;
	for (int i = 0; i <= (N / 4) - 1; i++) {
		result += ((4.0 * step) / 90.0) * (7 * f((4 * i) * step) + 32 * 
				  f((4 * i + 1) * step) + 12 * f((4 * i + 2) * step) + 32 * f((4 * i + 3) * step)
		 		  + 7 * f((4 * i + 4) * step));
	} 
	return result;
}

void solve(FILE* file, bool version) {

	double** D;
	createMatrix(&D, LOWERCASE_N + 1);
	int N;
	double step;

	for (int i = 0; i <= LOWERCASE_N; i++) {
		N = pow(2.0, i + (version ? 1 : 2));
		step = (STOP - START) / (double)N;
		D[i][0] = (version ? Simpson(step, N) : Milne(step, N));
	}
	puts("");

	for (int j = 1; j <= LOWERCASE_N; j++) {
		for (int i = j; i <= LOWERCASE_N; i++) {
			D[i][j] = ((pow(4.0, j) * D[i][j - 1]) - D[i - 1][j - 1]) / (pow(4.0, j) - 1.0);
		}
	}

	fprintf(file, "Pierwsza Kolumna:\n");
	for (int i = 0; i <= LOWERCASE_N; i++) {
		fprintf(file, "D[%d][0] = %.12lf\n", i, D[i][0]);
	}

	fprintf(file, "Diagonala:\n");
	for (int i = 0; i <= LOWERCASE_N; i++) {
		fprintf(file, "D[%d][%d] = %.12lf\n", i, i, D[i][i]);
	}

	deleteMatrix(&D, LOWERCASE_N + 1);
}
