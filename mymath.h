#include "float.h"
#include <math.h>
#include <stdlib.h>
#include "stdio.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

double *SFactors(int N);
double *SSFactors(int N);
double *SSSFactors(int N);


void vector_matrix_mult(double *result, double **m, double *v, int N, int skip);
void createSN(double **un, double N);
void createTN(double **un, double N);
void createUN(double **un, double N);
void t2(double *S, double t2[3][3]);
void u2(double *S, double u2[3][3]);
