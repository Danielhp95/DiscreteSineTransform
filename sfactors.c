#include "mymath.h"

double *SFactors(int N) {
  int den, nom, index = 1;
  double *res = (double *) calloc(N, sizeof(double)); res[0] = 1.0;

  for (den = 2; den <= N; den <<= 1) {
    for (nom = 1; nom <= (den/2); nom += 2) {
        res[index++] = sin(((double)nom) * PI / (double)den);
    }
  }
  return res;
}

double *SSFactors(int N) {
  int den, nom, index = 1;
  double *res = (double *) calloc(N, sizeof(double)); res[0] = 1.0;

  res[index++] = sin(PI/2.0); res[index++] = sin(PI/3.0);
  res[index++] = sin(PI/6.0);
  for (den = 12; den <= N; den <<= 1) {
    for (nom = 1; nom <= (den/2); nom += 2) {
        res[index++] = sin(((double)nom) * PI / (double)den);
    }
  }
  return res;
}

double *SSSFactors(int N) {
  int den, nom, index = 1;
  double *res = (double *) calloc(N, sizeof(double)); res[0] = 1.0;

  res[index++] = sin(PI/2.0);
  for (den = 5; den <= N; den <<= 1) {
    for (nom = 1; nom <= (den/2); nom += 2) {
        res[index++] = sin(((double)nom) * PI / (double)den);
    }
  }
  return res;
}
