#include "mymath.h"
#include <stdio.h>


double sign(int j) {
  return (j % 2 == 0) ? -1.0 : 1.0;
}

void printVector(double *x, int N, const char *name) {
  int j;
  for (j = 1; j <= N; j++) {
    printf("%s[%d] = %f\n", name, j, x[j]);
  }
}

int FastSN(double *x, double *y, double *w, double *S, int N, int skip);
int FastUN(double *x, double *y, double *w, double *S, int N, int skip);
int FastTN(double *x, double *y, double *w, double *S, int N, int skip);

int SlowSN(double *x, double *y, double *w, double *S, int N, int skip);
int SlowTN(double *x, double *y, double *w, double *S, int N, int skip);


int FastSN(double *x, double *y, double *w, double *S, int N, int skip) {
  if (N == 2) {

    x[skip] = y[skip]; // sin(pi/2) = 1

  } else if (N == 3) {

  } else if (N % 2 == 0) {

    // Vectors S and a are calculated.
    // s_j = y_j + y_{N-j}       a_j = y_j - y_{N-j}
    int j, subindex = 1;
    for (j = 1; j < N - 1; j++) {
      w[j*skip] = y[subindex*skip] + y[(N-subindex) * skip] * sign(j);
      if (sign(j) == -1.0) subindex++;
    }

    w[(N-1)*skip] = y[N/2];
  //  printVector(w, N-1, "w");

    int newSkip = skip << 1; int newN = N >> 1;

    FastSN(x, w, w, S, newN, newSkip);
    FastTN(x - skip, w - skip, w, S, newN, newSkip);
  }
  return 0;
}

int FastTN(double *x, double *y, double *w, double *S, int N, int skip) {
  if (N == 2) {
    double t2[3][3];
    t2[1][1] = 0.707106781186547524; t2[1][2] = 1; // sin(3*PI/4);
    t2[2][1] = 0.707106781186547524; t2[2][2] = -1.0; // sin(3*PI/2)

    int i, j;
    for (i = 1; i <= N; i++) {
      x[i*skip] = 0;
      for (j = 1; j <= N; j++) {
        x[i*skip] += t2[i][j] * y[j*skip];
      }
    }

  } else if (N == 3) {

  } else if (N % 2 == 0) {
    int newSkip = skip << 1; int newN = N >> 1;

    FastUN(x - skip, y - skip, w, S, newN, newSkip);
    FastTN(x, y, w, S, newN, newSkip);
  }
  return 0;
}


int FastUN(double *x, double *y, double *w, double *S, int N, int skip) {
  if (N == 2) {
    double u2[3][3];
    u2[1][1] = 0.3826834323650897717; u2[1][2] = 0.9238795325112867561;
    u2[2][1] = 0.9238795325112867561; u2[2][2] =-0.3826834323650897717;

    int i, j;
    for (i = 1; i <= N; i++) {
      x[i*skip] = 0;
      for (j = 1; j <= N; j++) {
        x[i*skip] += u2[i][j] * y[j*skip];
      }
    }

  } else if (N == 3) {

  } else if (N % 2 == 0) {
    int newSkip = skip << 1; int newN = N >> 1;

    // Calculations for vectors u
    int j;
    for (j = 1; j < (N/2); j++) {
      w[j*skip] = y[(N+1-2*j)*skip] - y[(N-2*j)*skip];
    }

    FastTN(x - skip, w - skip, w, S, newN, newSkip);

    for (j = 1; j < (N/2); j++) {
      w[j*skip] = y[(2*j)*skip] + y[(2*j+1)*skip];
    }

    //Calculations for vector v
    FastTN(x, w, w, S, newN, newSkip);
    //Now we have a and b intertwined in X


  }
  return 0;
}


int main() {
    int N    = 4;
    int skip = 1;
    double *x = malloc(N*skip*sizeof(double));
    double *y = malloc(N*skip*sizeof(double));
    int i;
    for (i = 0; i < N; i++) {
      y[i] = (double) i;
    }
    double *w = malloc(N * skip * sizeof(double));
    double *S = SFactors(N);

    FastSN(x, y, w, S, N, skip);

    printf("\nFinal result:\n");
    printVector(x, N-1, "x");
    SlowSN(x, y, w, S, N, skip);
    printf("\n");
    printVector(x, N-1, "x2");

}


int SlowTN(double *x, double *y, double *w, double *S, int N, int skip) {
  int i, j;
  for (i = 1; i <= N; i++) {
    x[i*skip] = 0.0;
    for (j = 1; j <= N; j++) {
      x[i*skip] += y[j*skip] *
                    sin((double) (2i-1) * (double) j * PI / 2*N);
    }
  }
  return 0;
}

int SlowUN(double *x, double *y, double *w, double *S, int N, int skip) {
  int i, j;
  for (i = 1; i <= N; i++) {
    x[i*skip] = 0.0;
    for (j = 1; j <= N; j++) {
      x[i*skip] += y[j*skip] *
                    sin((double) (2i-1) * (double) (2j - 1) * PI / 4*N);
    }
  }
  return 0;
}

int SlowSN(double *x, double *y, double *w, double *S, int N, int skip) {
  int i, j;
  for (i = 1; i < N; i++) {
    x[i*skip] = 0.0;
    for (j = 1; j < N; j++) {
      x[i*skip] += y[j*skip] * sin((double) i * (double) j * PI / N);
    }
  }
}
