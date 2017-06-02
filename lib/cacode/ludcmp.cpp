// Public domain code from Yi-Kuo Yu & Stephen Altschul, NCBI

#include <math.h>

#define TINY 1.0e-20;

double *dvector(int, int);

void nrerror(const char *);

void free_dvector(double *, int, int);

void ludcmp(double **a, int n, int *indx, double *d) {
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;

    vv = dvector(1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(a[i][j])) > big) big = temp;
        if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        if (j != n) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv, 1, n);
}

#undef TINY
