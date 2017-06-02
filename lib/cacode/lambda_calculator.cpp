// Copyright 2008 Michiaki Hamada
// Adapted from public domain code by Yi-Kuo Yu, NCBI

/**
 * See lambda_calculator.h
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

int Alphsize;

#include "nrutil.cpp"
#include "ludcmp.cpp"
#include "lubksb.cpp"
#include "lambda_calculator.h"

#define Epsilon 1.0e-36
#define E_bound 1.0e-12
#define Infty   1000000.0
#define min(A, B) ((A) > (B) ? (B) : (A) )
#define max(A, B) ((A) > (B) ? (A) : (B) )
#define bool int
#define true 1
#define false 0
double Lambda_UB; //lambda upper bound
double r_max_m, c_max_m; //min of each row's (column's) max

void makematrix(const double **mat_b, double **a, double lambda);

typedef struct Lambda {
    double min;
    double max;
    int flag;    // 1 means there is a range, -1 means no solution possible.
} Lambda;
typedef struct Sum {
    double value;
    int flag;   // 1 means no negative bg_freq, -1 means there is negative bg_freq
} Sum;

Lambda Find_JP(const double **mat_b, double la_min, double la_max, double **JP, double *p_in, double *q_in);

Sum Check_root(const double **mat_b, double **a, double lambda, double *p, double *q);

double Check_det(const double **mat_b, double **a, double lambda);

Sum Nail_lambda(const double **mat_b, int flag_sign, double lambda_min, double lambda_max, double *p, double *q,
                double *la_add);

double Nail_det(const double **mat_b, int flag_sign, double lambda_min, double lambda_max);

bool Check_range(const double **mat_b);

double *Locate_det_zero(const double **mat_b, int *); //pointer to root list are returned with how many of them by int


double calculate_lambda(const double **mat_b, int alpha_size,
                        double *p, double *q) {
    double **JP/*, *q, *p*/;
    int k;
    double *root_location;
    int N_root;
    Lambda Lambda_local;

    Alphsize = alpha_size;

    if (!Check_range(mat_b)) return -1.0;

    root_location = Locate_det_zero(mat_b, &N_root);
    if (root_location == NULL && N_root > 0) return -1.0;

    //q=dvector(1,Alphsize);
    //p=dvector(1,Alphsize);
    JP = dmatrix(1, Alphsize, 1, Alphsize);

    if (N_root == 0) {
        Lambda_local = Find_JP(mat_b, 0, Lambda_UB, JP, p, q);
        if (1 == Lambda_local.flag) { // sensible solution found
            // Remember to find the right place to free the vectors
            //free_dvector(p, 1,Alphsize);
            //free_dvector(q, 1,Alphsize);
            free(root_location);
            free_dmatrix(JP, 1, Alphsize, 1, Alphsize);
            return (Lambda_local.min + Lambda_local.max) / 2.0;
        } else if (-1 == Lambda_local.flag) {
            //printf("matrix pass first screening but no sensible solution found. :-( \n");
        }
    } else if (N_root > 0) {
        //printf("N_root = %d for this matirx \n", N_root);
        //for (i=0;i<N_root;i++) printf("N_root[%d] = %lf\n",i,root_location[i]);
        for (k = 0; k <= N_root; k++) {
            if (k == 0) {
                Lambda_local.min = 0;
                Lambda_local.max = root_location[0];
            }
            else if (k == N_root) {
                Lambda_local.min = root_location[N_root - 1];
                Lambda_local.max = Lambda_UB + Epsilon;
            }
            else {
                Lambda_local.min = root_location[k - 1];
                Lambda_local.max = root_location[k];
            }
            Lambda_local = Find_JP(mat_b, Lambda_local.min, Lambda_local.max, JP, p, q);
            if (1 == Lambda_local.flag) { // sensible solution found
                //free_dvector(p, 1,Alphsize);
                //free_dvector(q, 1,Alphsize);
                free(root_location);
                free_dmatrix(JP, 1, Alphsize, 1, Alphsize);
                return (Lambda_local.min + Lambda_local.max) / 2.0;
            } else if (-1 == Lambda_local.flag) {
                //printf("matrix pass first screening but still no sensible solution found. :-( \n");
            }
        }
    }
    // Remember to find the right place to free the vectors
    //free_dvector(p, 1,Alphsize);
    //free_dvector(q, 1,Alphsize);
    free(root_location);
    free_dmatrix(JP, 1, Alphsize, 1, Alphsize);
    return -1.0;
}

bool Check_range(const double **mat_b) {
    int pos_flag_r, neg_flag_r;
    int pos_flag_c, neg_flag_c;
    double r_max, c_max; //max of each row (or column)

    int L_r = 0, L_c = 0; // number of zero-score rows (columns)

    // First make sure each row and column have both pos and neg entries
    r_max_m = c_max_m = 100000000000.0;
    for (int i = 1; i <= Alphsize; i++) {
        r_max = 0;
        c_max = 0;
        pos_flag_r = -1;
        neg_flag_r = -1;
        pos_flag_c = -1;
        neg_flag_c = -1;
        for (int j = 1; j <= Alphsize; j++) {
            if (mat_b[i][j] > 0) {
                if (mat_b[i][j] > r_max) r_max = mat_b[i][j];
                pos_flag_r = 1;
            } else if (mat_b[i][j] < 0) neg_flag_r = 1;
            if (mat_b[j][i] > 0) {
                if (mat_b[j][i] > c_max) c_max = mat_b[j][i];
                pos_flag_c = 1;
            } else if (mat_b[j][i] < 0) neg_flag_c = 1;
        }
        if ((pos_flag_r == -1) || (neg_flag_r == -1) || (pos_flag_c == -1) || (neg_flag_c == -1)) {
            if ((pos_flag_r == -1) && (neg_flag_r == -1)) {
                printf("only zero score at row  %d\n", i);
                L_r++;
            } else if ((pos_flag_c == -1) && (neg_flag_c == -1)) {
                printf("only zero score at column %d\n", i);
                L_c++;
            } else {
                //printf("all positive or all negative at row or column %d\n", i);
                //printf("therefore invalid matrix. exit now. \n");
                return false;
                //exit(1);
            }
        }
        if ((r_max < r_max_m) && (r_max > 0)) r_max_m = r_max;
        if ((c_max < c_max_m) && (c_max > 0)) c_max_m = c_max;
    }


    // Find the upper bound for lambda
    if (r_max_m > c_max_m) {
        Lambda_UB = 1.1 * log(1.0 * Alphsize - L_r) / r_max_m;
    } else {
        Lambda_UB = 1.1 * log(1.0 * Alphsize - L_c) / c_max_m;
    }
    //printf("the upper bound for lambda is %lf\n", Lambda_UB);
    return true;
}


double Check_det(const double **mat_b, double **a, double lambda) {
    double d;
    int i, /*j,*/ *indx;

    indx = ivector(1, Alphsize);
    makematrix(mat_b, a, lambda);
    ludcmp(a, Alphsize, indx, &d);
    for (i = 1; i <= Alphsize; i++) d *= a[i][i];
    free_ivector(indx, 1, Alphsize);
    return d;  //returning the determinant
}


Sum Check_root(const double **mat_b, double **a, double lambda, double *p, double *q) {
    double **y, /* *col,*/ d;
    //double sum = 0.0;
    int i, j;//, *indx;
    Sum Sum_here;

    y = dmatrix(1, Alphsize, 1, Alphsize);
    //indx = ivector(1,Alphsize);
    int indx[Alphsize + 1];
    //col = dvector(1,Alphsize);
    double col[Alphsize + 1];

    makematrix(mat_b, a, lambda);
    ludcmp(a, Alphsize, indx, &d);
    Sum_here.value = 0.0;
    for (i = 1; i <= Alphsize; i++) q[i] = 0.0;
    for (j = 1; j <= Alphsize; j++) {
        for (i = 1; i <= Alphsize; i++){
            col[i] = 0.0;
        }
        col[j] = 1.0;
        lubksb(a, Alphsize, indx, col);
        p[j] = 0.0;
        for (i = 1; i <= Alphsize; i++) {
            y[i][j] = col[i];
            Sum_here.value += y[i][j];
            p[j] += y[i][j];
            q[i] += y[i][j];
        }
    }

    Sum_here.flag = 1;
    for (i = 1; i < Alphsize; i++) {
        if ((p[i] < 0) || (q[i] < 0)) {
            Sum_here.flag = -1;
            //printf("problematic freq. p[%d] = %.4f q[%d]=%.4f\n",i,p[i],i,q[i]);
        }
    }
    free_dmatrix(y, 1, Alphsize, 1, Alphsize);
    return Sum_here;
}


double *Locate_det_zero(const double **mat_b, int *N_root_add) {
    double **a/*,  *q, *p */; // a is the exponentiated matrix of socres, p and q are bg_freqs
    int i/*,j,k*/;
    int N;  // number of points for first round
    int flag_sign;
    double lambda/*, l_tmp, sum,  sum_min, sum_max */;
    double lambda_root, dlambda /*, dsum=0.5 */;
    //double *l_here, *s_here;
    double root[5000];
    double *root_temp;
    //double error=0.000000000001;
    int zero_monitor = 0;  // record number of zeros found in the range
    //int flag;

    a = dmatrix(1, Alphsize, 1, Alphsize);
    //Sum_local = (Sum *)malloc(sizeof(Sum));
    //Lambda_local = (Lambda *)malloc(sizeof(Lambda));

    N = 2 + max(400, ((int) (Lambda_UB - 0) / 0.005));
    //printf("N = %d in Locate_det_zero\n", N);
    dlambda = (Lambda_UB) / (N * 1.0);
    //l_here = (double *)malloc((N+1)*sizeof(double));
    //s_here = (double *)malloc((N+1)*sizeof(double));
    double l_here[N + 1];
    double s_here[N + 1];

    for (i = 0; i < N; i++) {
        lambda = (i + 1) * dlambda;
        s_here[i] = Check_det(mat_b, a, lambda);
        l_here[i] = lambda;
    }

    if (s_here[0] < 0.0) flag_sign = -1;
    if (s_here[0] > 0.0) flag_sign = 1;
    if (fabs(s_here[0]) / exp(l_here[0] * (r_max_m + c_max_m) / 2.0) <= Epsilon) {
        root[zero_monitor++] = l_here[0];
        flag_sign = 0;
    }

    for (i = 1; i < N; i++) {
        if ((flag_sign != 0) && (fabs(s_here[i]) > Epsilon)) {
            if (s_here[i - 1] * s_here[i] < 0) {
                //printf("occurring at regular places\n");
                lambda_root = Nail_det(mat_b, flag_sign, l_here[i - 1], l_here[i]);
                root[zero_monitor++] = lambda_root;
                flag_sign = -flag_sign;  // the flag switch sign after one sol found
                //printf("a (regular) root of det found at %12.10f, i= %d\n", lambda_root,i);
            }
        } else {
            if (s_here[i] < 0.0) flag_sign = -1;
            if (s_here[i] > 0.0) flag_sign = 1;
            if (fabs(s_here[i]) / exp(l_here[i] * (r_max_m + c_max_m) / 2.0) <= Epsilon) {
                root[zero_monitor++] = l_here[i];
            }
        }
    }
    //printf("total number of solution found in range is %d\n", i_monitor);
    root_temp = (double *) malloc(zero_monitor * sizeof(double));
    *N_root_add = zero_monitor;
    if (zero_monitor > 0) {
        if (zero_monitor >= N / 4) {
            //printf("It is likely that uniform zero determinant is occurring.\n");
            //printf("number of small det points = %d out of %d, exit now....\n",zero_monitor, N);
            free(root_temp);
            return NULL;
            //exit(1);
        }
        for (i = 0; i < zero_monitor; i++) {
            root_temp[i] = root[i];
            //printf("root_location[%d] = %lf\n",i,root_temp[i]);
        }
    }
    free_dmatrix(a, 1, Alphsize, 1, Alphsize);
    return root_temp;

}


Lambda Find_JP(const double **mat_b, double la_min, double la_max, double **JP, double *p_in, double *q_in) {
    double **a, *q, *p; // a is the exponentiated matrix of socres, p and q are bg_freqs
    int i, j/*,k*/;
    int N;  // number of points for first round
    double lambda/*, l_tmp, sum, sum_min, sum_max*/;
    double lambda_max, lambda_min, lambda_final, dlambda/*, dsum=0.5*/;
    //double *l_here, *s_here;
    //double error=0.000000000000001;
    //int validity_flag; // 1 means valid, -1 means not valid.
    int flag_sign;       // 1 means small lambda sum > 1, -1 means otherwise
    int flag_done = -1;       // 1 means find sensible solution, -1 means sensible not found
    int i_monitor = 0;          // record number of solution found in the range, including nonsense ones
    int j_monitor;

    Lambda Lambda_local;
    //Sum *Sum_local;
    Sum Sum_local;

    lambda_min = la_min;
    lambda_max = la_max;
    q = q_in;
    p = p_in;
    a = dmatrix(1, Alphsize, 1, Alphsize);
    //Sum_local = (Sum *)malloc(sizeof(Sum));
    //Lambda_local = (Lambda *)malloc(sizeof(Lambda));

    N = 2 + max(400, ((int) (lambda_max - lambda_min) / 0.005));
    //printf("N = %d in Find_JP\n", N);
    dlambda = (lambda_max - lambda_min) / (N * 1.0);
    //l_here = (double *)malloc((N+1)*sizeof(double));
    //s_here = (double *)malloc((N+1)*sizeof(double));
    double l_here[N + 1];
    double s_here[N + 1];
    //printf("lambda_min enter = %12.10e, lambda_max = %12.10f\n", lambda_min, lambda_max);
    for (i = 0; i < N - 1; i++) {
        lambda = lambda_min + (i + 1) * dlambda;
        makematrix(mat_b, a, lambda);
        Sum_local = Check_root(mat_b, a, lambda, p, q);
        l_here[i] = lambda;
        s_here[i] = Sum_local.value - 1.0;
        //printf("scan %d th time in Find_JP\n",i );
    }
    //printf("finish first time scanining in Find_JP\n");
    if (s_here[0] < 0.0) flag_sign = -1;
    else if (s_here[0] > 0.0) flag_sign = 1;
    else if (s_here[0] == 0.0) {  //needs refined definition on flag_sign
        printf("enter the exact hit, rarely occurs other than when lambda = 0 \n");
        j_monitor = 1;
        flag_sign = 0;
        while ((flag_sign == 0) && (j_monitor < N)) {
            Sum_local = Check_root(mat_b, a, l_here[0] + j_monitor * dlambda / N, p, q);
            if (Sum_local.value > 1.0) {
                flag_sign = 1;
            } else if (Sum_local.value < 1.0) {
                flag_sign = -1;
            }
            j_monitor++;
        }
    }

    for (i = 1; i < N; i++) {  // should be N-1 ???
        if (flag_sign == 0) {
            printf("flag_sign = 0 \n");
            exit(1);
        }
        if (s_here[i - 1] * s_here[i] < 0) {
            lambda_min = l_here[i - 1];
            lambda_max = l_here[i];
            Sum_local = Nail_lambda(mat_b, flag_sign, lambda_min, lambda_max, p, q, &lambda_final);
            if (Sum_local.flag == 1) {
                i = N;
                flag_done = 1;
                Lambda_local.flag = 1;
                Lambda_local.min = lambda_final, Lambda_local.max = lambda_final;
            }
            flag_sign = -flag_sign;  // the flag switch sign after one sol found
            i_monitor++;
        }
    }

    if (flag_done == 1) {
        // Write correct JP to the matrix
        makematrix(mat_b, a, lambda_final);
        for (i = 1; i <= Alphsize; i++) {
            for (j = 1; j <= Alphsize; j++) {
                JP[i][j] = a[i][j] * p[i] * q[j];
            }
        }
        free_dmatrix(a, 1, Alphsize, 1, Alphsize);
        return Lambda_local;
    } else if (flag_done == -1) {
        //printf("no sensible solution in the plausible x range: (%lf,%lf)\n", la_min, la_max);
        Lambda_local.flag = -1;
        Lambda_local.min = 0;
        Lambda_local.max = Infty;
        return Lambda_local;
    }
    // never come here
    return Lambda_local;
}


Sum Nail_lambda(const double **mat_b, int flag_sign, double lambda_min, double lambda_max, double *p, double *q,
                double *lam_add) {
    double **a;
    double lambda;

    //Sum *Sum_local;
    Sum Sum_local;
    a = dmatrix(1, Alphsize, 1, Alphsize);
    //Sum_local = (Sum *)malloc(sizeof(Sum));

    lambda = (lambda_min + lambda_max) / 2.0;
    Sum_local = Check_root(mat_b, a, lambda, p, q);
    while (fabs(Sum_local.value - 1.0) > E_bound) {
        if (flag_sign * (Sum_local.value - 1.0) < 0) lambda_max = lambda;
        else if (flag_sign * (Sum_local.value - 1.0) > 0) lambda_min = lambda;

        // Added by MCF to avoid infinite loop:
        if (lambda == (lambda_min + lambda_max) / 2.0) {
            Sum_local.flag = -1;
            break;
        }

        lambda = (lambda_min + lambda_max) / 2.0;
        Sum_local = Check_root(mat_b, a, lambda, p, q);
    }
    free_dmatrix(a, 1, Alphsize, 1, Alphsize);
    *lam_add = lambda;
    return Sum_local;
}


double Nail_det(const double **mat_b, int flag_sign, double lambda_min, double lambda_max) {
    double **a;
    double lambda;
    double value;

    a = dmatrix(1, Alphsize, 1, Alphsize);

    lambda = (lambda_min + lambda_max) / 2.0;
    value = Check_det(mat_b, a, lambda);
    while ((fabs(value) > E_bound) && (lambda > 0)) {
        if (flag_sign * (value) < 0) lambda_max = lambda;
        else if (flag_sign * (value) > 0) lambda_min = lambda;
        lambda = (lambda_min + lambda_max) / 2.0;
        value = Check_det(mat_b, a, lambda);
    }
    free_dmatrix(a, 1, Alphsize, 1, Alphsize);
    return lambda;
}

void makematrix(const double **mat_b, double **a, double lambda) {

    int i, j;

    for (i = 1; i <= Alphsize; i++)
        for (j = 1; j <= Alphsize; j++) {
            *(*(a + i) + j) = exp(lambda * mat_b[i][j]);
        }
}



