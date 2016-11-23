
#include <math.h>
#include <stdio.h>

#define n 3 // maxN

void Matrix_print(double A[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%8.4lf\t", A[i][j]);
        printf("\n");
    }
}

double Matrix_norm(double a[n]) {
    double d = 0;
    for (int i = 0; i < n; i++)
        d += a[i] * a[i];
    return sqrt(d);
}

void Matrix_multiply(double A[n][n], double B[n][n], double C[n][n]) {

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int t = 0; t < n; t++)
                C[i][j] += A[i][t] * B[t][j];
        }
}

void Matrix_copy(double A[n][n], double B[n][n]) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = B[i][j];
}

void householder_trans(double A[n][n], int k, double Q[n][n]) {
    double a[n];
    for (int i = 0; i < n - k; i++)
        a[i] = 0;
    for (int i = n - k; i < n; i++)
        a[i] = A[i][n - k];
    a[n - k] -= Matrix_norm(a); // ?e1
    double d = Matrix_norm(a);
    for (int i = 0; i < n; i++)
        a[i] = a[i] / d;
    double H[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            H[i][j] = -2 * a[i] * a[j];

        H[i][i]++;
    } //��?H=I-2vvT
    if(k==2){
        int e=3;
    }
    double temp[n][n];
    Matrix_multiply(H, A, temp);
    Matrix_copy(A, temp);
    Matrix_multiply(H,Q, temp);
    Matrix_copy(Q, temp);
}

void Matrix_input(double A[n][n]) {
//    A[0][0] = 0;
//    A[0][1] = -20;
//    A[0][2] = -14;
//    A[1][0] = 3;
//    A[1][1] = 27;
//    A[1][2] = -4;
//    A[2][0] = 4;
//    A[2][1] = 11;
//    A[2][2] = -2;
    A[0][0] = 3;
    A[0][1] = 14;
    A[0][2] = 9;
    A[1][0] = 6;
    A[1][1] = 43;
    A[1][2] = 3;
    A[2][0] = 6;
    A[2][1] = 22;
    A[2][2] = 15;
}

int main4() {
    double Q[n][n];
    double A[n][n];

    Matrix_input(A);
    printf("A: \n");
    Matrix_print(A);
    int i;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Q[i][j] = 0;
        }
        Q[i][i] = 1;
    }

    for (i = n; i >= 2; i--) householder_trans(A, i, Q);
    printf("R: \n");
    Matrix_print(A);
    printf("Q: \n");
    Matrix_print(Q);
}
