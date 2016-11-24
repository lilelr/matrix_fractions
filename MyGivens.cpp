#include <math.h>
#include <stdio.h>

#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

using namespace std;
#define LL long long
const int maxn = 5;
int a_m, a_n;

void print_m(double a[maxn][maxn]) {
    cout << "测试" << endl;
    for (int i = 0; i < maxn; ++i) {
        for (int j = 0; j < maxn; ++j) {
            printf("%8.4lf\t", a[i][j]);
        }
        cout << endl;
    }
    cout << endl;
}

//计算向量a，长度为a_m, 的模
double norm(double a[maxn]) {
    double sum = 0;
    for (int i = 0; i < a_m; ++i) {
        sum += a[i] * a[i];
    }
    return sqrt(sum);
}

/**
 * 两个矩阵A和B相乘
 * @param A
 * @param A_r A的行数
 * @param A_c  A的列数
 * @param B
 * @param B_r  B的行数
 * @param B_c   C的列数
 * @param C   矩阵相乘后的积C
 */
void multiply(double A[maxn][maxn], int A_r, int A_c, double B[maxn][maxn], int B_r, int B_c, double C[maxn][maxn]) {
    for (int i = 0; i < A_r; ++i) {
        for (int j = 0; j < B_c; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < A_c; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

}

// m*n 的矩阵a 和一个向量u相乘,得到v向量
void multipy_vector(double a[maxn][maxn], int a_m, int a_n, double u[maxn], double v[maxn]) {
    for (int i = 0; i < a_m; ++i) {
        v[i] = 0;
        for (int j = 0; j < a_n; ++j) {
            v[i] += a[i][j] * u[j];
        }
    }
}

//向量a复制到向量b
void copy_vector(double a[maxn], int m, double b[maxn]) {
    for (int i = 0; i < m; ++i) {
        b[i] = a[i];
    }
}

// 把A矩阵复制到B矩阵
void copy(double ta[maxn][maxn], int m, int n, double tb[maxn][maxn]) {
//    cout<<"把A矩阵复制到B矩阵"<<endl;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            tb[i][j] = ta[i][j];
        }
    }
}


void givens(double A[maxn][maxn], int k, double Q[maxn][maxn]) {
    double u[maxn] = {0};
    double v[maxn] = {0};
    double c = 0;
    double s = 0;
    double temp[maxn][maxn] = {0};

    //记录一次旋转变换的T矩阵
    double T[maxn][maxn] = {0};
    for (int i = 0; i < a_m; ++i) {
        T[i][i] = 1;
    }

    // 求u向量
    for (int i = k; i < a_m; ++i) {
        u[i] = A[i][k];
    }

    // 旋转变换
    for (int i = k + 1; i < a_m; ++i) {

        c = u[k];
        s = u[i];
        double tnorm = sqrt(c * c + s * s);
        c = c / tnorm;
        s = s / tnorm;
        double H[maxn][maxn] = {0};
        for (int j = 0; j < a_m; ++j) {
            H[j][j] = 1;
        }
        H[k][k] = c;
        H[k][i] = s;
        H[i][k] = -s;
        H[i][i] = c;

        //计算u,v向量
        multipy_vector(H, a_m, a_m, u, v);
        copy_vector(v, a_m, u);

        //计算旋转矩阵的连乘，得到T矩阵
        multiply(H, a_m, a_m, T, a_m, a_m, temp);
        copy(temp, a_m, a_m, T);

    }

    //计算总的旋转矩阵Q
    multiply(T, a_m, a_m, Q, a_m, a_m, temp);
    copy(temp, a_m, a_m, Q);

    //计算新的A矩阵
    double temp2[maxn][maxn] = {0};
    multiply(T, a_m, a_m, A, a_m, a_n, temp2);
    copy(temp2, a_m, a_n, A);

}

int main() {
    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/in", "r", stdin);
    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/out", "w", stdout);

    cout << "请输入行数和列数，以空格分开：" << endl;
    cin >> a_m >> a_n;
    double A[maxn][maxn] = {0};
    cout << "请输入" << a_m << "行" << a_n << "列的A矩阵，元素以空格分开：" << endl;

    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_n; ++j) {
            cin >> A[i][j];
        }
    }

    // 打印A矩阵
    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_n; ++j) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    double Q[maxn][maxn] = {0};
    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_m; ++j) {
            Q[i][j] = 0;
            if (i == j) Q[i][j] = 1;
        }
    }

    int cyl = 0;
    if (a_m <= a_n) {
        cyl = a_m - 1;
    } else {
        cyl = a_n;
    }
    for (int k = 0; k < cyl; ++k) {
        givens(A, k, Q);
    }

    // 打印Q矩阵
    cout << "Q矩阵" << endl;
    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_m; ++j) {
            printf("%8.4lf\t", Q[j][i]);
        }
        cout << endl;
    }

    // 打印R矩阵
    cout << "R矩阵" << endl;

    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_n; ++j) {
            printf("%8.4lf\t", A[i][j]);
        }
        cout << endl;
    }


    return 0;
}
