//
//  main.cpp
//  Matrix_Fracorization
//
//  Created by YuXiao
//  Copyright © 2016 YuXiao. All rights reserved.
//
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

using namespace std;
const int maxn = 10;       //矩阵的最大维数
double nums[maxn][maxn];     //输入矩阵
double matrix_L[maxn][maxn]; //分解矩阵L
double matrix_U[maxn][maxn]; //分解矩阵U
int matrix_P[maxn][maxn];    //矩阵P
int row_interchange[maxn];
int m, n; //矩阵维数
bool ok;  //是否可分解

// LU 分解代码开始————————————————
//交换两个数
void swapA_B(double &a, double &b) {
    double t = a;
    a = b;
    b = t;
}

//寻找该列绝对值最大的元素进行交换
void interchange(int i) {
    int k = i + 1;
    int ans = i;
    double max_pivot = abs(nums[i][i]);
    for (; k < n; k++) {
        if (max_pivot < abs(nums[k][i])) {
            max_pivot = abs(nums[k][i]);
            ans = k;
        }
    }
    if (ans != i) {
        for (int j = 0; j < n; j++) {
            swapA_B(nums[i][j], nums[ans][j]);
        }
    }

    int temp = row_interchange[i];
    row_interchange[i] = row_interchange[ans];
    row_interchange[ans] = temp;
    return;
}

// LU 分解程序入口
// haha
void LU() {
    cout << "矩阵LU 分解开始，所输入方阵维数不能超过" << maxn << "维--------------------" << endl;
    cout << "请输入方阵的维数？" << endl;
    cin >> n;
    ok = true;
    if (n <= 0 || n > maxn) {
        cout << "输入的矩阵维数不合法！" << endl;
        return;
    } else {
        cout << "请输入方阵,元素以空格分开:" << endl;
        for (int i = 0; i < n; i++) {
            row_interchange[i] = i;
            for (int j = 0; j < n; j++) {
                cin >> nums[i][j];
            }
        }
    }

    cout << endl;

    cout << "matrix_A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << nums[i][j] << " ";
        }
        cout << endl;
    }
    //  int rank = gauss_jordan(nums, n, n);
    //  cout << "矩阵的秩" << rank << endl;
    memset(matrix_L, 0, sizeof(matrix_L));
    memset(matrix_U, 0, sizeof(matrix_U));
    memset(matrix_P, 0, sizeof(matrix_P));

    //矩阵LU分解
    for (int i = 0; i < n; i++) {
        interchange(i);
        matrix_P[i][row_interchange[i]] = 1;

        double pivot = nums[i][i];
        if (pivot == 0) {
            ok = false;
            break;
        } else {
            for (int k = i + 1; k < n; k++) {
                int j = i;
                double times = nums[k][j] / pivot;
                //        matrix_L[k][j] = times;
                nums[k][j] = times;
                j++;
                for (; j < n; j++) {
                    nums[k][j] -= times * nums[i][j];
                }
            }
        }
    }

    if (!ok) {
        cout << "输入矩阵不能进行LU分解" << endl;
    }

    //矩阵 L，U求解
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            matrix_L[i][j] = nums[i][j];
        }
        matrix_L[i][i] = 1;
        for (int j = i; j < n; j++) {
            matrix_U[i][j] = nums[i][j];
        }
    }

    cout << endl;
    cout << "输入矩阵A的LU分解如下:" << endl;

    cout << "矩阵P:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix_P[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl;
    cout << "矩阵L:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix_L[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "矩阵U:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix_U[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
// LU 分解代码结束————————————————

// 施密特正交分解代码开始————————————————

//计算向量a，长度为m, 的模
double norm(vector<double> a, int m) {
    double sum = 0;
    for (int i = 0; i < m; i++) {
        sum += a[i] * a[i];
    }
    sum = sqrt(sum);
    return sum;
}

// 计算两个向量组的点积累
double inner_product(vector<double> a, vector<double> b, int n) {
    double product = 0;
    for (int i = 0; i < n; i++) {
        product += a[i] * b[i];
    }
    return product;
}

// 施密特正交分解程序入口
void Schmidt() {
    cout << "-----------------------施密特正交分解，所输入矩阵维数不能超过" << maxn << "维-------------------" << endl;
    cout << "请输入待处理矩阵A的维数(m行n列,以空格分开)？" << endl;

    cin >> m >> n;
    ok = true;
    if (n <= 0 || n > maxn || m <= 0 || m > maxn) {
        cout << "输入的矩阵维数不合法！" << endl;
    } else {
        cout << "请输入矩阵A,元素以空格分开:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cin >> nums[i][j];
            }
        }

        cout<<"\n施密特正交分解结果如下："<<"--------------------"<<endl;
        cout << "输入矩阵matrix_A:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout << nums[i][j] << " ";
            }
            cout << endl;
        }

        //记录列向量组
        vector<vector<double>> input_m;
        for (int j = 0; j < n; j++) {
            vector<double> t_c;
            for (int i = 0; i < m; i++) {
                t_c.push_back(nums[i][j]);
            }
            input_m.push_back(t_c);
        }
        //    cout << "input_m.size" << input_m.size() << endl;

        // 结果Q，R矩阵
        vector<vector<double>> Q_matrix;
        vector<vector<double>> R_matrix;
        for (int i = 0; i < n; i++) {
            R_matrix.push_back(vector<double>(n, 0));
            Q_matrix.push_back(vector<double>(m));
        }

        for (int col = 0; col < n; col++) {
            for (int i = 0; i < col; i++) {
                R_matrix[i][col] = inner_product(Q_matrix[i], input_m[col], m);
            }

            Q_matrix[col] = input_m[col];

            cout << endl;
            for (int i = 0; i < col; i++) {
                for (int x = 0; x < m; x++) {
                    Q_matrix[col][x] -= R_matrix[i][col] * Q_matrix[i][x];
                }
            }

            double temp_norm = norm(Q_matrix[col], m);
            R_matrix[col][col] = temp_norm;
            for (int i = 0; i < m; i++) {
                Q_matrix[col][i] /= temp_norm;
            }
        }

        //输出Q矩阵
        cout << "输出Q矩阵:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.3lf\t", Q_matrix[j][i]);
            }
            cout << endl;
        }


        cout << "输出R矩阵:" << endl;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.3lf\t", R_matrix[i][j]);
            }
            cout << endl;
        }
    }

    cout << endl;
}
// 施密特正交分解代码结束————————————————

// HouseHold 分解代码开始---------------
int a_m, a_n;// 记录输入矩阵的维数

void print_m(double a[maxn][maxn]) {
    cout << "测试" << endl;
    for (int i = 0; i < maxn; ++i) {
        for (int j = 0; j < maxn; ++j) {
            printf("%.3lf\t", a[i][j]);
        }
        cout << endl;
    }
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
//    print_m(C);

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

// houseHolder 分解迭代过程
void householder(double A[maxn][maxn], int k, double Q[maxn][maxn]) {
    double u[maxn] = {0};
    for (int i = 0; i < k; ++i) {
        u[i] = 0;
    }
    for (int i = k; i < a_m; ++i) {
        u[i] = A[i][k];
    }
    u[k] -= norm(u);
    double u_norm = norm(u);

    // 求u向量
    for (int i = 0; i < a_m; ++i) {
        u[i] /= u_norm;
    }

    double H[maxn][maxn] = {0};
    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_m; ++j) {
            H[i][j] = -2 * u[i] * u[j];
        }
        H[i][i] += 1;
    }

    double temp[maxn][maxn] = {0};
    multiply(H, a_m, a_m, A, a_m, a_n, temp);

    //求每一步迭代后的新的R矩阵
    copy(temp, a_m, a_n, A);

    //求每一步迭代后的新的Q矩阵
    double temp2[maxn][maxn] = {0};
    multiply(H, a_m, a_m, Q, a_m, a_m, temp2);
    copy(temp2, a_m, a_m, Q);
}

// HouseHold 分解程序入口
void houseHold_main() {
//    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/in", "r", stdin);
//    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/out", "w", stdout);
    cout << "// HouseHold 分解开始，所输入矩阵维数不能超过" << maxn << "维--------------" << endl;
    cout << "请输入要处理矩阵A的行数m和列数n，以空格分开：" << endl;
    cin >> a_m >> a_n;
    if (a_m <= 0 || a_m > maxn || a_n <= 0 || a_n > maxn) {
        cout << "输入的矩阵维数不合法！" << endl;
        return;
    }
    double A[maxn][maxn] = {0};
    cout << "请输入" << a_m << "行" << a_n << "列的矩阵A，元素以空格分开：" << endl;

    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_n; ++j) {
            cin >> A[i][j];
        }
    }

    // 打印A矩阵
    cout << "\nHouseHold分解结果如下：" << endl;

    cout << "矩阵A:" << endl;

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
    //houselder 迭代开始
    for (int k = 0; k < cyl; ++k) {
        householder(A, k, Q);
    }

    // 打印Q矩阵
    cout << "Q矩阵" << endl;
    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_m; ++j) {
            printf("%.3lf\t", Q[j][i]);
        }
        cout << endl;
    }

    // 打印R矩阵
    cout << "R矩阵" << endl;

    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_n; ++j) {
            printf("%.3lf\t", A[i][j]);
        }
        cout << endl;
    }
}
// houseHold 分解代码结束--------

// Givens 分解代码开始-------------
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

// givens 分解迭代过程
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

    //计算新的R矩阵
    double temp2[maxn][maxn] = {0};
    multiply(T, a_m, a_m, A, a_m, a_n, temp2);
    copy(temp2, a_m, a_n, A);
}

// givens 分解程序入口
void givens_main() {
//    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/in", "r", stdin);
//    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/out", "w", stdout);
    cout << "Givens分解开始,所输入矩阵维数不能超过" << maxn << "维------------------------------" << endl;
    cout << "请输入要处理矩阵A的行数m和列数n，以空格分开：" << endl;
    cin >> a_m >> a_n;
    if (a_m <= 0 || a_m > maxn || a_n <= 0 || a_n > maxn) {
        cout << "输入的矩阵维数不合法！" << endl;
        return;
    }

    double A[maxn][maxn] = {0};
    cout << "请输入" << a_m << "行" << a_n << "列的A矩阵，元素以空格分开：" << endl;

    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_n; ++j) {
            cin >> A[i][j];
        }
    }

    // 打印A矩阵
    cout << "\nGivens分解结果如下：" << endl;
    cout << "矩阵A:" << endl;
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
            printf("%.3lf\t", Q[j][i]);
        }
        cout << endl;
    }

    // 打印R矩阵
    cout << "R矩阵" << endl;

    for (int i = 0; i < a_m; ++i) {
        for (int j = 0; j < a_n; ++j) {
            printf("%.3lf\t", A[i][j]);
        }
        cout << endl;
    }
}
// givens 分解代码结束-----------

int main() {
    cout<<"欢迎使用矩阵分解程序，输入1：LU分解，2：施密特正交分解，3：HouseHold分解，4：Givens分解，0：结束"<<endl;

    int opt;
    while (cin>>opt){
        switch (opt){
            case 1:
                LU();
                break;
            case 2:
                Schmidt();
                break;
            case 3:
                houseHold_main();
                break;
            case 4:
                givens_main();
                break;
            case 0:
                cout<<"谢谢使用！"<<endl;
                return 0;
        }
        cout<<"\n请继续选择分解算法，输入1：LU分解，2：施密特正交分解，3：HouseHold分解，4：Givens分解，0：结束"<<endl;
    }
    return 0;
}
