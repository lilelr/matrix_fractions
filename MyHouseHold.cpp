//#include <math.h>
//#include <stdio.h>
//
//#include <iostream>
//#include <map>
//#include <queue>
//#include <set>
//#include <string>
//#include <vector>
//
//using namespace std;
//#define LL long long
//const int maxn = 5;
//int a_m, a_n;
//
//void print_m(double a[maxn][maxn]){
//    cout<<"测试"<<endl;
//    for (int i = 0; i < maxn; ++i) {
//        for (int j = 0; j < maxn; ++j) {
//            printf("%8.4lf\t", a[i][j]);
//        }
//        cout<<endl;
//    }
//}
//
////计算向量a，长度为a_m, 的模
//double norm(double a[maxn]) {
//    double sum = 0;
//    for (int i = 0; i < a_m; ++i) {
//        sum += a[i] * a[i];
//    }
//    return sqrt(sum);
//}
//
///**
// * 两个矩阵A和B相乘
// * @param A
// * @param A_r A的行数
// * @param A_c  A的列数
// * @param B
// * @param B_r  B的行数
// * @param B_c   C的列数
// * @param C   矩阵相乘后的积C
// */
//void multiply(double A[maxn][maxn], int A_r, int A_c, double B[maxn][maxn], int B_r, int B_c, double C[maxn][maxn]) {
//    for (int i = 0; i < A_r; ++i) {
//        for (int j = 0; j < B_c; ++j) {
//            C[i][j] = 0;
//            for (int k = 0; k < A_c; ++k) {
//                C[i][j] += A[i][k] * B[k][j];
//            }
//        }
//    }
////    print_m(C);
//
//}
//
//// 把A矩阵复制到B矩阵
//void copy(double ta[maxn][maxn], int m, int n, double tb[maxn][maxn]) {
////    cout<<"把A矩阵复制到B矩阵"<<endl;
//    for (int i = 0; i < m; ++i) {
//        for (int j = 0; j < n; ++j) {
//            tb[i][j] = ta[i][j];
//        }
//    }
//}
//
//
//void householder(double A[maxn][maxn], int k, double Q[maxn][maxn]) {
//    double u[maxn] = {0};
//    for (int i = 0; i < k; ++i) {
//        u[i] = 0;
//    }
//    for (int i = k; i < a_m; ++i) {
//        u[i] = A[i][k];
//    }
//    u[k] -= norm(u);
//    double u_norm = norm(u);
//
//    // 求u向量
//    for (int i = 0; i < a_m; ++i) {
//        u[i] /= u_norm;
//    }
//
//    double H[maxn][maxn] = {0};
//    for (int i = 0; i < a_m; ++i) {
//        for (int j = 0; j < a_m; ++j) {
//            H[i][j] = -2 * u[i] * u[j];
//        }
//        H[i][i] += 1;
//    }
//
//    if(k==1){
////        print_m(A);
//        int ew=1;
//    }
//
//    double temp[maxn][maxn] = {0};
//    multiply(H, a_m, a_m,  A, a_m, a_n,  temp);
////    cout<<"测试A1"<<endl;
////    print_m(temp);
//
//    copy(temp, a_m, a_n,  A);
//    double temp2[maxn][maxn] = {0};
//    multiply( H, a_m, a_m,  Q, a_m, a_m,  temp2);
//    copy(temp2,a_m,a_m,Q);
//}
//
//int main() {
//    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/in", "r", stdin);
//    freopen("/Users/yuxiao/ClionProjects/Matrix_Fractorization/out", "w", stdout);
//
//    cout << "请输入行数和列数，以空格分开：" << endl;
//    cin >> a_m >> a_n;
//    double A[maxn][maxn] = {0};
//    cout << "请输入"<<a_m <<"行"<<a_n<<"列的A矩阵，元素以空格分开："<< endl;
//
//    for (int i = 0; i < a_m; ++i) {
//        for (int j = 0; j < a_n; ++j) {
//            cin >> A[i][j];
//        }
//    }
//
//    // 打印A矩阵
//    for (int i = 0; i < a_m; ++i) {
//        for (int j = 0; j < a_n; ++j) {
//            cout<<A[i][j]<<" ";
//        }
//        cout<<endl;
//    }
//
//    double Q[maxn][maxn]={0};
//    for (int i = 0; i < a_m; ++i) {
//        for (int j = 0; j < a_m; ++j) {
//            Q[i][j] = 0;
//            if (i == j) Q[i][j] = 1;
//        }
//    }
//
//    int cyl = 0;
//    if(a_m<=a_n){
//        cyl=a_m-1;
//    }else{
//        cyl = a_n;
//    }
//    for (int k = 0; k < cyl; ++k) {
//        householder(A,k,Q);
//    }
//
//    // 打印Q矩阵
//    cout<<"Q矩阵"<<endl;
//    for (int i = 0; i < a_m; ++i) {
//        for (int j = 0; j < a_m; ++j) {
//            printf("%8.4lf\t", Q[j][i]);
//        }
//        cout<<endl;
//    }
//
//    // 打印R矩阵
//    cout<<"R矩阵"<<endl;
//
//    for (int i = 0; i < a_m; ++i) {
//        for (int j = 0; j < a_n; ++j) {
//            printf("%8.4lf\t", A[i][j]);
//        }
//        cout<<endl;
//    }
//
//
//    return 0;
//}
