//
//  main.cpp
//  Matrix_LU
//
//  Created by YuXiao on 10/4/16.
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
const int maxn = 1000;       //矩阵的最大维数
double nums[maxn][maxn];     //输入矩阵
double matrix_L[maxn][maxn]; //分解矩阵L
double matrix_U[maxn][maxn]; //分解矩阵U
int matrix_P[maxn][maxn];    //矩阵P
int row_interchange[maxn];
int m, n; //矩阵维数
bool ok;  //是否可分解

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

int LU() {

  cout << "请输入方阵的维数？" << endl;
  cin >> n;
  ok = true;
  if (n <= 0 || n > 1000) {
    cout << "输入的矩阵维数不合法！" << endl;
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
    return 0;
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

  return 0;
}

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

// 施密特正交化
int Schmidt() {
  //  cout << "请输入矩阵的维数(m行n列)？" << endl;
  freopen(
      "/Users/yuxiao/ClionProjects/Matrix_Fractorization/in",
      "r", stdin);
  freopen(
      "/Users/yuxiao/ClionProjects/Matrix_Fractorization/out",
      "w", stdout);
  //    cout << "fe" << endl;
  cin >> m >> n;
  ok = true;
  if (n <= 0 || n > 1000 || m <= 0 || m > 1000) {
    cout << "输入的矩阵维数不合法！" << endl;
  } else {
    cout << "请输入矩阵,元素以空格分开:" << endl;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        cin >> nums[i][j];
      }
    }

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
      //      int ta = Q_matrix[col].size();
      //      int tb = input_m[col].size();
      //      //      for (int i = 0; i < m; i++) {
      //      //        cout << Q_matrix[col][i] << " ";
      //      //      }
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
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        cout << Q_matrix[j][i] << " ";
      }
      cout << endl;
    }

    cout << "输出R矩阵:" << endl;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cout << R_matrix[i][j] << " ";
      }
      cout << endl;
    }
  }

  cout << endl;
  return 1;
}


vector<vector<double >> house_R; //m*n的矩阵
vector<vector<double>> house_Q; // m*m 的矩阵
vector<vector<double>> matrix_S;
/**
 *
 * @param cur 第cur轮迭代
 * @param a
 * @param b
 * @param m1
 * @param n1
 * @param m2
 * @param n2
 * @param ret
 */
void multiple(int cur,vector<vector<double >>& a,vector<vector<double >>& b,int m1,int n1,int m2,int n2, vector<vector<double>>& ret ){
    vector<vector<double >> ans ;
    int m = m1,n =n2;
  for (int i = 0; i < m; i++) {
    ans.push_back(vector<double>(n,0));
  }

  double item=0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k <a[i].size(); ++k) {
        item += a[i][k]*b[k][j];
      }
      ans[i][j] = item;
      // 返回下一次迭代的矩阵ret
      if(i!=0 && j!=0){
        ret[i-1].push_back(item);
      }
    }
  }
  // 得到R的第i行,即ans的第一行
  vector<double> tv(cur,0);
  for (int i = 0; i < n;  ++i) {
    tv.push_back(ans[0][i]);
  }

  // 得到S矩阵
//  for (int i = cur ; i < m; ++i) {
//    for (int j = 0; j <; ++j) {
//
//    }
//  }


}

int main1() {
  Schmidt();
  return 0;
}
