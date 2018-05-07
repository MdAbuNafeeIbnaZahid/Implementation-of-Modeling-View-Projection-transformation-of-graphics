#include <iostream>
#include <bits/stdc++.h>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using namespace std;

int main()
{

}


int cmain()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  cout << m << std::endl;

  MatrixXd mul = m * m;
  cout << mul << endl;

  getchar();
}
