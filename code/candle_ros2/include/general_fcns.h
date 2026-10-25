#ifndef GENERAL_FCNS_H
#define GENERAL_FCNS_H

#include <Eigen/Dense>
#include <iostream>
#include <algorithm>  // Per std::min e std::max
#include <cmath>      // Per std::acos
#include <limits>
using namespace std;
using namespace Eigen;

MatrixXd adjmap_fcn(const MatrixXd& mat);

MatrixXd coadjmap_fcn(const MatrixXd& mat);

MatrixXd ginv(const MatrixXd& g);

MatrixXd adjop_fcn(const VectorXd& v);

MatrixXd coadjop_fcn(const VectorXd& v);

MatrixXd coadjopbar_fcn(const VectorXd& v);

Matrix3d skewMat_fcn(const Vector3d& v);

Vector3d skewVec_fcn(const Matrix3d& sks_mat);

VectorXd VecFromHat_fcn(const MatrixXd& hat_mat);

MatrixXd HatFromVec_fcn(const VectorXd& v);

MatrixXd GaussQuadrature(int N_GQ);

MatrixXd loadMatrixFromTxt(const string& filename);


#endif // GENERAL_FCNS_H
