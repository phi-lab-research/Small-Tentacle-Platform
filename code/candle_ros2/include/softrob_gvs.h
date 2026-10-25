#ifndef SOFTROB_GVS_H
#define SOFTROB_GVS_H

#include "general_fcns.h"

MatrixXd compute_Bq1(const double& X,const double& L);
MatrixXd compute_Bq2(const double& X,const double& L);
MatrixXd compute_Bq3(const double& X,const double& L);
MatrixXd compute_Bq4(const double& X,const double& L);

VectorXd compute_xi(const VectorXd& q, const double& X, const int& sel, const VectorXi& n, const double& L);
VectorXd compute_xidot(const VectorXd& qdot, const double& X, const int& sel, const VectorXi& n, const double& L);

MatrixXd compute_J(const MatrixXd& J0, const VectorXd &q, const VectorXd &h, int h_pt, double L, int sel, const VectorXi& n);
MatrixXd compute_T(double theta, const MatrixXd &adjop_omega);
MatrixXd compute_Z(double h, const VectorXd &xi1, const VectorXd &xi2, int sel, double X1, double X2, double L);
MatrixXd direct_kin(const MatrixXd &g0, const VectorXd &q, const VectorXd &h, int h_pt, double L,int sel, const VectorXi &n);
MatrixXd Mag_fcn(double h, const VectorXd &xi1, const VectorXd &xi2);

MatrixXd compute_intGQ(const int& n, const MatrixXd& IntArg, const VectorXd& X, const VectorXd& W, const double& L);
                            
MatrixXd compute_AT(const MatrixXd& Btau_i, const MatrixXd& Bq, const VectorXd& X, const VectorXd& W, const double& L);
VectorXd compute_F(const MatrixXd& J, const MatrixXd& Febar, const VectorXd& X, const VectorXd& W, const VectorXi& n, const double& L);
MatrixXd compute_K(const MatrixXd& Sigma, const int& sel, const VectorXd& X, const VectorXd& W, const VectorXi& n, const MatrixXd& Bq, const double& L);


double compute_PassiveOutput(const MatrixXd& Btau_i, const VectorXd& xi, const VectorXd& d_first_i, const VectorXd& X, const VectorXd& W, const double& L);



#endif // SOFTROB_GVS_H
