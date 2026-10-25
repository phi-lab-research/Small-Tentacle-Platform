#ifndef SOFTROB_CCandCLIK_H
#define SOFTROB_CCandCLIK_H

#include "general_fcns.h"
#include "softrob_gvs.h"

void input_collocatedCtrlPsIDFFd_MIXED(
    int n_divisions,
    const VectorXi& n,
    int na,
    const VectorXd& L,
    double tsample,
    const MatrixXd& dist,
    const MatrixXd& m,
    const VectorXd& alpha,
    const VectorXd& qd,
    VectorXd& intThadTha,
    const VectorXd& Rb, 
    const VectorXd& Rt,
    double E,
    double G,
    double rho,
    const VectorXd& mu,
    VectorXd& u,
    const VectorXd& Lc,
    const VectorXd& Lcdot,
    VectorXd& Lcd,
    VectorXd& y_bar,
    VectorXd& yd_bar
);

VectorXd clik_rt(
    int n_divisions,
    const VectorXi& n,
    int na,
    const VectorXd& L,
    const VectorXd& qd_old,
    const VectorXd& xd,
    const VectorXd& xe,
    double tsample,
    const MatrixXd& dist,
    const MatrixXd& m,
    const VectorXd& alpha,
    const VectorXd& Rb, 
    const VectorXd& Rt,
    double E,
    double G,
    double rho,
    const VectorXd& mu,
    bool change_g
);

VectorXd clik_rt_analytic(
    int n_divisions,
    const VectorXi& n,
    int na,
    const VectorXd& L,
    const VectorXd& qd_old,
    const VectorXd& xd,
    const VectorXd& xe,
    double tsample,
    const MatrixXd& dist,
    const MatrixXd& m,
    const VectorXd& alpha,
    const VectorXd& Rb, 
    const VectorXd& Rt,
    double E,
    double G,
    double rho,
    const VectorXd& mu,
    bool change_g
);

MatrixXd compute_dAdq(
    const MatrixXd& t,
    const MatrixXd& dtilde,
    const MatrixXd& tc_normINV,
    const MatrixXd& friction,
    const MatrixXd& Bq,
    int nd,
    const VectorXd& X,
    const VectorXd& Ws,
    const VectorXi& n,
    int na,
    double L);

    MatrixXd compute_dFdq(
    const VectorXd& qd_old,
    const MatrixXd& g,
    const VectorXd& X,
    const VectorXd& Ws,
    const VectorXd& h,
    int nd,
    const VectorXi& n,
    double L,
    const MatrixXd& Mbar,
    const VectorXi& g_pt);

#endif // SOFTROB_CCandCLIK_H
