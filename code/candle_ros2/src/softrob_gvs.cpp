#include "softrob_gvs.h"

 /*MatrixXd compute_Bq1(const double& X,const double& L) {
    MatrixXd Bq1(6, 4);
    Bq1 << 1/L, 0, 0, 0,
             0, 1/L, 0, 0,
             0, 0, 1/L, 0,
             0, 0, 0, 1,
             0, 0, 0, 0,
             0, 0, 0, 0;            
    return Bq1;
}*/

MatrixXd compute_Bq1(const double& X,const double& L) {
    MatrixXd Bq1(6, 3);
    Bq1 <<  0, 0, 0,
            1/L, 0, 0,
            0, 1/L, 0,
            0, 0, 1,
            0, 0, 0,
            0, 0, 0;            
    return Bq1;
}

 MatrixXd compute_Bq2(const double& X,const double& L) {
    MatrixXd Bq2(6, 8);
    Bq2 << 1/L, X/(L*L), 0, 0, 0, 0, 0, 0,
            0, 0, 1/L, X/(L*L), 0, 0, 0, 0,
            0, 0, 0, 0, 1/L, X/(L*L), 0, 0,
            0, 0, 0, 0, 0, 0, 1, X/L,
            0, 0, 0,0, 0, 0, 0, 0,
            0, 0,0, 0, 0, 0, 0, 0;            
    return Bq2;
}
 
 MatrixXd compute_Bq3(const double& X,const double& L) {
    MatrixXd Bq3(6, 4);
    Bq3 << 1/L, 0, 0, 0,
            0, 1/L, 0, 0,
            0, 0, 1/L, 0,
            0, 0, 0, 1,
            0, 0, 0, 0,
            0, 0, 0, 0;              
    return Bq3;
}
 
 MatrixXd compute_Bq4(const double& X,const double& L) {
    MatrixXd Bq4(6, 2);
    Bq4 << 0, 0,
           1/L, 0,
           0, 1/L,
           0, 0,
           0, 0,
           0, 0;            
    return Bq4;
}


VectorXd compute_xi(const VectorXd& q, const double& X, const int& sel, const VectorXi& n, const double& L) {

    VectorXd xi_ref(6);
    xi_ref << 0, 0, 0, 1, 0, 0;

    VectorXd xi = VectorXd::Zero(6);  // Result vector
    MatrixXd Bq;

    Bq.resize(6,n[sel]);
    
    switch (sel) {
        case 0:
            Bq = compute_Bq1(X, L);
            xi = Bq * q.head(n[0]) + xi_ref;  // Extract first n(0) elements of q
            break;
        case 1:
            Bq = compute_Bq2(X, L);
            xi = Bq * q.segment(n[0], n[1]) + xi_ref;  // Extract q from n(0) to n(0)+n(1)
            break;
        case 2:
            Bq = compute_Bq3(X, L);
            xi = Bq * q.segment(n[0] + n[1], n[2]) + xi_ref;  // Extract q from n(0)+n(1) to n(0)+n(1)+n(2)
            break;
        case 3:
            Bq = compute_Bq4(X, L);
            xi = Bq * q.segment(n[0] + n[1] + n[2], n[3]) + xi_ref;  // Extract q from n(0)+n(1)+n(2) to n.sum()
            break;
        default:
            cout<< "Invalid 'sel' value. It must be between 1 and 4." <<endl;
    }

    return xi;
}



VectorXd compute_xidot(const VectorXd& qdot, const double& X, const int& sel, const VectorXi& n, const double& L) {
    // xi_ref initialization
    VectorXd xi_ref(6);
    xi_ref << 0, 0, 0, 1, 0, 0;

    VectorXd xi_dot = VectorXd::Zero(6);  // Result vector
    MatrixXd Bq;
    Bq.resize(6,n[sel]);

    switch (sel) {
        case 0:
            Bq = compute_Bq1(X, L);
            xi_dot = Bq * qdot.head(n[0]);  // Extract first n(0) elements of qdot
            break;
        case 1:
            Bq = compute_Bq2(X, L);
            xi_dot = Bq * qdot.segment(n[0], n[1]);  // Extract qdot from n(0) to n(0)+n(1)
            break;
        case 2:
            Bq = compute_Bq3(X, L);
            xi_dot = Bq * qdot.segment(n[0] + n[1], n[2]);  // Extract qdot from n(0)+n(1) to n(0)+n(1)+n(2)
            break;
        case 3:
            Bq = compute_Bq4(X, L);
            xi_dot = Bq * qdot.segment(n[0] + n[1] + n[2], n[3]);  // Extract qdot from n(0)+n(1)+n(2) to n.sum()
            break;
        default:
            cout<< "Invalid 'sel' value. It must be between 1 and 4." <<endl;
    }

    return xi_dot;
}


MatrixXd compute_T(double theta, const MatrixXd &adjop_omega) {
    double eps_min = 1e-4;
    double eps_max = 1e8;
    
    MatrixXd T(6, 6);
    
    if (theta > eps_max) {
        T = MatrixXd::Identity(6, 6); 
    }
    else if (theta < eps_min) {
        T = MatrixXd::Identity(6, 6) + (1.0 / 2.0) * adjop_omega
            + (1.0 / 6.0) * adjop_omega * adjop_omega
            + (1.0 / 24.0) * adjop_omega * adjop_omega * adjop_omega
            + (1.0 / 120.0) * adjop_omega * adjop_omega * adjop_omega * adjop_omega;
    } else {
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);
        
        double factor1 = (4 - 4 * cos_theta - theta * sin_theta) / (2 * theta * theta);
        double factor2 = (4 * theta - 5 * sin_theta + theta * cos_theta) / (2 * theta * theta * theta);
        double factor3 = (2 - 2 * cos_theta - theta * sin_theta) / (2 * theta * theta * theta * theta);
        double factor4 = (2 * theta - 3 * sin_theta + theta * cos_theta) / (2 * theta * theta * theta * theta * theta);

        T = MatrixXd::Identity(6, 6)
            + factor1 * adjop_omega
            + factor2 * adjop_omega * adjop_omega
            + factor3 * adjop_omega * adjop_omega * adjop_omega
            + factor4 * adjop_omega * adjop_omega * adjop_omega * adjop_omega;
    }
    
    return T;
}

MatrixXd compute_Z(double h, const VectorXd &xi1, const VectorXd &xi2, int sel, double X1, double X2, double L) {
    // Calculate adjoint operators for xi1 and xi2
    MatrixXd adjop_xi1 = adjop_fcn(xi1);
    MatrixXd adjop_xi2 = adjop_fcn(xi2);

    MatrixXd Bq_xi1, Bq_xi2;

    // Switch case for selection (sel) to calculate Bq matrices
    switch (sel) {
        case 0:
            Bq_xi1 = compute_Bq1(X1, L);
            Bq_xi2 = compute_Bq1(X2, L);
            break;
        case 1:
            Bq_xi1 = compute_Bq2(X1, L);
            Bq_xi2 = compute_Bq2(X2, L);
            break;
        case 2:
            Bq_xi1 = compute_Bq3(X1, L);
            Bq_xi2 = compute_Bq3(X2, L);
            break;
        case 3:
            Bq_xi1 = compute_Bq4(X1, L);
            Bq_xi2 = compute_Bq4(X2, L);
            break;
        default:
            cerr << "Invalid selection for sel." << endl;
            return MatrixXd::Zero(6, 6);  // Return a zero matrix in case of error
    }

    // Compute Z based on formula
    MatrixXd Z = (h / 2.0) * (Bq_xi1 + Bq_xi2) + sqrt(3.0) / 12.0 * h * h * (adjop_xi1 * Bq_xi2 - adjop_xi2 * Bq_xi1);
    
    return Z;
}




MatrixXd Mag_fcn(double h, const VectorXd &xi1, const VectorXd &xi2) {

    MatrixXd adjop_xi1 = adjop_fcn(xi1);


    VectorXd omega(6);
    omega = h / 2 * (xi1 + xi2) + sqrt(3.0) * h * h / 12 * (adjop_xi1 * xi2);
  
    MatrixXd omega_hat = HatFromVec_fcn(omega);

    double theta = omega.head(3).norm();  // Norm of the first 3 components of omega
    double eps_min = 1e-4;
    double eps_max = 1e8;

    MatrixXd Mag_exp(4, 4);
    Mag_exp.setZero();  

    if (theta < eps_min) {
        Mag_exp = MatrixXd::Identity(4,4) + omega_hat + 0.5 * omega_hat * omega_hat + (1.0 / 6.0) * omega_hat * omega_hat * omega_hat;
    } else if (theta > eps_max) {
        Mag_exp = MatrixXd::Identity(4,4) + omega_hat;
    } else {
        Mag_exp = MatrixXd::Identity(4,4) + omega_hat + (1.0 - cos(theta)) / (theta * theta) * omega_hat * omega_hat
                             + (theta - sin(theta)) / (theta * theta * theta) * omega_hat * omega_hat * omega_hat;
    }

    return Mag_exp;
}




MatrixXd direct_kin(const MatrixXd &g0, const VectorXd &q, const VectorXd &h, int h_pt, double L,int sel, const VectorXi &n) {
    MatrixXd g = MatrixXd::Zero(4*(h_pt+2),g0.cols());
    g.topRows(4) = g0;  
    double X = 0.0;

    for (int ii = 0; ii < h_pt + 1; ii++) {
        double X1 = X + h[ii] / 2.0 - sqrt(3.0) * h[ii] / 6.0;
        double X2 = X + h[ii] / 2.0 + sqrt(3.0) * h[ii] / 6.0;
        
        // Calculate xi1 and xi2 using compute_xi
        VectorXd xi1 = VectorXd::Zero(6), xi2 = VectorXd::Zero(6);
        xi1 = compute_xi(q, X1 , sel, n, L);
        xi2 = compute_xi(q, X2 , sel, n, L);
        
        // Compute Mag_exp using Mag_fcn
        MatrixXd Mag_exp = Mag_fcn(h[ii], xi1, xi2);

        // Update g using matrix multiplication
        g.block(4 * (ii + 1), 0, 4, g0.cols()) = g.block(4 * (ii), 0, 4, g0.cols()) * Mag_exp;

        // Update X
        X = X + h[ii];
    }

    return g;
}





MatrixXd compute_J(const MatrixXd& J0, const VectorXd &q, const VectorXd &h, int h_pt, double L, int sel, const VectorXi& n) {

    double X = 0.0;
    MatrixXd J = MatrixXd::Zero(6*(h_pt+2),n.sum());
    J.topRows(6) = J0; // Initialize J with J0

    for (int ii = 0; ii <= h_pt; ii++) {
        
        double X1 = X + h[ii] / 2.0 - sqrt(3.0) * h[ii] / 6.0;
        double X2 = X + h[ii] / 2.0 + sqrt(3.0) * h[ii] / 6.0;
        
        VectorXd xi1 = compute_xi(q, X1, sel, n, L);
        VectorXd xi2 = compute_xi(q, X2, sel, n, L);

        MatrixXd Mag_exp = Mag_fcn(h[ii], xi1, xi2);
        MatrixXd adjop_xi1 = adjop_fcn(xi1);
        
        VectorXd omega = h[ii] / 2.0 * (xi1 + xi2) + sqrt(3.0) * pow(h[ii], 2) * (adjop_xi1 * xi2) / 12.0;
        double theta = omega.head(3).norm();  

        MatrixXd adjop_omega = adjop_fcn(omega);
        MatrixXd TOm = compute_T(theta, adjop_omega); // 6x6 matrix
        
        MatrixXd Z = compute_Z(h[ii], xi1, xi2, sel, X1, X2, L); // 6xni matrix

        MatrixXd Meinv = ginv(Mag_exp);
        MatrixXd adjmap_expOmInv = adjmap_fcn(Meinv);


        MatrixXd temp = MatrixXd::Zero(6,n.sum());
        // Depending on `sel`, update `J`
        switch (sel) {
            case 0:
                temp.block(0,0,6,n[0]) = TOm*Z;
                J.block(6 * (ii+1), 0, 6, J.cols()) = adjmap_expOmInv * (J.block(6 * (ii), 0, 6, J.cols()) + temp);
                break;
            case 1:
                temp.block(0,n[0],6,n[1]) = TOm*Z;
                J.block(6 * (ii+1), 0, 6, J.cols()) = adjmap_expOmInv * (J.block(6 * (ii), 0, 6, J.cols()) + temp);
                break;
            case 2:
                temp.block(0,n[0]+n[1],6,n[2]) = TOm*Z;
                J.block(6 * (ii+1), 0, 6, J.cols()) = adjmap_expOmInv * (J.block(6 * (ii), 0, 6, J.cols()) + temp);
                break;
            case 3:
                temp.block(0,n[0]+n[1]+n[2],6,n[3]) = TOm*Z;
                J.block(6 * (ii+1), 0, 6, J.cols()) = adjmap_expOmInv * (J.block(6 * (ii), 0, 6, J.cols()) + temp);
                break;
            default:
                cerr << "Invalid value for sel." << endl;
                break;
        }
        X += h[ii];
    }
    
    return J;
}


MatrixXd compute_intGQ(const int& n, const MatrixXd& IntArg, const VectorXd& X, const VectorXd& W, const double& L) {
    int numCount = X.size();  // size of X (assumed to be the same as size of W)
    MatrixXd Mat = MatrixXd::Zero(n, IntArg.cols());  // Initialize Mat with zeros
    
    for (int count = 0; count < numCount; count++) {
        // Extracting the relevant rows from IntArg
        MatrixXd IntArgSubset = IntArg.block(n * count, 0, n, IntArg.cols());  // block of size n x cols
        Mat += IntArgSubset * W[count] * L;  // Update Mat
    }

    return Mat;
}


MatrixXd compute_AT(const MatrixXd& Btau_i, const MatrixXd& Bq, const VectorXd& X, const VectorXd& W, const double& L) {
    MatrixXd Ai_T = MatrixXd :: Zero(1,Bq.cols());
    MatrixXd Ldfun(X.size(),Bq.cols());  
    int numSections = X.size();

    for (int ii = 0; ii < numSections; ii++) {
        VectorXd Btau_section = Btau_i.block(6 * ii, 0, 6, Btau_i.cols()); 
        MatrixXd Bq_section = Bq.block(6 * ii, 0, 6, Bq.cols());  
        
        Ldfun.row(ii) = Btau_section.transpose() * Bq_section;
    }
    Ai_T = compute_intGQ(1, Ldfun, X, W, L);
    return Ai_T;
}

VectorXd compute_F(const MatrixXd& J, const MatrixXd& Febar, const VectorXd& X, const VectorXd& W, const VectorXi& n, const double& L) {
    int numSections = X.size();

    MatrixXd Ffun(X.size()*n.sum(), Febar.cols());
    Ffun.setZero(); 

    for (int ii = 0; ii < numSections; ii++) {
        Ffun.block(n.sum()*ii, 0, n.sum(), Febar.cols()) = J.block(6*ii, 0, 6, J.cols()).transpose() * Febar.block(6*ii, 0, 6, Febar.cols());
    }

    VectorXd F = compute_intGQ(n.sum(), Ffun, X, W, L);

    return F;
}

MatrixXd compute_K(const MatrixXd& Sigma, const int& sel, const VectorXd& X, const VectorXd& W, const VectorXi& n, const MatrixXd& Bq, const double& L) {
    
    int numSections = X.size();   

    MatrixXd Kfun(n[sel]*X.size(), n[sel]);
    Kfun.setZero(); 
    //cout << "Kfun" << endl << Kfun <<endl; 


    for (int ii = 0; ii < numSections; ii++) {

        MatrixXd partSigma = Sigma.block(6*ii, 6*ii, 6, 6);
        //cout << "partSigma" << endl << partSigma <<endl;

        MatrixXd partBq = Bq.block(6*ii, 0, 6, Bq.cols());
       // cout << "partBq" << endl << partBq <<endl;

        Kfun.block(n[sel]*ii, 0, n[sel], Kfun.cols()) = partBq.transpose() * partSigma * partBq;    
      //  cout << "Kfun block" << endl << Kfun.block(n[sel]*ii, 0, n[sel], Kfun.cols()) <<endl;
    
    }

    MatrixXd K = compute_intGQ(n[sel], Kfun, X, W, L);
    return K;
}

double compute_PassiveOutput(const MatrixXd& Btau_i, 
                            const VectorXd& xi, 
                            const VectorXd& d_first_i, 
                            const VectorXd& X, 
                            const VectorXd& W, 
                            const double& L) {
    
    VectorXd Lfun(X.size()); 
    Lfun.setZero();
    int numSections = X.size();  
    
    for (int ii = 0; ii < numSections; ii++) {
        VectorXd Btau_section = Btau_i.block(6 * ii, 0, 6, Btau_i.cols());
        //cout << "Btau_section " << endl << Btau_section << endl;
    
        VectorXd xi_section = xi.segment(6 * ii, 6);  
        VectorXd d_first_section = d_first_i.segment(3 * (ii), 3); 
        VectorXd temp = VectorXd::Zero(6);
        temp.tail(3) = d_first_section;
        Lfun[ii] = Btau_section.transpose() * (xi_section + temp);
    }

    // Chiamata alla funzione di integrazione
    MatrixXd Lact_mat = compute_intGQ(1, Lfun, X, W, L); 
    double Lact = Lact_mat(0,0);
    
    return Lact;
}
