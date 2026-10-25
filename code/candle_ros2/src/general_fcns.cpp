#include "general_fcns.h"
#include <fstream>
#include <sstream>
#include <vector>


MatrixXd adjmap_fcn(const MatrixXd& mat) {
    MatrixXd ad_m(6, 6);
    ad_m << mat(0,0), mat(0,1), mat(0,2), 0, 0, 0,   
            mat(1,0), mat(1,1), mat(1,2), 0, 0, 0,   
            mat(2,0), mat(2,1), mat(2,2), 0, 0, 0,   
            mat(2,0)*mat(1,3) - mat(1,0)*mat(2,3), mat(2,1)*mat(1,3) - mat(1,1)*mat(2,3), mat(2,2)*mat(1,3) - mat(1,2)*mat(2,3), mat(0,0), mat(0,1), mat(0,2),   
            mat(0,0)*mat(2,3) - mat(2,0)*mat(0,3), mat(0,1)*mat(2,3) - mat(2,1)*mat(0,3), mat(0,2)*mat(2,3) - mat(2,2)*mat(0,3), mat(1,0), mat(1,1), mat(1,2),
            mat(1,0)*mat(0,3) - mat(0,0)*mat(1,3), mat(1,1)*mat(0,3) - mat(0,1)*mat(1,3), mat(1,2)*mat(0,3) - mat(0,2)*mat(1,3), mat(2,0), mat(2,1), mat(2,2);            
    return ad_m;
}

MatrixXd coadjmap_fcn(const MatrixXd& mat) {
    MatrixXd coad_m(6, 6);
    coad_m << mat(0,0), mat(0,1), mat(0,2), mat(2,0)*mat(1,3) - mat(1,0)*mat(2,3), mat(2,1)*mat(1,3) - mat(1,1)*mat(2,3), mat(2,2)*mat(1,3) - mat(1,2)*mat(2,3),   
              mat(1,0), mat(1,1), mat(1,2), mat(0,0)*mat(2,3) - mat(2,0)*mat(0,3), mat(0,1)*mat(2,3) - mat(2,1)*mat(0,3), mat(0,2)*mat(2,3) - mat(2,2)*mat(0,3),   
              mat(2,0), mat(2,1), mat(2,2), mat(1,0)*mat(0,3) - mat(0,0)*mat(1,3), mat(1,1)*mat(0,3) - mat(0,1)*mat(1,3), mat(1,2)*mat(0,3) - mat(0,2)*mat(1,3),   
              0, 0, 0, mat(0,0), mat(0,1), mat(0,2),   
              0, 0, 0, mat(1,0), mat(1,1), mat(1,2),
              0, 0, 0, mat(2,0), mat(2,1), mat(2,2);            
    return coad_m;
}

MatrixXd ginv(const MatrixXd& g) {
    MatrixXd g_inv(4, 4);
    g_inv << g(0,0), g(1,0), g(2,0), -g(0,0)*g(0,3)-g(1,0)*g(1,3)-g(2,0)*g(2,3),
             g(0,1), g(1,1), g(2,1), -g(0,1)*g(0,3)-g(1,1)*g(1,3)-g(2,1)*g(2,3),
             g(0,2), g(1,2), g(2,2), -g(0,2)*g(0,3)-g(1,2)*g(1,3)-g(2,2)*g(2,3),
             0, 0, 0, 1;
    return g_inv;
}


MatrixXd adjop_fcn(const VectorXd& v) {
MatrixXd adj(6,6);  
adj << 0, -v[2], v[1], 0, 0, 0,
       v[2], 0, -v[0], 0, 0, 0,
       -v[1], v[0], 0, 0, 0, 0,
       0, -v[5], v[4], 0, -v[2], v[1],
       v[5], 0, -v[3], v[2], 0, -v[0],
       -v[4], v[3], 0, -v[1], v[0], 0;
  return adj;
}

MatrixXd coadjop_fcn(const VectorXd& v) {
MatrixXd coadj(6,6);  
coadj << 0, -v[2], v[1],0, -v[5], v[4],
         v[2], 0, -v[0],  v[5], 0, -v[3],
         -v[1], v[0], 0, -v[4], v[3], 0,
         0, 0, 0, 0, -v[2], v[1],
         0, 0, 0, v[2], 0, -v[0],
         0, 0, 0, -v[1], v[0], 0;
  return coadj;
}

MatrixXd coadjopbar_fcn(const VectorXd& v) {
    Vector3d v1 = v.segment(0,3);
    Vector3d v2 = v.segment(3,3);

    Matrix3d s11 = skewMat_fcn(v1);
    Matrix3d s12 = skewMat_fcn(v2);
    Matrix3d s21 = s12;
    Matrix3d s22 = Matrix3d::Zero();

    MatrixXd M(6,6);
    M.block(0,0,3,3) = s11;
    M.block(0,3,3,3) = s12;
    M.block(3,0,3,3) = s21;
    M.block(3,3,3,3) = s22;

    return -M;
}


Matrix3d skewMat_fcn(const Vector3d& v) {
    Matrix3d sks_mat;
    sks_mat << 0, -v[2], v[1],
               v[2], 0, -v[0],
               -v[1], v[0], 0;
    return sks_mat;
}

Vector3d skewVec_fcn(const Matrix3d& sks_mat) {
    Vector3d v;
    v << sks_mat(2, 1),  // v(1)
         sks_mat(0, 2),  // v(2)
         sks_mat(1, 0);  // v(3)
    return v;
}

VectorXd VecFromHat_fcn(const MatrixXd&  hat_mat) {
  VectorXd v,v_ang,v_lin;
  Matrix3d sks_mat;
  sks_mat << hat_mat(0,0),hat_mat(0,1),hat_mat(0,2),
             hat_mat(1,0),hat_mat(1,1),hat_mat(1,2),
             hat_mat(2,0),hat_mat(2,1),hat_mat(2,2);
  v_ang=skewVec_fcn(sks_mat);
  v_lin << hat_mat(0,3),hat_mat(1,3),hat_mat(2,3);
  v.resize(6);
  v<<v_ang,v_lin;
  return v;
}


MatrixXd HatFromVec_fcn(const VectorXd& v) {
    MatrixXd hat_mat(4, 4);

    Matrix3d sks_mat = skewMat_fcn(v.head(3));

    hat_mat << sks_mat(0,0),sks_mat(0,1),sks_mat(0,2), v[3],
               sks_mat(1,0),sks_mat(1,1),sks_mat(1,2), v[4],
               sks_mat(2,0),sks_mat(2,1),sks_mat(2,2), v[5],
               0, 0, 0, 0;

    return hat_mat;
}


MatrixXd GaussQuadrature(int N_GQ) {
    double a = 0.0;
    double b = 1.0;

    // Compute the number of Legendre nodes (N-1)
    int N = N_GQ - 1;
    int N1 = N + 1;
    int N2 = N + 2;

    // Initial guess for the Legendre nodes (cosine-based)
    VectorXd xu = VectorXd::LinSpaced(N1, -1.0, 1.0);  // linspace(-1, 1, N1)

    VectorXd y = VectorXd::Zero(N1);
    for (double i = 0; i <= N; i++) {
      y[i] = cos((2*i+1)*M_PI/(2*N+2))+(0.27/N1)*(sin(xu[i]*M_PI*N/N2)); 
    }
    
    // Initialize the Legendre-Gauss Vandermonde Matrix and its derivative
    MatrixXd L = MatrixXd::Zero(N1, N2);
    MatrixXd Lp = VectorXd::Zero(N1);

    // Perform Newton-Raphson iteration to refine the guess of Legendre nodes
    VectorXd y0 = VectorXd::Constant(N1, 2.0);  // Initial "previous" values for convergence
    
    while ((y - y0).cwiseAbs().maxCoeff() > numeric_limits<double>::epsilon()) {


        // First two columns of the Legendre Vandermonde Matrix
        L.col(0) = VectorXd::Ones(N1); 
        L.col(1) = y;
   
        // Compute the Legendre polynomial and its derivative
        for (int k = 1; k < N1; k++) {
            L.col(k+1) = (((2*(k+1) - 1) * y.array() * L.col(k).array() - ((k+1) - 1) * L.col(k-1).array())/(k+1)).matrix();
        }

        ArrayXd Lp_array = N2 * (L.col(N1-1).array() -( y.array() * L.col(N2-1).array())) / ( VectorXd::Ones(N1).array() - y.array().square());
        Lp = Lp_array.matrix();

        // Update y with Newton-Raphson
        y0 = y;
        y = y0 - ( L.col(N2-1).array() / Lp_array).matrix();
  
    }

    // Linear map from [-1, 1] to [a, b]
    VectorXd Xs = ((a * (VectorXd::Ones(N1).array() - y.array())) + (b * (VectorXd::Ones(N1).array()  + y.array()))).matrix() / 2;
    
    // Reverse the order of nodes
    reverse(Xs.data(), Xs.data() + Xs.size());
    //Xs.reverse();

    // Add boundary points (0, 1)
    VectorXd temp=Xs;
    Xs.resize(N1 + 2);
    Xs<<0,temp,1;

    // Compute the weights
    ArrayXd Ws_array = (b - a) / ((VectorXd::Ones(N1).array() - y.array().square()) * Lp.array().square()) * pow(N2/(N1*1.0), 2);
    VectorXd Ws= Ws_array.matrix();
    // Add boundary weights (0)

    temp=Ws;
    Ws.resize(N1 + 2);
    Ws<<0,temp,0;

    // Return results as a matrix containing nodes and weights
    MatrixXd result(N1 + 2, 2);
    result << Xs, Ws;

    return result;
}



MatrixXd loadMatrixFromTxt(const string& filename) {
    ifstream file(filename);
    string line;
    vector<vector<double>> data;

    while (getline(file, line)) {
        stringstream ss(line);
        vector<double> row;
        double val;
        while (ss >> val) {
            row.push_back(val);
        }
        data.push_back(row);
    }

    if (data.empty()) return MatrixXd(); // empty matrix

    int rows = data.size();
    int cols = data[0].size();
    MatrixXd mat(rows, cols);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            mat(i, j) = data[i][j];

    return mat;
}
