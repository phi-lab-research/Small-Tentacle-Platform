#include "softrob_CCandCLIK.h"


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
) {

    double gamma = 1;
    MatrixXd Kp = 800 * MatrixXd::Identity(na, na);
    MatrixXd Kd = 1 * MatrixXd::Identity(na, na);
    MatrixXd Ki = 800 * MatrixXd::Identity(na, na); //1500 ctrl, 800 clik
    VectorXd yd = VectorXd::Zero(na);
    Lcd = VectorXd::Zero(na);

    MatrixXd sum_frcd = MatrixXd::Zero(na,n_divisions);
    MatrixXd sum_mixd = MatrixXd::Zero(na,n_divisions);
   
    MatrixXd frictiond;
    double frcd = 0.0;

    Matrix3d Rot_in;
    Rot_in << 0, 0, 1, 0, 1, 0, -1, 0 ,0;
    MatrixXd g_in = MatrixXd::Zero(4,4);
    g_in.block(0,0,3,3) = Rot_in;
    g_in.block(0,3,3,1) = MatrixXd::Zero(3,1);
    g_in.block(3,0,1,3) = MatrixXd::Zero(1,3);
    g_in(3,3) = 1;
    MatrixXd g0d=MatrixXd::Identity(4, 4)*g_in;

    //MatrixXd g0d=MatrixXd::Identity(4, 4);

    MatrixXd J0d = MatrixXd::Zero(6, n.sum());
    
    VectorXi g_pt(4);
     g_pt << 5, 8, 5, 5;

    // Initialize structures to hold the values
    VectorXd X,Ws,h;
    VectorXd r,A,Jx,Jy,Jz;
    MatrixXd Bq;
    MatrixXd Sigma,Mbar,Fed;
    MatrixXd d_first;
    MatrixXd dtilde;
    VectorXd xid;
    MatrixXd xid_hat;
    VectorXd tcd,td;
    MatrixXd Btaud,Btaud_nofrc;

    MatrixXd Jd= MatrixXd::Zero(6 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2),n.sum());
    MatrixXd gd = MatrixXd::Zero(4 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), 4);
    MatrixXd phi_d = MatrixXd::Zero(g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2, na);
    MatrixXd d = MatrixXd::Zero(3 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), na);

    MatrixXd Ad_T(na,n.sum());
    Ad_T.setZero();


    MatrixXd K(n.sum(),n.sum());
    K.setZero();
    VectorXd Fd(n.sum());
    Fd.setZero();

    VectorXd gravity(6);
    gravity << 0,0,0,0,0,-9.81;


    // Main loop
    for (int nd = 0; nd < n_divisions; nd++) {
        // Gauss quadrature and scaling
        X.resize(g_pt[nd]+2);
        Ws.resize(g_pt[nd]+2);
        h.resize(g_pt[nd]+1);

        MatrixXd res = MatrixXd::Zero(g_pt[nd] + 2, 2);

        res = GaussQuadrature(g_pt[nd]);
        X = res.col(0);
        Ws = res.col(1);
        
        // Scale X by L
        for (int i = 0; i < g_pt[nd]+2; i++) {
            X[i] = X[i] * L[nd];
        }


        for (int i = 0; i < g_pt[nd]+1; i++) {
            h[i] = X[i+1] - X[i];
        }

        switch(nd) {
            case 0:

                gd.block(0,0,4*(g_pt[nd]+2),4)=direct_kin(g0d,qd,h,g_pt[nd],L[nd],nd,n);
                g0d=gd.block(0,0,4*(g_pt[nd]+2),4).bottomRows(4);

                Jd.block(0,0,6*(g_pt[nd]+2),n.sum())=compute_J(J0d,qd,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0d=Jd.block(0,0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
        
                break;
            case 1:
  
                gd.block(4*(g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0d,qd,h,g_pt[nd],L[nd],nd,n);
                g0d=gd.block(4*(g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);   
                
                Jd.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0d,qd,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0d=Jd.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
    
                break;
            case 2:

                gd.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0d,qd,h,g_pt[nd],L[nd],nd,n);
                g0d=gd.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);    
                
                Jd.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0d,qd,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0d=Jd.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
    
                break;
            case 3:
    
                gd.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0d,qd,h,g_pt[nd],L[nd],nd,n);
                g0d=gd.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);
                
                Jd.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0d,qd,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0d=Jd.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);

                break;
        }
        for (int ii = 0; ii < g_pt[nd] + 2; ii++) {
           
            switch (nd) {
                case 0:
                    d.block(ii*3,0,3,na) << 0,0,0,0,
                    (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                    (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
             
                // left-> right:[x;y;z] for each actuator,
                // up->down: i repeat for each point of the section
                break;
                
                case 1:
                    d.block(3*(g_pt[nd-1] + 2)+ii*3,0,3,na) << 0,0,0,0,
                    (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                    (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),(dist(1,nd)+m(0,nd)*X[ii])*sin(alpha[1]),(dist(1,nd)+m(0,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
                 break;
                case 2:
                    d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,0,3,na) << 0,0,0,0,
                    (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                    (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),-(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[1]),-(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
                
                    break;
                case 3:
                    d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,0,3,na) = MatrixXd::Zero(3,na);

                    break;
            }
        }  
    }        
    for (int nd = 0; nd < n_divisions; nd++) {
        // Gauss quadrature and scaling
        X.resize(g_pt[nd]+2);
        Ws.resize(g_pt[nd]+2);
        r.resize(g_pt[nd]+2);
        A.resize(g_pt[nd]+2);
        Jx.resize(g_pt[nd]+2);
        Jy.resize(g_pt[nd]+2);
        Jz.resize(g_pt[nd]+2);
        h.resize(g_pt[nd]+1);
        
        
        MatrixXd res = MatrixXd::Zero(g_pt[nd] + 2, 2);

        res = GaussQuadrature(g_pt[nd]);
        X = res.col(0);
        Ws = res.col(1);
        
        // Scale X by L
        for (int i = 0; i < g_pt[nd]+2; i++) {
            X[i] = X[i] * L[nd];

            r[i]= Rb[nd] + (Rt[nd] - Rb[nd]) * (X[i]/ L[nd]); //X not normalized
            A[i]=M_PI*pow(r[i],2);
            Jx[i]=M_PI/2*pow(r[i],4);
            Jy[i]=M_PI/4*pow(r[i],4);
            Jz[i]=Jy[i];
        }


        for (int i = 0; i < g_pt[nd]+1; i++) {
            h[i] = X[i+1] - X[i];
        }



       
       Bq.resize(6*(g_pt[nd] + 2),n[nd]);
       d_first.resize(3*(g_pt[nd] + 2),na);
       dtilde.resize(3*(g_pt[nd] + 2),3*na);
       xid.resize(6*(g_pt[nd] + 2));
       Btaud.resize(6*(g_pt[nd] + 2),na);
       Btaud_nofrc.resize(6*(g_pt[nd] + 2),na);
       frictiond.resize(1*(g_pt[nd] + 2),na);
  

   
       Sigma.resize(6*(g_pt[nd] + 2),6*(g_pt[nd] + 2));
       Mbar.resize(6*(g_pt[nd] + 2),6*(g_pt[nd] + 2));
       Fed.resize(6*(g_pt[nd] + 2),1);

       Bq.setZero();
       d_first.setZero();
       dtilde.setZero();
       xid.setZero();
       Btaud.setZero();
       Sigma.setZero();
       Mbar.setZero();
       Fed.setZero();
       Btaud_nofrc.setZero();
       frictiond.setZero();

       


        for (int ii = 0; ii < g_pt[nd] + 2; ii++) {    
            VectorXd sigma_el(6),mbar_el(6);
        sigma_el << G*Jx[ii],E*Jy[ii],E*Jz[ii],E*A[ii],G*A[ii],G*A[ii];
        Sigma.block(6*ii,6*ii,6,6)=sigma_el.asDiagonal();
        mbar_el << Jx[ii],Jy[ii],Jz[ii],A[ii],A[ii],A[ii];
        Mbar.block(6*ii,6*ii,6,6)=rho*mbar_el.asDiagonal();

        MatrixXd gnewd = MatrixXd::Zero(4,4);

       
        switch (nd) {
            case 0:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq1(X[ii], L[nd]);
                
                gnewd=gd.block(4*ii,0,4,4);
                
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),m(1,nd)*sin(alpha[1]),m(1,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
          

                break;
            case 1:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq2(X[ii], L[nd]);
                
                gnewd=gd.block(4*(g_pt[nd-1]+2 + ii),0,4,4);
                
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),m(0,nd)*sin(alpha[1]),m(0,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
              
      
                break;
            case 2:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq3(X[ii], L[nd]);

                gnewd=gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,4,4);
                
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),-m(1,nd)*sin(alpha[1]),-m(1,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
              
      
                break;
            case 3:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq4(X[ii], L[nd]);
                
                gnewd=gd.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,4,4);
                
                d_first.block(ii*3,0,3,na) = MatrixXd::Zero(3,na);
                   
                break;
        }

        gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);   
        Fed.block(6*ii,0,6,1)=Mbar.block(6*ii,6*ii,6,6)*coadjmap_fcn(gnewd.inverse())*gravity;    


            
            for (int kk = 0; kk < na; kk++) {
                VectorXd d_homg = VectorXd::Zero(4);
                VectorXd d_first_homg = VectorXd::Zero(4);
            
                switch(nd) {
                    case 0:
                        dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(ii*3,kk,3,1)); 
                        // left-> right:3x3 dtilde for each actuator,
                        // up->down: i repeat for each point of the section
                    
                        //cout <<"dtilde.block(ii*3,kk*3,3,3):"<<endl << dtilde.block(ii*3,kk*3,3,3) << endl;
                    
                        d_homg << d.block(ii*3,kk,3,1),1;                 
                        break;
                    case 1:
                        dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                        d_homg << d.block(3*(g_pt[nd-1] + 2)+ii*3,kk,3,1),1;                 
                        break;
                    case 2:
                        dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                        d_homg << d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1),1;              
                        break;
                    case 3:
                        dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                        d_homg << d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1),1;
                        break;
                }             

			    xid.segment(6*ii,6)=compute_xi(qd,X[ii],nd,n,L[nd]);
               //cout <<"xid:"<<endl << xid << endl;
                xid_hat=HatFromVec_fcn(xid.segment(6*ii,6));
                d_first_homg << d_first.block(ii*3,kk,3,1),1;  
               
                tcd = (xid_hat*d_homg + d_first_homg);
                td= tcd.head(3) / tcd.head(3).norm();

              
                double  th_d,dot_prod;
                Vector3d pd,pnd,fd,fnd;

                //friction model
                if (nd == (n_divisions - 1)){
                        th_d = 0;                    
                }
                else{
                    if(ii == 0){
                        switch(nd){
                            case 0:    
                                pd << 1,0,0;
                                pnd = gd.block(4*(ii+1),3,3,1) - gd.block(4*(ii),3,3,1) + gd.block(4*(ii+1),0,3,3) * d.block(3*(ii+1),kk,3,1) - gd.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);
                                break;
                            case 1:
                                pd =  gd.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - gd.block(4*(g_pt[nd-1]),3,3,1) + gd.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - gd.block(4*(g_pt[nd-1]),0,3,3) * d.block(3*(g_pt[nd-1]),kk,3,1);
                                pnd = gd.block(4*(g_pt[nd-1]+2 + ii+1),3,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + gd.block(4*(g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                                break;
                            case 2:
                                pd =  gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1)   - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),3,3,1) + gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]),kk,3,1);
                                pnd = gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) + gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                           
                                break;                            
                        }
                                    
                        fd=pd/pd.norm();
                        fnd=pnd/pnd.norm();
                        dot_prod = fd.dot(fnd);
                        th_d = acos(max(-1.0, min(dot_prod, 1.0)));
                    }
                    else if (ii == (g_pt[nd] + 1)){
                        switch(nd){
                            case 0:
                                pd = gd.block(4*(ii),3,3,1) - gd.block(4*(ii-1),3,3,1) + gd.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1) - gd.block(4*(ii-1),0,3,3) * d.block(3*(ii-1),kk,3,1);
                                pnd = gd.block(4*(g_pt[nd]+2 +1),3,3,1) - gd.block(4*(ii),3,3,1) + gd.block(4*(g_pt[nd]+2 +1),0,3,3) * d.block(3*(g_pt[nd]+2 +1),kk,3,1) - gd.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);
                                
                                fd=pd/pd.norm();
                                fnd=pnd/pnd.norm();
                                dot_prod = fd.dot(fnd);
                                th_d = acos(max(-1.0, min(dot_prod, 1.0)));
                                break;
                            case 1:
                                pd = gd.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii-1),3,3,1) + gd.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                                pnd = gd.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),3,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + gd.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),kk,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                        
                                fd=pd/pd.norm();
                                fnd=pnd/pnd.norm();
                                dot_prod = fd.dot(fnd);
                                th_d = acos(max(-1.0, min(dot_prod, 1.0)));
                                break;    
                            case 2:
                                pd = gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3,3,1) + gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                               
                                fd=pd/pd.norm();
                                fnd=VectorXd::Zero(3);
                                dot_prod = fd.dot(fnd);
                                th_d = acos(max(-1.0, min(dot_prod, 1.0)));      
                                break;
                        }
                    }
                    else{
                        switch(nd){
                            case 0:
                                pd = gd.block(4*(ii),3,3,1) - gd.block(4*(ii-1),3,3,1) + gd.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1) - gd.block(4*(ii-1),0,3,3) * d.block(3*(ii-1),kk,3,1);
                                pnd = gd.block(4*(ii+1),3,3,1) - gd.block(4*(ii),3,3,1) + gd.block(4*(ii+1),0,3,3) * d.block(3*(ii+1),kk,3,1) - gd.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);
                
                                break;
                            case 1:
                                pd = gd.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii-1),3,3,1) + gd.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                                pnd = gd.block(4*(g_pt[nd-1]+2 + ii+1),3,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + gd.block(4*(g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - gd.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                                
                                break;    
                            case 2:
                                pd = gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3,3,1) + gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                                pnd = gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) + gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - gd.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                           
                                
                                break;
                        }
                        fd=pd/pd.norm();
                        fnd=pnd/pnd.norm();
                        dot_prod = fd.dot(fnd);
                        th_d = acos(max(-1.0, min(dot_prod, 1.0)));
                    }
                }     

               //phi
                  switch(nd) {
                    case 0:
                        phi_d(ii,kk)=th_d;
                        break;
                    case 1:
                        phi_d(g_pt[nd-1]+2 + ii,kk)=th_d;
                        break;
                    case 2:
                        phi_d(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk)=th_d;
                        break;
                    case 3:
                        phi_d(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk)=th_d;
                        break;
                }


                switch(nd) {
                    case 0:
                        frcd = phi_d.col(kk).topRows(ii+1).sum();
                        break;
                    case 1:
                        frcd = phi_d.col(kk).topRows(g_pt[nd-1]+2 + ii+1).sum()- phi_d(g_pt[nd-1]+2,kk);
                        break;
                    case 2:
                        frcd = phi_d.col(kk).topRows(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1).sum()- phi_d(g_pt[nd-2]+2,kk)- phi_d(g_pt[nd-2]+2 + g_pt[nd-1]+2,kk);
                        break;
                    case 3:
                        frcd = 0;
                        break;
                }

                // easier alternative version  
                //  switch(nd) {
                //    case 0:
                //        frcd = X[ii]/(L[0]+L[1]+L[2]);
                //        break;
                //    case 1:
                //        frcd = (L[0]+X[ii])/(L[0]+L[1]+L[2]);
                //        break;
                //    case 2:
                //        frcd = (L[0]+L[1]+X[ii])/(L[0]+L[1]+L[2]);
                //        break;
                //    case 3:
                //        frcd = 0;
                //        break;
                //}
              

                if (nd != (n_divisions - 1)) {
                    frictiond(ii,kk) = exp(-mu[kk] * frcd);
                    Btaud.block(6*ii, kk, 3, 1) = dtilde.block(3*ii,3*kk,3,3)*td*frictiond(ii,kk);
                    Btaud.block(3 + 6*ii, kk, 3, 1) = td*frictiond(ii,kk);
                    Btaud_nofrc.block(6*ii, kk, 3, 1) = dtilde.block(3*ii,3*kk,3,3)*td;
                    Btaud_nofrc.block(3 + 6*ii, kk, 3, 1) = td;
                } else {
                    frictiond(ii,kk) = exp(-mu[kk] * frcd);

                    Btaud.block(6*ii, kk, 6, 1).setZero();
                    Btaud_nofrc.block(6*ii, kk, 6, 1).setZero();
                }
            
            }
                
                
        }
   

        for (int kk = 0; kk < na; kk++) {
            yd [kk] = yd [kk] +  compute_PassiveOutput(Btaud.col(kk),xid,d_first.col(kk),X,Ws,L[nd]);



            for (int ii = 0; ii < g_pt[nd] + 2; ii++) {   
                for (int jj = 0; jj < g_pt[nd] + 2; jj++) {
                    if(ii!=jj){
                        VectorXd d_first_section = d_first.block(3 * (ii),kk, 3,1); 
                        VectorXd zdf_vec = VectorXd::Zero(6);
                        zdf_vec.tail(3) = d_first_section;
                        sum_mixd(kk,nd) = sum_mixd(kk,nd) + (Btaud_nofrc.block(6 * ii, kk, 6, 1).transpose()*( xid.segment(6*ii,6)+zdf_vec)).value()*Ws[ii]*frictiond(jj,kk)*L[nd];
                    }   
                }
            }

            sum_frcd(kk,nd) = frictiond.col(kk).sum();

            Lcd [kk] = Lcd [kk] + (compute_PassiveOutput(Btaud.col(kk),xid,d_first.col(kk),X,Ws,L[nd])+ sum_mixd(kk,nd))/sum_frcd(kk,nd);
             // Lcd [kk] = Lcd [kk] +  compute_PassiveOutput(Btaud.col(kk),xid,d_first.col(kk),X,Ws,L[nd]);
            switch (nd) {
                case 0:
                    Ad_T.block(kk,0,1,n[nd]) = compute_AT(Btaud.col(kk),Bq,X,Ws,L[nd]);                    
                    break;
                case 1:
                    Ad_T.block(kk,n[0],1,n[nd]) = compute_AT(Btaud.col(kk),Bq,X,Ws,L[nd]);                    
                    break;
                case 2:
                    Ad_T.block(kk,n[0]+n[1],1,n[nd]) = compute_AT(Btaud.col(kk),Bq,X,Ws,L[nd]);                    
                    break;
                case 3:
                    Ad_T.block(kk,n[0]+n[1]+n[2],1,n[nd]) = compute_AT(Btaud.col(kk),Bq,X,Ws,L[nd]);                    
                    break;
            }
        }
        switch (nd) {
            case 0:
                K.block(0,0,n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(Jd.block(0,0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                break;
            case 1:
                K.block(n[0],n[0],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(Jd.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                break;
            case 2:
                K.block(n[0]+n[1],n[0]+n[1],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(Jd.block(6*(g_pt[nd-2]+2 + g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                break;
            case 3:
                K.block(n[0]+n[1]+n[2],n[0]+n[1]+n[2],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(Jd.block(6*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                break;
        }
        
    }
    VectorXd F_load = VectorXd::Zero(n.sum());
    VectorXd F_cup = VectorXd::Zero(n.sum());
    VectorXd F_rings = VectorXd::Zero(n.sum());

   
     VectorXd grav_cup = VectorXd::Zero(6);
    //   grav_cup << 0,0,0,0,0,-9.81*0.0016;
     grav_cup << 0,0,0,0,0,-9.81*0.0006;

    //cout << "grav_cup" << grav_cup << endl;
   
    VectorXd grav_ring1 = VectorXd::Zero(6);
    //grav_ring1 << 0,0,0,0,0,-9.81*0.0013;
    grav_ring1 << 0,0,0,0,0,-9.81*0.001;

    VectorXd grav_ring2 = VectorXd::Zero(6);
    //grav_ring2 << 0,0,0,0,0,-9.81*0.0013;
    grav_ring2 << 0,0,0,0,0,-9.81*0.001;
   
    VectorXd grav_ring3 = VectorXd::Zero(6);
    //grav_ring3 << 0,0,0,0,0,-9.81*0.0013;
    grav_ring3 << 0,0,0,0,0,-9.81*0.001;


    MatrixXd gnewd = MatrixXd::Zero(4,4);
    gnewd=gd.bottomRows(4);
    gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
   //F_cup=Jd.bottomRows(6)*coadjmap_fcn(gnewd.inverse())*grav_cup;    
   //cout << "Jd" << endl << Jd.bottomRows(6) << endl;
   //cout << "coadjmap_fcn(gnewd.inverse())"<< endl << coadjmap_fcn(gnewd.inverse()) << endl;
   // 
   
   F_cup=Jd.bottomRows(6).transpose()*coadjmap_fcn(gnewd.inverse())*grav_cup;    

    gnewd=gd.block(4*(g_pt[0]+2 + 4),0,4,4);
    gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
    F_rings = F_rings + Jd.block(6*(g_pt[0]+2 + 4),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring1;    
 
    gnewd=gd.block(4*(g_pt[0]+2 + 6),0,4,4);
    gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
    F_rings = F_rings + Jd.block(6*(g_pt[0]+2 + 6),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring2;    
 
    gnewd=gd.block(4*(g_pt[0]+2 + 9),0,4,4);
    gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
    F_rings = F_rings + Jd.block(6*(g_pt[0]+2 + 9),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring3;    
 
    Fd = Fd + F_load + F_rings + F_cup;

    MatrixXd Jhd = MatrixXd::Zero(Ad_T.rows() + n.sum() - na, Ad_T.cols());
   
    Jhd.block(0,0,Ad_T.rows(),Ad_T.cols()) = Ad_T;
    Jhd.block(Ad_T.rows(),na,n.sum() - na,n.sum()-na) = MatrixXd::Identity(n.sum()-na,n.sum()-na);
   
    VectorXd kthetad = (Jhd.transpose()).inverse() * K * qd;
    VectorXd gthetad = (Jhd.transpose()).inverse() * Fd;

// Calcolo di FFWd
    VectorXd FFWd = kthetad.head(na) - gthetad.head(na);


    // Final terms for u
	 //MatrixXd a1 = MatrixXd::Zero(na, na);
    //a1<< -1.7174,2.4177,0.6246,-1.5001,    
    //-1.8849,2.6948,0.9134,-2.0754,
    //-1.4454,2.5125,0.2784,-1.5384,
    //-1.5297,2.5762,0.7970,-2.0395;
    
    MatrixXd a2 = MatrixXd::Zero(na, na);
    a2<< -0.4685,1.0606,-0.1367,-0.3339,
    0.0999,0.9766,-0.5903,-0.2138,
   -0.1544,0.6599,-0.2230,-0.2185,
   -0.3849,1.0973,0.0567,-0.7438;
    
    MatrixXd b2 = MatrixXd::Zero(na, 1);
    b2 << -0.0692,-0.1520,-0.0378,-0.0353;

    MatrixXd a_dot = MatrixXd::Zero(na, na);
    a_dot<<0.0782,   -0.1791,   -0.2315,    0.2513,
    0.5048,   -0.5460,   -0.5685,    0.4377,
    0.1771,   -0.0913,   -0.2972,    0.1716,
    0.2310,   -0.0506,   -0.2959,    0.0424;

    // y_bar = Lc + a1 * Lc;

    y_bar = Lc + a2 * Lc + b2;
    yd_bar = Lcd + a2 * Lcd + b2;
    VectorXd ydot_bar = Lcdot + a_dot * Lcdot;

    //VectorXd theta_ad=yd_bar;
    //VectorXd theta_a=y_bar;
    //VectorXd thetadot_a=ydot_bar;

    VectorXd theta_ad=Lcd;
    VectorXd theta_a=Lc;
    //VectorXd thetadot_a=Lcdot;
    
    intThadTha = intThadTha + (tanh(theta_ad.array() - theta_a.array()) * tsample).matrix();
  //  u = Kp * (theta_ad - theta_a) - Kd * thetadot_a + Ki / gamma * intThadTha + FFWd;
    u =  Ki / gamma * intThadTha + FFWd;
} 



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
){
    VectorXd qd = VectorXd::Zero(n.sum());

    double deltaq = 0.001;
    MatrixXd Gamma_Z = 0.1*MatrixXd::Identity(na, na); //1 slow, 0.1 fast
    MatrixXd Gamma_eps = 0.1*MatrixXd::Identity(n.sum()-na,n.sum()- na);  //1 slow, 0.1 fast


    if (change_g)
    {
        Gamma_Z = 0.3*MatrixXd::Identity(na, na);
        Gamma_eps = 0.3*MatrixXd::Identity(n.sum()-na,n.sum()- na);
    }
  
   
   // Initialize structures to hold the values
   VectorXd X,Ws,h;
   VectorXd r, A, Jx, Jy, Jz;
   MatrixXd Bq;
   MatrixXd d_first,dtilde;
   VectorXd xi;
   MatrixXd xi_hat;
   VectorXd tc;
   VectorXd t;
   MatrixXd xi_dq,xi_hat_dq;
   MatrixXd Btau,Btau_dq;
   MatrixXd friction,friction_dq;
   MatrixXd Sigma, Mbar, Fed, Fed_dq;

    VectorXd y = VectorXd::Zero(na);
    MatrixXd y_dq = MatrixXd::Zero(na, n.sum());

    MatrixXd t_dq = MatrixXd::Zero(3, n.sum());     
    MatrixXd tc_dq = MatrixXd::Zero(4, n.sum());

   double frc = 0.0;
   MatrixXd frc_dq = MatrixXd::Zero(1, n.sum());


   Matrix3d Rot_in;
   Rot_in << 0, 0, 1, 0, 1, 0, -1, 0 ,0;
   MatrixXd g_in = MatrixXd::Zero(4,4);
   g_in.block(0,0,3,3) = Rot_in;
   g_in.block(0,3,3,1) = MatrixXd::Zero(3,1);
   g_in.block(3,0,1,3) = MatrixXd::Zero(1,3);
   g_in(3,3) = 1;
   MatrixXd g0=MatrixXd::Identity(4, 4)*g_in;
   MatrixXd J0=MatrixXd::Zero(6, n.sum());
   
   MatrixXd g0_dq=MatrixXd::Zero(4, 4*n.sum());
   MatrixXd J0_dq=MatrixXd::Zero(6, n.sum()*n.sum());

   for (int comp=0;comp<n.sum();comp++){
       g0_dq.block(0,4*comp,4,4)=g0;
       J0_dq.block(0,comp*n.sum(),6,n.sum())=J0;
   }


   VectorXi g_pt(4);
   g_pt << 5, 8, 5, 5;

   MatrixXd g = MatrixXd::Zero(4 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), 4);
   MatrixXd phi = MatrixXd::Zero(g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2, na);
   MatrixXd J = MatrixXd::Zero(6 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), n.sum());
   MatrixXd g_dq = MatrixXd::Zero(4 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), 4*n.sum());
   MatrixXd phi_dq = MatrixXd::Zero(g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2, na*n.sum());
   MatrixXd J_dq = MatrixXd::Zero(6 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), n.sum()*n.sum());
   MatrixXd d = MatrixXd::Zero(3 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), na);

   MatrixXd A_T(na, n.sum());
   A_T.setZero();
   MatrixXd A_dq_T(na, n.sum()*n.sum());
   A_dq_T.setZero();

   MatrixXd K(n.sum(), n.sum());
   K.setZero();
   VectorXd Fd(n.sum());
   Fd.setZero();
   MatrixXd Fd_dq(n.sum(),n.sum());
   Fd_dq.setZero();

   VectorXd gravity(6);
   gravity << 0, 0, 0, 0, 0, -9.81;

   MatrixXd comp_sel = MatrixXd::Identity(n.sum(),n.sum());
    

   for (int nd = 0; nd < n_divisions; nd++) {
    X.resize(g_pt[nd]+2);
    Ws.resize(g_pt[nd]+2);
    h.resize(g_pt[nd]+1);
    
    
    MatrixXd res = MatrixXd::Zero(g_pt[nd] + 2, 2);

    res = GaussQuadrature(g_pt[nd]);
    X = res.col(0);
    Ws = res.col(1);
    
    // Scale X by L
    for (int i = 0; i < g_pt[nd]+2; i++) {
        X[i] = X[i] * L[nd];

    }


    for (int i = 0; i < g_pt[nd]+1; i++) {
        h[i] = X[i+1] - X[i];
    }

    switch(nd) {
        case 0:
            g.block(0,0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(0,0,4*(g_pt[nd]+2),4).bottomRows(4);                
            
            J.block(0,0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(0,0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
    
            for (int comp=0;comp<n.sum();comp++){       
                g_dq.block(0,4*comp,4*(g_pt[nd]+2),4)=direct_kin(g0_dq.block(0,4*comp,4,4),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n);
                g0_dq.block(0,4*comp,4,4)=g_dq.block(0,4*comp,4*(g_pt[nd]+2),4).bottomRows(4);
              
                J_dq.block(0,comp*n.sum(),6*(g_pt[nd]+2),n.sum())=compute_J(J0_dq.block(0,comp*n.sum(),6,n.sum()),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0_dq.block(0,comp*n.sum(),6,n.sum())=J_dq.block(0,comp*n.sum(),6*(g_pt[nd]+2),n.sum()).bottomRows(6);
            }
            break;


        case 1:
            g.block(4*(g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(4*(g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);    
            
            J.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);


            for (int comp=0;comp<n.sum();comp++){       
                g_dq.block(4*(g_pt[nd-1]+2),4*comp,4*(g_pt[nd]+2),4)=direct_kin(g0_dq.block(0,4*comp,4,4),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n);
                g0_dq.block(0,4*comp,4,4)=g_dq.block(4*(g_pt[nd-1]+2),4*comp,4*(g_pt[nd]+2),4).bottomRows(4);
            
                J_dq.block(6*(g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum())=compute_J(J0_dq.block(0,comp*n.sum(),6,n.sum()),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0_dq.block(0,comp*n.sum(),6,n.sum())=J_dq.block(6*(g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum()).bottomRows(6);
            }

            break;
        case 2:
            g.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);         

            J.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
            
            for (int comp=0;comp<n.sum();comp++){       
                g_dq.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),4*comp,4*(g_pt[nd]+2),4)=direct_kin(g0_dq.block(0,4*comp,4,4),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n);
                g0_dq.block(0,4*comp,4,4)=g_dq.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),4*comp,4*(g_pt[nd]+2),4).bottomRows(4);
            
                J_dq.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum())=compute_J(J0_dq.block(0,comp*n.sum(),6,n.sum()),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0_dq.block(0,comp*n.sum(),6,n.sum())=J_dq.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum()).bottomRows(6);
            }

            break;
        case 3:

            g.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);     
   
            
            J.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
            for (int comp=0;comp<n.sum();comp++){       
                g_dq.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),4*comp,4*(g_pt[nd]+2),4)=direct_kin(g0_dq.block(0,4*comp,4,4),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n);
                g0_dq.block(0,4*comp,4,4)=g_dq.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),4*comp,4*(g_pt[nd]+2),4).bottomRows(4);
            
                J_dq.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum())=compute_J(J0_dq.block(0,comp*n.sum(),6,n.sum()),qd_old+comp_sel.col(comp)*deltaq,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
                J0_dq.block(0,comp*n.sum(),6,n.sum())=J_dq.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum()).bottomRows(6);
            }
            break;
    }
    for (int ii = 0; ii < g_pt[nd] + 2; ii++) {
           
        switch (nd) {
            case 0:
                d.block(ii*3,0,3,na) << 0,0,0,0,
                (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
         
            // left-> right:[x;y;z] for each actuator,
            // up->down: i repeat for each point of the section
            break;
            
            case 1:
                d.block(3*(g_pt[nd-1] + 2)+ii*3,0,3,na) << 0,0,0,0,
                (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),(dist(1,nd)+m(0,nd)*X[ii])*sin(alpha[1]),(dist(1,nd)+m(0,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
             break;
            case 2:
                d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,0,3,na) << 0,0,0,0,
                (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),-(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[1]),-(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
            
                break;
            case 3:
                d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,0,3,na) = MatrixXd::Zero(3,na);

                break;
        }
    }     
}  
for (int nd = 0; nd < n_divisions; nd++) {
    // Gauss quadrature and scaling
    X.resize(g_pt[nd]+2);
    Ws.resize(g_pt[nd]+2);
    r.resize(g_pt[nd]+2);
    A.resize(g_pt[nd]+2);
    Jx.resize(g_pt[nd]+2);
    Jy.resize(g_pt[nd]+2);
    Jz.resize(g_pt[nd]+2);
    h.resize(g_pt[nd]+1);
    
    
    MatrixXd res = MatrixXd::Zero(g_pt[nd] + 2, 2);

    res = GaussQuadrature(g_pt[nd]);
    X = res.col(0);
    Ws = res.col(1);
    
    // Scale X by L
    for (int i = 0; i < g_pt[nd]+2; i++) {
        X[i] = X[i] * L[nd];

        r[i]= Rb[nd] + (Rt[nd] - Rb[nd]) * (X[i]/ L[nd]); //X not normalized
        A[i]=M_PI*pow(r[i],2);
        Jx[i]=M_PI/2*pow(r[i],4);
        Jy[i]=M_PI/4*pow(r[i],4);
        Jz[i]=Jy[i];
    }


    for (int i = 0; i < g_pt[nd]+1; i++) {
        h[i] = X[i+1] - X[i];
    }

    
    Bq.resize(6*(g_pt[nd] + 2),n[nd]);
    d_first.resize(3*(g_pt[nd] + 2),na);
    dtilde.resize(3*(g_pt[nd] + 2),3*na);
    xi.resize(6*(g_pt[nd] + 2));
    Btau.resize(6*(g_pt[nd] + 2),na);
    friction.resize(1*(g_pt[nd] + 2),na);
    
    Sigma.resize(6*(g_pt[nd] + 2),6*(g_pt[nd] + 2));
    Mbar.resize(6*(g_pt[nd] + 2),6*(g_pt[nd] + 2));
    Fed.resize(6*(g_pt[nd] + 2),1);
    Fed_dq.resize(6*(g_pt[nd] + 2),n.sum());
    friction_dq.resize(1*(g_pt[nd] + 2),na*n.sum());
    Btau_dq.resize(6*(g_pt[nd] + 2),na*n.sum());
    xi_dq.resize(6*(g_pt[nd] + 2),n.sum());

    Bq.setZero();
    d_first.setZero();
    dtilde.setZero();
    xi.setZero();
    Btau.setZero();
    Sigma.setZero();
    Mbar.setZero();
    Fed.setZero();
    friction.setZero();
    friction_dq.setZero();
    Btau_dq.setZero();
    xi_dq.setZero();
    Fed_dq.setZero();
    
    
    // Loop for each Gauss point and division
    for (int ii = 0; ii < g_pt[nd] + 2; ii++) {
        VectorXd sigma_el(6),mbar_el(6);
        sigma_el << G*Jx[ii],E*Jy[ii],E*Jz[ii],E*A[ii],G*A[ii],G*A[ii];
        Sigma.block(6*ii,6*ii,6,6)=sigma_el.asDiagonal();
        mbar_el << Jx[ii],Jy[ii],Jz[ii],A[ii],A[ii],A[ii];
        Mbar.block(6*ii,6*ii,6,6)=rho*mbar_el.asDiagonal();

        MatrixXd gnew = MatrixXd::Zero(4,4);
        MatrixXd gnew_dq = MatrixXd::Zero(4,4*n.sum());

        switch (nd) {
            case 0:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq1(X[ii], L[nd]);
                
                gnew=g.block(4*ii,0,4,4);


                for (int comp=0;comp<n.sum();comp++){       
                    gnew_dq.block(0,4*comp,4,4)=g_dq.block(4*ii,4*comp,4,4);
                }
                
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),m(1,nd)*sin(alpha[1]),m(1,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
        
          

                break;
            case 1:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq2(X[ii], L[nd]);
                
                gnew=g.block(4*(g_pt[nd-1]+2 + ii),0,4,4);

                for (int comp=0;comp<n.sum();comp++){       
                    gnew_dq.block(0,4*comp,4,4)=g_dq.block(4*(g_pt[nd-1]+2 + ii),4*comp,4,4);
                }
                
                            
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),m(0,nd)*sin(alpha[1]),m(0,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
            
                break;
            case 2:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq3(X[ii], L[nd]);

                gnew=g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,4,4);
                for (int comp=0;comp<n.sum();comp++){       
                    gnew_dq.block(0,4*comp,4,4)=g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),4*comp,4,4);
                }
                 
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),-m(1,nd)*sin(alpha[1]),-m(1,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
           
          
                break;
            case 3:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq4(X[ii], L[nd]);
                
                gnew=g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,4,4);
                for (int comp=0;comp<n.sum();comp++){       
                    gnew_dq.block(0,4*comp,4,4)=g_dq.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),4*comp,4,4);
                }
                
                d_first.block(ii*3,0,3,na) = MatrixXd::Zero(3,na);
                   
                break;
        }

        gnew.block(0,3,3,1)=VectorXd::Zero(3,1);   

        Fed.block(6*ii,0,6,1)=Mbar.block(6*ii,6*ii,6,6)*coadjmap_fcn(gnew.inverse())*gravity;    

        for (int comp=0;comp<n.sum();comp++){       
            gnew_dq.block(0,3+4*comp,3,1)=VectorXd::Zero(3,1);
            Fed_dq.block(6*ii,comp,6,1)=Mbar.block(6*ii,6*ii,6,6)*coadjmap_fcn(gnew_dq.block(0,4*comp,4,4).inverse())*gravity;    
        }


        for (int kk = 0; kk < na; kk++) {

            VectorXd d_homg = VectorXd::Zero(4);
            VectorXd d_first_homg = VectorXd::Zero(4);

            switch(nd) {
                case 0:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(ii*3,kk,3,1)); 
                    // left-> right:3x3 dtilde for each actuator,
                    // up->down: i repeat for each point of the section
                
                    //cout <<"dtilde.block(ii*3,kk*3,3,3):"<<endl << dtilde.block(ii*3,kk*3,3,3) << endl;
                
                    d_homg << d.block(ii*3,kk,3,1),1;                 
                    break;
                case 1:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                    d_homg << d.block(3*(g_pt[nd-1] + 2)+ii*3,kk,3,1),1;                 
                    break;
                case 2:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                    d_homg << d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1),1;              
                    break;
                case 3:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                    d_homg << d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1),1;
                    break;
            }
            xi.segment(6*ii,6)=compute_xi(qd_old,X[ii],nd,n,L[nd]);    
            xi_hat=HatFromVec_fcn(xi.segment(6*ii,6));

            d_first_homg << d_first.block(ii*3,kk,3,1),1;
        
            tc = (xi_hat*d_homg + d_first_homg);
            t = tc.head(3)/tc.head(3).norm();     

            for(int comp = 0; comp < n.sum(); comp++) {
                xi_dq.block(6*ii,comp,6,1)=compute_xi(qd_old+comp_sel.col(comp)*deltaq,X[ii],nd,n,L[nd]);
                xi_hat_dq=HatFromVec_fcn(xi_dq.block(6*ii,comp,6,1));

                tc_dq.col(comp) = (xi_hat_dq*d_homg + d_first_homg);
                t_dq.col(comp) = tc_dq.col(comp).head(3)/tc_dq.col(comp).head(3).norm();       
            }

            double  th,dot_prod;
            Vector3d p,pn,f,fn;
            MatrixXd th_dq = MatrixXd::Zero(1, n.sum());
            MatrixXd p_dq = MatrixXd::Zero(3, n.sum());
            MatrixXd pn_dq = MatrixXd::Zero(3, n.sum()); 
            MatrixXd f_dq = MatrixXd::Zero(3, n.sum());
            MatrixXd fn_dq = MatrixXd::Zero(3, n.sum());
             

            //friction model
            if (nd == (n_divisions - 1)){
                for(int comp = 0; comp < n.sum(); comp++) {
                    th_dq(0,comp) = 0;       
                }                 
                th =0; 
            }
            else{
                if(ii == 0){
                    switch(nd){
                        case 0:    
                            p << 1,0,0;
                            pn = g.block(4*(ii+1),3,3,1) - g.block(4*(ii),3,3,1) + g.block(4*(ii+1),0,3,3) * d.block(3*(ii+1),kk,3,1) - g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);
                        
                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) << 1,0,0;   
                                pn_dq.col(comp) = g_dq.block(4*(ii+1),3+4*comp,3,1) - g_dq.block(4*(ii),3+4*comp,3,1) + g_dq.block(4*(ii+1),4*comp,3,3) * d.block(3*(ii+1),kk,3,1) - g_dq.block(4*(ii),4*comp,3,3) * d.block(3*(ii),kk,3,1);    
                            } 

                            break;
                        case 1:
                            p =  g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-1]),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-1]),0,3,3) * d.block(3*(g_pt[nd-1]),kk,3,1);
                            pn = g.block(4*(g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                            
                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) = g_dq.block(4*(g_pt[nd-1]+2 + ii),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-1]),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g_dq.block(4*(g_pt[nd-1]),4*comp,3,3) * d.block(3*(g_pt[nd-1]),kk,3,1);
                                pn_dq.col(comp) = g_dq.block(4*(g_pt[nd-1]+2 + ii+1),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-1]+2 + ii+1),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);   
                            } 
                             break;
                        case 2:
                            p =  g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1)   - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]),kk,3,1);
                            pn = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                           
                            
                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) = g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3+4*comp,3,1)   - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),comp*4,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]),kk,3,1);
                                pn_dq.col(comp) = g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),4*comp,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),comp*4,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                              
                            } 
                            break;                            
                    }

                    f=p/p.norm();
                    fn=pn/pn.norm();
                    dot_prod = f.dot(fn);
                    th = acos(max(-1.0, min(dot_prod, 1.0)));


                    for(int comp = 0; comp < n.sum(); comp++) {
                        f_dq.col(comp) = p_dq.col(comp)/p_dq.col(comp).norm();
                        fn_dq.col(comp) = pn_dq.col(comp)/pn_dq.col(comp).norm();
                        dot_prod = f_dq.col(comp).dot(fn_dq.col(comp));
                        th_dq(0,comp) = acos(max(-1.0, min(dot_prod, 1.0)));
                    }
                }
                else if (ii == (g_pt[nd] + 1)){
                    
                    switch(nd){
                        case 0:
                          
                            for(int comp = 0; comp < n.sum(); comp++) {

                                p_dq.col(comp) = g_dq.block(4*(ii),3+4*comp,3,1) - g_dq.block(4*(ii-1),3+4*comp,3,1) + g_dq.block(4*(ii),4*comp,3,3) * d.block(3*(ii),kk,3,1) - g_dq.block(4*(ii-1),4*comp,3,3) * d.block(3*(ii-1),kk,3,1);
                                pn_dq.col(comp) = g_dq.block(4*(g_pt[nd]+2 +1),3+4*comp,3,1) - g_dq.block(4*(ii),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd]+2 +1),4*comp,3,3) * d.block(3*(g_pt[nd]+2 +1),kk,3,1) - g_dq.block(4*(ii),4*comp,3,3) * d.block(3*(ii),kk,3,1);    
                          
                                f_dq.col(comp) = p_dq.col(comp)/p_dq.col(comp).norm();
                                fn_dq.col(comp) = pn_dq.col(comp)/pn_dq.col(comp).norm();
                                dot_prod = f_dq.col(comp).dot(fn_dq.col(comp));
                                th_dq(0,comp) = acos(max(-1.0, min(dot_prod, 1.0)));
                            }

                            p = g.block(4*(ii),3,3,1) - g.block(4*(ii-1),3,3,1) + g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1) - g.block(4*(ii-1),0,3,3) * d.block(3*(ii-1),kk,3,1);
                            pn = g.block(4*(g_pt[nd]+2 +1),3,3,1) - g.block(4*(ii),3,3,1) + g.block(4*(g_pt[nd]+2 +1),0,3,3) * d.block(3*(g_pt[nd]+2 +1),kk,3,1) - g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);
                            
                            f=p/p.norm();
                            fn=pn/pn.norm();
                            dot_prod = f.dot(fn);
                            th = acos(max(-1.0, min(dot_prod, 1.0)));
                        
                            break;
                        case 1:

                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) = g_dq.block(4*(g_pt[nd-1]+2 + ii),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii-1),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii-1),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                                pn_dq.col(comp) = g_dq.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),kk,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);   
                                f_dq.col(comp) = p_dq.col(comp)/p_dq.col(comp).norm();
                                fn_dq.col(comp) = pn_dq.col(comp)/pn_dq.col(comp).norm();
                                dot_prod = f_dq.col(comp).dot(fn_dq.col(comp));
                                th_dq(0,comp) = acos(max(-1.0, min(dot_prod, 1.0)));
                            } 


                            p = g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                            pn = g.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                    
                            f=p/p.norm();
                            fn=pn/pn.norm();
                            dot_prod = f.dot(fn);
                            th = acos(max(-1.0, min(dot_prod, 1.0)));
                            break;    
                        case 2:

                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) = g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),4*comp,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                          
                                f_dq.col(comp) = p_dq.col(comp)/p_dq.col(comp).norm();
                                fn_dq.col(comp) = VectorXd::Zero(3);
                                dot_prod = f_dq.col(comp).dot(fn_dq.col(comp));
                                th_dq(0,comp) = acos(max(-1.0, min(dot_prod, 1.0)));
                           
                           
                            } 

                            p = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                        
                            f=p/p.norm();
                            fn=VectorXd::Zero(3);
                            dot_prod = f.dot(fn);
                            th = acos(max(-1.0, min(dot_prod, 1.0)));   
                            break;
                    }
                }
                else{
                    switch(nd){
                        case 0:

                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) = g_dq.block(4*(ii),3+4*comp,3,1) - g_dq.block(4*(ii-1),3+4*comp,3,1) + g_dq.block(4*(ii),4*comp,3,3) * d.block(3*(ii),kk,3,1) - g_dq.block(4*(ii-1),4*comp,3,3) * d.block(3*(ii-1),kk,3,1);
                                pn_dq.col(comp) = g_dq.block(4*(ii+1),3+4*comp,3,1) - g_dq.block(4*(ii),3+4*comp,3,1) + g_dq.block(4*(ii+1),4*comp,3,3) * d.block(3*(ii+1),kk,3,1) - g_dq.block(4*(ii),4*comp,3,3) * d.block(3*(ii),kk,3,1);    
                            } 

                            p = g.block(4*(ii),3,3,1) - g.block(4*(ii-1),3,3,1) + g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1) - g.block(4*(ii-1),0,3,3) * d.block(3*(ii-1),kk,3,1);
                            pn = g.block(4*(ii+1),3,3,1) - g.block(4*(ii),3,3,1) + g.block(4*(ii+1),0,3,3) * d.block(3*(ii+1),kk,3,1) - g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);

                            break;
                        case 1:
                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) = g_dq.block(4*(g_pt[nd-1]+2 + ii),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii-1),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii-1),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                                pn_dq.col(comp) = g_dq.block(4*(g_pt[nd-1]+2 + ii+1),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-1]+2 + ii+1),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - g_dq.block(4*(g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);   
                            }   
                            p = g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                            pn = g.block(4*(g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                            
                            break;    
                        case 2:
                            for(int comp = 0; comp < n.sum(); comp++) {
                                p_dq.col(comp) = g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),4*comp,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),4*comp,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                                pn_dq.col(comp) = g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3+4*comp,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3+4*comp,3,1) + g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),4*comp,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - g_dq.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),comp*4,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                              
                            } 
                           

                            p = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                            pn = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                           
                            
                            break;
                    }
                    for(int comp = 0; comp < n.sum(); comp++) {
                        f_dq.col(comp) = p_dq.col(comp)/p_dq.col(comp).norm();
                        fn_dq.col(comp) = pn_dq.col(comp)/pn_dq.col(comp).norm();
                        dot_prod = f_dq.col(comp).dot(fn_dq.col(comp));
                        th_dq(0,comp) = acos(max(-1.0, min(dot_prod, 1.0)));
                    }
                
                    f=p/p.norm();
                    fn=pn/pn.norm();
                    dot_prod = f.dot(fn);
                    th = acos(max(-1.0, min(dot_prod, 1.0)));
                }
            }    


            //phi
            switch(nd) {
                case 0:
                    phi(ii,kk)=th;
                    for(int comp=0; comp < n.sum(); comp++) {
                        phi_dq(ii,kk+na*comp)=th_dq(0,comp);
                    }
                    break;
                case 1:
                    phi(g_pt[nd-1]+2 + ii,kk)=th;
                    for (int comp=0; comp < n.sum(); comp++) {
                        phi_dq(g_pt[nd-1]+2 + ii,kk+na*comp)=th_dq(0,comp);
                    }
                    break;
                case 2:
                    phi(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk)=th;
                    for(int comp=0; comp < n.sum(); comp++) {
                        phi_dq(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk+na*comp)=th_dq(0,comp);
                    }
                    
                    break;
                case 3:
                    phi(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk)=th;
                    for(int comp=0; comp < n.sum(); comp++) {
                        phi_dq(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk+na*comp)=th_dq(0,comp);
                    }
                    break;
            }


            switch(nd) {
                case 0:
                    frc = phi.col(kk).topRows(ii+1).sum();
                    for(int comp=0; comp < n.sum(); comp++) {
                        frc_dq(0,comp) = phi_dq.col(kk+na*comp).topRows(ii+1).sum();                       
                    }

                    break;
                case 1:
                    frc = phi.col(kk).topRows(g_pt[nd-1]+2 + ii+1).sum()-phi(g_pt[nd-1]+2,kk);  
                    for (int comp=0; comp < n.sum(); comp++) {
                        frc_dq(0,comp) = phi_dq.col(kk+na*comp).topRows(g_pt[nd-1]+2 + ii+1).sum()-phi_dq(g_pt[nd-1]+2,kk+na*comp);                       
                    }
                    break;
                case 2:
                    frc = phi.col(kk).topRows(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1).sum()-phi(g_pt[nd-2]+2,kk)- phi(g_pt[nd-2]+2 + g_pt[nd-1]+2,kk);
                    for (int comp=0; comp < n.sum(); comp++) {
                        frc_dq(0,comp) = phi_dq.col(kk+na*comp).topRows(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1).sum()-phi_dq(g_pt[nd-2]+2,kk+na*comp)- phi_dq(g_pt[nd-2]+2 + g_pt[nd-1]+2,kk+na*comp);                       
                    }
                    break;
                case 3:
                    frc = 0;
                    for (int comp=0; comp < n.sum(); comp++) {
                        frc_dq(0,comp) = 0;                       
                    }
                    break;
            }

            // easier alternative version  
            //  switch(nd) {
            //    case 0:
            //        frc = X[ii]/(L[0]+L[1]+L[2]);
            //        break;
            //    case 1:
            //        frc = (L[0]+X[ii])/(L[0]+L[1]+L[2]);
            //        break;
            //    case 2:
            //        frc = (L[0]+L[1]+X[ii])/(L[0]+L[1]+L[2]);
            //        break;
            //    case 3:
            //        frc = 0;
            //        break;
            //}
            //for (int comp=0; comp < n.sum(); comp++) {
            //    frcd_dq(0,comp) = frc;
            //}

            if (nd != (n_divisions - 1)) {
                // Calcolo dell'attrito
                friction(ii,kk) = exp(-mu[kk] * frc);
                for (int comp=0; comp < n.sum(); comp++) {
                    friction_dq(ii,kk+na*comp) = exp(-mu[kk] * frc_dq(0,comp));
                }

                Btau.block(6*ii, kk, 3, 1) = dtilde.block(3*ii,3*kk,3,3)*t*friction(ii,kk);
                Btau.block(3 + 6*ii, kk, 3, 1) = t*friction(ii,kk);
                 
                for (int comp=0; comp < n.sum(); comp++) {
                    Btau_dq.block(6*ii, kk+na*comp, 3, 1) = dtilde.block(3*ii,3*kk,3,3)*t_dq.col(comp)*friction_dq(ii,kk+na*comp);
                    Btau_dq.block(3 + 6*ii, kk+na*comp, 3, 1) = t_dq.col(comp)*friction_dq(ii,kk+na*comp);
                }
             

            } else { 
                friction(ii,kk) = exp(-mu[kk] * frc);
                for (int comp=0; comp < n.sum(); comp++) {
                    friction_dq(ii,kk+na*comp) = exp(-mu[kk] * frc_dq(0,comp));
                }

                Btau.block(6*ii, kk, 6, 1).setZero();
                for (int comp=0; comp < n.sum(); comp++) {
                    Btau_dq.block(6*ii, kk+na*comp, 6, 1).setZero();
                }
            }
        }
    }
        for (int kk = 0; kk < na; kk++) {
            y [kk] = y [kk] + compute_PassiveOutput(Btau.col(kk),xi,d_first.col(kk),X,Ws,L[nd]);
           
            for (int comp=0; comp < n.sum(); comp++) {
                y_dq(kk,comp) = y_dq(kk,comp) + compute_PassiveOutput(Btau_dq.col(kk+na*comp),xi_dq.col(comp),d_first.col(kk),X,Ws,L[nd]);
            }
   

            switch (nd) {
                case 0:
                    A_T.block(kk,0,1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);

                    for (int comp=0; comp < n.sum(); comp++) {
                       A_dq_T.block(kk,comp*n.sum(),1,n[nd]) = compute_AT(Btau_dq.col(kk+na*comp),Bq,X,Ws,L[nd]);
                    }                  
                
                  break;
                case 1:
                    A_T.block(kk,n[0],1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);
                    for (int comp=0; comp < n.sum(); comp++) {
                        A_dq_T.block(kk,n[0]+comp*n.sum(),1,n[nd]) = compute_AT(Btau_dq.col(kk+na*comp),Bq,X,Ws,L[nd]);
                    }                
                   
                    break;
                case 2:
                    A_T.block(kk,n[0]+n[1],1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);
                    for (int comp=0; comp < n.sum(); comp++) {
                        A_dq_T.block(kk,n[0]+n[1]+comp*n.sum(),1,n[nd]) = compute_AT(Btau_dq.col(kk+na*comp),Bq,X,Ws,L[nd]);
                    }
                    break;
                case 3:
                    A_T.block(kk,n[0]+n[1]+n[2],1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);
                    for (int comp=0; comp < n.sum(); comp++) {
                        A_dq_T.block(kk,n[0]+n[1]+n[2]+comp*n.sum(),1,n[nd]) = compute_AT(Btau_dq.col(kk+na*comp),Bq,X,Ws,L[nd]);
                    }
                    break;
            }
        }
        switch (nd) {
            case 0:
                K.block(0,0,n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(0,0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                for (int comp=0; comp < n.sum(); comp++) {
                    Fd_dq.col(comp) = Fd_dq.col(comp) + compute_F(J_dq.block(0,comp*n.sum(),6*(g_pt[nd]+2),n.sum()),Fed_dq.col(comp),X,Ws,n,L[nd]); 
                }
                break;
            case 1:
                K.block(n[0],n[0],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                for (int comp=0; comp < n.sum(); comp++) {
                    Fd_dq.col(comp) = Fd_dq.col(comp) + compute_F(J_dq.block(6*(g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum()),Fed_dq.col(comp),X,Ws,n,L[nd]); 
                }
                break;
            case 2:
                K.block(n[0]+n[1],n[0]+n[1],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(6*(g_pt[nd-2]+2 + g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                for (int comp=0; comp < n.sum(); comp++) {
                    Fd_dq.col(comp) = Fd_dq.col(comp) + compute_F(J_dq.block(6*(g_pt[nd-2]+2 + g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum()),Fed_dq.col(comp),X,Ws,n,L[nd]); 
                }
                break;
            case 3:
                K.block(n[0]+n[1]+n[2],n[0]+n[1]+n[2],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(6*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                for (int comp=0; comp < n.sum(); comp++) {
                    Fd_dq.col(comp) = Fd_dq.col(comp) + compute_F(J_dq.block(6*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2),comp*n.sum(),6*(g_pt[nd]+2),n.sum()),Fed_dq.col(comp),X,Ws,n,L[nd]); 
                }
                break;
        }
    }   
    
        MatrixXd Jhd = MatrixXd::Zero(A_T.rows() + n.sum() - na, A_T.cols());
    
        Jhd.block(0,0,A_T.rows(),A_T.cols()) = A_T;
        Jhd.block(A_T.rows(),na,n.sum() - na,n.sum()-na) = MatrixXd::Identity(n.sum()-na,n.sum()-na);
        

        MatrixXd Jh_dq = MatrixXd::Zero(A_T.rows() + n.sum() - na, A_T.cols()*n.sum());


        for (int comp=0; comp < n.sum(); comp++) {
            Jh_dq.block(0,comp*n.sum(),A_T.rows(),A_T.cols()) = A_dq_T.block(0,comp*n.sum(),na,n.sum());
            Jh_dq.block(A_T.rows(),na+comp*n.sum(),n.sum()-na,n.sum()-na) = MatrixXd::Identity(n.sum()-na,n.sum()-na);
        }

        VectorXd F_load = VectorXd::Zero(n.sum());
        VectorXd F_cup = VectorXd::Zero(n.sum());
        VectorXd F_rings = VectorXd::Zero(n.sum());
        
        MatrixXd F_load_dq = MatrixXd::Zero(n.sum(),n.sum());
        MatrixXd F_cup_dq = MatrixXd::Zero(n.sum(),n.sum());
        MatrixXd F_rings_dq = MatrixXd::Zero(n.sum(),n.sum());

       
        VectorXd grav_cup = VectorXd::Zero(6);
    //   grav_cup << 0,0,0,0,0,-9.81*0.0016;
     grav_cup << 0,0,0,0,0,-9.81*0.0006;

    //cout << "grav_cup" << grav_cup << endl;
   
    VectorXd grav_ring1 = VectorXd::Zero(6);
    //grav_ring1 << 0,0,0,0,0,-9.81*0.0013;
    grav_ring1 << 0,0,0,0,0,-9.81*0.001;

    VectorXd grav_ring2 = VectorXd::Zero(6);
    //grav_ring2 << 0,0,0,0,0,-9.81*0.0013;
    grav_ring2 << 0,0,0,0,0,-9.81*0.001;
   
    VectorXd grav_ring3 = VectorXd::Zero(6);
    //grav_ring3 << 0,0,0,0,0,-9.81*0.0013;
    grav_ring3 << 0,0,0,0,0,-9.81*0.001;
    
        MatrixXd gnewd = MatrixXd::Zero(4,4);
        gnewd=g.bottomRows(4);
        gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
     
        
       F_cup=J.bottomRows(6).transpose()*coadjmap_fcn(gnewd.inverse())*grav_cup;    
    
        gnewd=g.block(4*(g_pt[0]+2 + 4),0,4,4);
        gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
        F_rings = F_rings + J.block(6*(g_pt[0]+2 + 4),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring1;    
     
        gnewd=g.block(4*(g_pt[0]+2 + 6),0,4,4);
        gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
        F_rings = F_rings + J.block(6*(g_pt[0]+2 + 6),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring2;    
     
        gnewd=g.block(4*(g_pt[0]+2 + 9),0,4,4);
        gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
        F_rings = F_rings + J.block(6*(g_pt[0]+2 + 9),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring3;    
     
        Fd = Fd + F_load + F_rings + F_cup;
        
        for (int comp=0; comp < n.sum(); comp++) {
            
            gnewd=g_dq.block(0,comp*4,g_dq.rows(),4).bottomRows(4);
            gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
            
            F_cup_dq.col(comp)=J_dq.block(0, n.sum()*comp,J_dq.rows(),n.sum()).bottomRows(6).transpose()*coadjmap_fcn(gnewd.inverse())*grav_cup;    
    

           gnewd=g_dq.block(4*(g_pt[0]+2 + 4),comp*4,4,4);
           gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
           F_rings_dq.col(comp) = F_rings_dq.col(comp) + J_dq.block(6*(g_pt[0]+2 + 4),n.sum()*comp,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring1;    
        
           gnewd=g_dq.block(4*(g_pt[0]+2 + 6),comp*4,4,4);
           gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
           F_rings_dq.col(comp) = F_rings_dq.col(comp) + J_dq.block(6*(g_pt[0]+2 + 6),n.sum()*comp,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring2;    
        
           gnewd=g_dq.block(4*(g_pt[0]+2 + 9),comp*4,4,4);
           gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
           F_rings_dq.col(comp) = F_rings_dq.col(comp) + J_dq.block(6*(g_pt[0]+2 + 9),n.sum()*comp,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring3;    

           Fd_dq.col(comp) = Fd_dq.col(comp) + F_load_dq.col(comp) + F_rings_dq.col(comp) + F_cup_dq.col(comp);
        }
            
        VectorXd kthetad = (Jhd.transpose()).inverse() * K * qd_old;
        VectorXd gthetad = (Jhd.transpose()).inverse() * Fd;
       
        VectorXd G_theta = kthetad - gthetad;
        VectorXd G_u =G_theta.tail(n.sum()-na);
        MatrixXd H_Guq = MatrixXd::Zero(n.sum()-na,n.sum());

        MatrixXd ktheta_dq = MatrixXd::Zero(n.sum(),n.sum());
        MatrixXd gtheta_dq = MatrixXd::Zero(n.sum(),n.sum());
        MatrixXd G_theta_dq = MatrixXd::Zero(n.sum(),n.sum());
        MatrixXd G_u_dq = MatrixXd::Zero(n.sum()-na,n.sum());

        for (int comp=0; comp < n.sum(); comp++) {
            ktheta_dq.col(comp) = Jh_dq.block(0,comp*n.sum(),n.sum(),n.sum()).transpose().inverse() * K * (qd_old+comp_sel.col(comp)*deltaq);
            gtheta_dq.col(comp) = Jh_dq.block(0,comp*n.sum(),n.sum(),n.sum()).transpose().inverse() * Fd_dq.col(comp);
            G_theta_dq.col(comp) = ktheta_dq.col(comp) - gtheta_dq.col(comp);
            G_u_dq.col(comp) = G_theta_dq.col(comp).tail(n.sum()-na);

            H_Guq.col(comp) = (G_u_dq.col(comp) - G_u) / deltaq;
        }
        
        MatrixXd H_ux = H_Guq*Jhd.inverse();
        MatrixXd H_ua = H_ux.block(0,0,n.sum()-na,na);
        MatrixXd H_uu = H_ux.block(0,na,n.sum()-na,n.sum()-na);

        //tip_ctrl
        //MatrixXd J_v = J.bottomRows(3);
        //MatrixXd R = g.bottomRows(4).block(0,0,3,3);
        //MatrixXd J_v_gl = R*J_v;


        //tipMarker_ctrl
        MatrixXd g_tip = MatrixXd::Zero(4,4);
        g_tip.block(0,0,3,3) = MatrixXd::Identity(3,3);
        VectorXd r_tip = VectorXd::Zero(3);
        //r_tip << 0.01,0,0;
        r_tip << 0.008,0,0;

        g_tip.block(0,3,3,1) = r_tip;
        g_tip.block(3,0,1,3) = MatrixXd::Zero(1,3);
        g_tip(3,3) = 1;
        MatrixXd J_t = adjmap_fcn(ginv(g_tip)) *  J.bottomRows(6);
        MatrixXd J_v = J_t.bottomRows(3);
        MatrixXd R = MatrixXd::Identity(3,3)*g.bottomRows(4).block(0,0,3,3);
        MatrixXd J_v_gl = R*J_v;

        MatrixXd g_tip2 = MatrixXd::Zero(4,4);
        g_tip2.block(0,0,3,3) = MatrixXd::Identity(3,3);
        VectorXd r_tip2 = VectorXd::Zero(3);
        //r_tip2 << 0,0,0.018;
        r_tip2 << 0,0,0.02;
        
        g_tip2.block(0,3,3,1) = r_tip2;
        g_tip2.block(3,0,1,3) = MatrixXd::Zero(1,3);
        g_tip2(3,3) = 1;
        MatrixXd J_t2 = adjmap_fcn(ginv(g_tip2)) *  J.block(6*(g_pt[0]+2 +g_pt[1]+2)-6,0,6,n.sum());
        MatrixXd J_v2 = J_t2.bottomRows(3);
        MatrixXd R2 = MatrixXd::Identity(3,3)*g.block(4*(g_pt[0]+2 +g_pt[1]+2 )-4,0,3,3);
        MatrixXd J_v2_gl = R2*J_v2;

        MatrixXd J_v_gl_all = MatrixXd::Zero(na, n.sum());
    
        //na=4, to be changed according to the desired points
        J_v_gl_all.topRows(3) = J_v_gl;
        J_v_gl_all.bottomRows(1) = J_v2_gl.row(1);

 
        MatrixXd J_hoha = J_v_gl_all*Jhd.inverse();
        MatrixXd J_a = J_hoha.block(0,0,na,na);
        MatrixXd J_u = J_hoha.block(0,na,na,n.sum()-na);


        MatrixXd temp_Peps = MatrixXd::Zero(n.sum(),na);
        temp_Peps.topRows(na) = MatrixXd::Identity(na,na);
        temp_Peps.bottomRows(n.sum()-na) = -H_uu.inverse() * H_ua;

        MatrixXd J_aINV = MatrixXd::Zero(na,na); 
        double lambda = 0.0001;
        if (J_a.determinant() < 1e-10){
            J_aINV = J_a.transpose()*(J_a*J_a.transpose()+lambda*lambda*MatrixXd::Identity(na,na) ).inverse();
        }else{
            J_aINV = J_a.inverse(); 
        }

        MatrixXd temp_PZ = MatrixXd::Zero(n.sum(),n.sum()-na);
        temp_PZ.topRows(na) = -J_aINV* J_u;
        temp_PZ.bottomRows(n.sum()-na) = MatrixXd::Identity(n.sum()-na,n.sum()-na);

        MatrixXd P_eps = temp_Peps * (MatrixXd::Identity(na,na) - J_aINV*J_u*H_uu.inverse()*H_ua).inverse();
        MatrixXd P_Z = temp_PZ * (MatrixXd::Identity(n.sum()-na,n.sum()-na) - H_uu.inverse()*H_ua*J_aINV*J_u).inverse();

        //tipMarker_ctrl
        VectorXd x = xe;

        VectorXd thetadot = P_eps* J_aINV* Gamma_Z * (xd - x) + P_Z * (-H_uu.inverse()*Gamma_eps*G_u);

        VectorXd qdot = Jhd.inverse()*thetadot;
        qd = qd_old + tsample * qdot;
     

    return qd;

}



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
){
    VectorXd qd = VectorXd::Zero(n.sum());

    MatrixXd Gamma_Z = 0.1*MatrixXd::Identity(na, na); //2
    MatrixXd Gamma_eps = 0.1*MatrixXd::Identity(n.sum()-na,n.sum()- na);  //2

    //if (change_g)
    //{
    //    Gamma_Z = 0.1*MatrixXd::Identity(na, na);
    //    Gamma_eps = 0.1*MatrixXd::Identity(n.sum()-na,n.sum()- na);
    //}
  
   // Initialize structures to hold the values
   VectorXd X,Ws,h;
   VectorXd r, A, Jx, Jy, Jz;
   MatrixXd Bq;
   MatrixXd d_first,dtilde;
   VectorXd xi;
   MatrixXd xi_hat;
   VectorXd tc;
   MatrixXd t;
   MatrixXd Btau;
   MatrixXd friction;
   MatrixXd Sigma, Mbar, Fed;
   MatrixXd tc_normINV;

    VectorXd y = VectorXd::Zero(na);

   double frc = 0.0;


    MatrixXd dA1_dq = MatrixXd::Zero(n[0],na*n.sum());
    MatrixXd dA2_dq = MatrixXd::Zero(n[1],na*n.sum());
    MatrixXd dA3_dq = MatrixXd::Zero(n[2],na*n.sum());
    MatrixXd dA4_dq = MatrixXd::Zero(n[3],na*n.sum());

   Matrix3d Rot_in;
   Rot_in << 0, 0, 1, 0, 1, 0, -1, 0 ,0;
   MatrixXd g_in = MatrixXd::Zero(4,4);
   g_in.block(0,0,3,3) = Rot_in;
   g_in.block(0,3,3,1) = MatrixXd::Zero(3,1);
   g_in.block(3,0,1,3) = MatrixXd::Zero(1,3);
   g_in(3,3) = 1;
   MatrixXd g0=MatrixXd::Identity(4, 4)*g_in;
   MatrixXd J0=MatrixXd::Zero(6, n.sum());

   VectorXi g_pt(4);
   g_pt << 5, 8, 5, 5;

   MatrixXd g = MatrixXd::Zero(4 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), 4);
   MatrixXd phi = MatrixXd::Zero(g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2, na);
   MatrixXd J = MatrixXd::Zero(6 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), n.sum());
   MatrixXd d = MatrixXd::Zero(3 * (g_pt[0] + 2 + g_pt[1] + 2 + g_pt[2] + 2 + g_pt[3] + 2), na);

   MatrixXd A_T(na, n.sum());
   A_T.setZero();


   MatrixXd K(n.sum(), n.sum());
   K.setZero();
   VectorXd Fd(n.sum());
   Fd.setZero();
   MatrixXd Fd_dq(n.sum(),n.sum());
   Fd_dq.setZero();

   VectorXd gravity(6);
   gravity << 0, 0, 0, 0, 0, -9.81;

   MatrixXd comp_sel = MatrixXd::Identity(n.sum(),n.sum());
    

   for (int nd = 0; nd < n_divisions; nd++) {
    X.resize(g_pt[nd]+2);
    Ws.resize(g_pt[nd]+2);
    h.resize(g_pt[nd]+1);
    
    
    MatrixXd res = MatrixXd::Zero(g_pt[nd] + 2, 2);

    res = GaussQuadrature(g_pt[nd]);
    X = res.col(0);
    Ws = res.col(1);
    
    // Scale X by L
    for (int i = 0; i < g_pt[nd]+2; i++) {
        X[i] = X[i] * L[nd];

    }


    for (int i = 0; i < g_pt[nd]+1; i++) {
        h[i] = X[i+1] - X[i];
    }

    switch(nd) {
        case 0:
            g.block(0,0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(0,0,4*(g_pt[nd]+2),4).bottomRows(4);                
            
            J.block(0,0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(0,0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
    
            break;


        case 1:
            g.block(4*(g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(4*(g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);    
            
            J.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);


            break;
        case 2:
            g.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(4*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);         

            J.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(6*(g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
            
        

            break;
        case 3:

            g.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4)=direct_kin(g0,qd_old,h,g_pt[nd],L[nd],nd,n);
            g0=g.block(4*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4).bottomRows(4);     
   
            
            J.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum())=compute_J(J0,qd_old,h,g_pt[nd],L[nd],nd,n); //from up to down we have the different Jac of the body points
            J0=J.block(6*(g_pt[nd-3]+2+g_pt[nd-2]+2+g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()).bottomRows(6);
        
            break;
    }
    for (int ii = 0; ii < g_pt[nd] + 2; ii++) {
           
        switch (nd) {
            case 0:
                d.block(ii*3,0,3,na) << 0,0,0,0,
                (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
         
            // left-> right:[x;y;z] for each actuator,
            // up->down: i repeat for each point of the section
            break;
            
            case 1:
                d.block(3*(g_pt[nd-1] + 2)+ii*3,0,3,na) << 0,0,0,0,
                (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),(dist(1,nd)+m(0,nd)*X[ii])*sin(alpha[1]),(dist(1,nd)+m(0,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
             break;
            case 2:
                d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,0,3,na) << 0,0,0,0,
                (dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[0]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[1]),(dist(1,nd)+m(1,nd)*X[ii])*cos(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*cos(alpha[3]),
                (dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[0]),-(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[1]),-(dist(1,nd)+m(1,nd)*X[ii])*sin(alpha[2]),(dist(0,nd)+m(0,nd)*X[ii])*sin(alpha[3]);
            
                break;
            case 3:
                d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,0,3,na) = MatrixXd::Zero(3,na);

                break;
        }
    }     
}  
for (int nd = 0; nd < n_divisions; nd++) {
    // Gauss quadrature and scaling
    X.resize(g_pt[nd]+2);
    Ws.resize(g_pt[nd]+2);
    r.resize(g_pt[nd]+2);
    A.resize(g_pt[nd]+2);
    Jx.resize(g_pt[nd]+2);
    Jy.resize(g_pt[nd]+2);
    Jz.resize(g_pt[nd]+2);
    h.resize(g_pt[nd]+1);
    
    
    MatrixXd res = MatrixXd::Zero(g_pt[nd] + 2, 2);

    res = GaussQuadrature(g_pt[nd]);
    X = res.col(0);
    Ws = res.col(1);
    
    // Scale X by L
    for (int i = 0; i < g_pt[nd]+2; i++) {
        X[i] = X[i] * L[nd];

        r[i]= Rb[nd] + (Rt[nd] - Rb[nd]) * (X[i]/ L[nd]); //X not normalized
        A[i]=M_PI*pow(r[i],2);
        Jx[i]=M_PI/2*pow(r[i],4);
        Jy[i]=M_PI/4*pow(r[i],4);
        Jz[i]=Jy[i];
    }


    for (int i = 0; i < g_pt[nd]+1; i++) {
        h[i] = X[i+1] - X[i];
    }

    t.resize(3*(g_pt[nd]+2),na);
    Bq.resize(6*(g_pt[nd] + 2),n[nd]);
    d_first.resize(3*(g_pt[nd] + 2),na);
    dtilde.resize(3*(g_pt[nd] + 2),3*na);
    xi.resize(6*(g_pt[nd] + 2));
    Btau.resize(6*(g_pt[nd] + 2),na);
    friction.resize(1*(g_pt[nd] + 2),na);
    tc_normINV.resize(1*(g_pt[nd] + 2),na);
    
    Sigma.resize(6*(g_pt[nd] + 2),6*(g_pt[nd] + 2));
    Mbar.resize(6*(g_pt[nd] + 2),6*(g_pt[nd] + 2));
    Fed.resize(6*(g_pt[nd] + 2),1);

    t.setZero();
    Bq.setZero();
    d_first.setZero();
    dtilde.setZero();
    xi.setZero();
    Btau.setZero();
    Sigma.setZero();
    Mbar.setZero();
    Fed.setZero();
    friction.setZero();
    tc_normINV.setZero();

    // Loop for each Gauss point and division
    for (int ii = 0; ii < g_pt[nd] + 2; ii++) {
        VectorXd sigma_el(6),mbar_el(6);
        sigma_el << G*Jx[ii],E*Jy[ii],E*Jz[ii],E*A[ii],G*A[ii],G*A[ii];
        Sigma.block(6*ii,6*ii,6,6)=sigma_el.asDiagonal();
        mbar_el << Jx[ii],Jy[ii],Jz[ii],A[ii],A[ii],A[ii];
        Mbar.block(6*ii,6*ii,6,6)=rho*mbar_el.asDiagonal();

        MatrixXd gnew = MatrixXd::Zero(4,4);

        switch (nd) {
            case 0:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq1(X[ii], L[nd]);
                
                gnew=g.block(4*ii,0,4,4);


                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),m(1,nd)*sin(alpha[1]),m(1,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
        
          

                break;
            case 1:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq2(X[ii], L[nd]);
                
                gnew=g.block(4*(g_pt[nd-1]+2 + ii),0,4,4);

                            
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),m(0,nd)*sin(alpha[1]),m(0,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
            
                break;
            case 2:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq3(X[ii], L[nd]);

                gnew=g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,4,4);

                 
                d_first.block(ii*3,0,3,na) << 0,0,0,0,
                m(0,nd)*cos(alpha[0]),m(1,nd)*cos(alpha[1]),m(1,nd)*cos(alpha[2]),m(0,nd)*cos(alpha[3]),
                m(0,nd)*sin(alpha[0]),-m(1,nd)*sin(alpha[1]),-m(1,nd)*sin(alpha[2]),m(0,nd)*sin(alpha[3]);
           
          
                break;
            case 3:
                Bq.block(6*ii,0,6,n[nd]) = compute_Bq4(X[ii], L[nd]);
                
                gnew=g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,4,4);

                d_first.block(ii*3,0,3,na) = MatrixXd::Zero(3,na);
                   
                break;
        }

        gnew.block(0,3,3,1)=VectorXd::Zero(3,1);   

        Fed.block(6*ii,0,6,1)=Mbar.block(6*ii,6*ii,6,6)*coadjmap_fcn(gnew.inverse())*gravity;    

        for (int kk = 0; kk < na; kk++) {

            VectorXd d_homg = VectorXd::Zero(4);
            VectorXd d_first_homg = VectorXd::Zero(4);

            switch(nd) {
                case 0:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(ii*3,kk,3,1)); 
                    // left-> right:3x3 dtilde for each actuator,
                    // up->down: i repeat for each point of the section
                
                    //cout <<"dtilde.block(ii*3,kk*3,3,3):"<<endl << dtilde.block(ii*3,kk*3,3,3) << endl;
                
                    d_homg << d.block(ii*3,kk,3,1),1;                 
                    break;
                case 1:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                    d_homg << d.block(3*(g_pt[nd-1] + 2)+ii*3,kk,3,1),1;                 
                    break;
                case 2:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                    d_homg << d.block(3*(g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1),1;              
                    break;
                case 3:
                    dtilde.block(ii*3,kk*3,3,3) = skewMat_fcn(d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1)); 
                    d_homg << d.block(3*(g_pt[nd-3] + 2 + g_pt[nd-2] + 2 + g_pt[nd-1] + 2)+ii*3,kk,3,1),1;
                    break;
            }
            xi.segment(6*ii,6)=compute_xi(qd_old,X[ii],nd,n,L[nd]);    
            xi_hat=HatFromVec_fcn(xi.segment(6*ii,6));

            d_first_homg << d_first.block(ii*3,kk,3,1),1;
        
            tc = (xi_hat*d_homg + d_first_homg);
            t.block(ii*3,kk,3,1) = tc.head(3)/tc.head(3).norm();     
            tc_normINV(ii,kk) = 1/tc.head(3).norm(); 

            double  th,dot_prod;
            Vector3d p,pn,f,fn;

            //friction model
            if (nd == (n_divisions - 1)){              
                th =0; 
            }
            else{
                if(ii == 0){
                    switch(nd){
                        case 0:    
                            p << 1,0,0;
                            pn = g.block(4*(ii+1),3,3,1) - g.block(4*(ii),3,3,1) + g.block(4*(ii+1),0,3,3) * d.block(3*(ii+1),kk,3,1) - g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);
                      

                            break;
                        case 1:
                            p =  g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-1]),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-1]),0,3,3) * d.block(3*(g_pt[nd-1]),kk,3,1);
                            pn = g.block(4*(g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                            
                         
                             break;
                        case 2:
                            p =  g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1)   - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]),kk,3,1);
                            pn = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                           
                            
                     
                            break;                            
                    }

                    f=p/p.norm();
                    fn=pn/pn.norm();
                    dot_prod = f.dot(fn);
                    th = acos(max(-1.0, min(dot_prod, 1.0)));


                }
                else if (ii == (g_pt[nd] + 1)){
                    
                    switch(nd){
                        case 0:
                          
                        

                            p = g.block(4*(ii),3,3,1) - g.block(4*(ii-1),3,3,1) + g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1) - g.block(4*(ii-1),0,3,3) * d.block(3*(ii-1),kk,3,1);
                            pn = g.block(4*(g_pt[nd]+2 +1),3,3,1) - g.block(4*(ii),3,3,1) + g.block(4*(g_pt[nd]+2 +1),0,3,3) * d.block(3*(g_pt[nd]+2 +1),kk,3,1) - g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);
                            
                            f=p/p.norm();
                            fn=pn/pn.norm();
                            dot_prod = f.dot(fn);
                            th = acos(max(-1.0, min(dot_prod, 1.0)));
                        
                            break;
                        case 1:

                        

                            p = g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                            pn = g.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + g_pt[nd]+2 +1),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                    
                            f=p/p.norm();
                            fn=pn/pn.norm();
                            dot_prod = f.dot(fn);
                            th = acos(max(-1.0, min(dot_prod, 1.0)));
                            break;    
                        case 2:

                            p = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                        
                            f=p/p.norm();
                            fn=VectorXd::Zero(3);
                            dot_prod = f.dot(fn);
                            th = acos(max(-1.0, min(dot_prod, 1.0)));   
                            break;
                    }
                }
                else{
                    switch(nd){
                        case 0:

                            p = g.block(4*(ii),3,3,1) - g.block(4*(ii-1),3,3,1) + g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1) - g.block(4*(ii-1),0,3,3) * d.block(3*(ii-1),kk,3,1);
                            pn = g.block(4*(ii+1),3,3,1) - g.block(4*(ii),3,3,1) + g.block(4*(ii+1),0,3,3) * d.block(3*(ii+1),kk,3,1) - g.block(4*(ii),0,3,3) * d.block(3*(ii),kk,3,1);

                            break;
                        case 1:
                           
                            p = g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii-1),kk,3,1);
                            pn = g.block(4*(g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-1]+2 + ii),kk,3,1);
                            
                            break;    
                        case 2:
                          
                           

                            p = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii-1),kk,3,1);                           
                            pn = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),3,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),3,3,1) + g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1),kk,3,1) - g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),0,3,3) * d.block(3*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii),kk,3,1);                           
                            
                            break;
                    }
                   
                
                    f=p/p.norm();
                    fn=pn/pn.norm();
                    dot_prod = f.dot(fn);
                    th = acos(max(-1.0, min(dot_prod, 1.0)));
                }
            }    


            //phi
            switch(nd) {
                case 0:
                    phi(ii,kk)=th;
                    
                    break;
                case 1:
                    phi(g_pt[nd-1]+2 + ii,kk)=th;
                    
                    break;
                case 2:
                    phi(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk)=th;
                    
                    
                    break;
                case 3:
                    phi(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii,kk)=th;
                   
                    break;
            }


            switch(nd) {
                case 0:
                    frc = phi.col(kk).topRows(ii+1).sum();
                   

                    break;
                case 1:
                    frc = phi.col(kk).topRows(g_pt[nd-1]+2 + ii+1).sum()-phi(g_pt[nd-1]+2,kk);  
                
                    break;
                case 2:
                    frc = phi.col(kk).topRows(g_pt[nd-2]+2 + g_pt[nd-1]+2 + ii+1).sum()-phi(g_pt[nd-2]+2,kk)- phi(g_pt[nd-2]+2 + g_pt[nd-1]+2,kk);
   
                    break;
                case 3:
                    frc = 0;

                    break;
            }

            // easier alternative version  
            //  switch(nd) {
            //    case 0:
            //        frc = X[ii]/(L[0]+L[1]+L[2]);
            //        break;
            //    case 1:
            //        frc = (L[0]+X[ii])/(L[0]+L[1]+L[2]);
            //        break;
            //    case 2:
            //        frc = (L[0]+L[1]+X[ii])/(L[0]+L[1]+L[2]);
            //        break;
            //    case 3:
            //        frc = 0;
            //        break;
            //}
 

            if (nd != (n_divisions - 1)) {
                // Calcolo dell'attrito
                friction(ii,kk) = exp(-mu[kk] * frc);

                Btau.block(6*ii, kk, 3, 1) = dtilde.block(3*ii,3*kk,3,3)* t.block(ii*3,kk,3,1)*friction(ii,kk);
                Btau.block(3 + 6*ii, kk, 3, 1) = t.block(ii*3,kk,3,1)*friction(ii,kk);
            } else { 
                friction(ii,kk) = exp(-mu[kk] * frc);
                Btau.block(6*ii, kk, 6, 1).setZero();

            }
        }
    }
        for (int kk = 0; kk < na; kk++) {
            y [kk] = y [kk] + compute_PassiveOutput(Btau.col(kk),xi,d_first.col(kk),X,Ws,L[nd]);

               switch (nd) {
                case 0:
                    A_T.block(kk,0,1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);              
                
                  break;
                case 1:
                    A_T.block(kk,n[0],1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);
                                 
                   
                    break;
                case 2:
                    A_T.block(kk,n[0]+n[1],1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);
                 
                    break;
                case 3:
                    A_T.block(kk,n[0]+n[1]+n[2],1,n[nd]) = compute_AT(Btau.col(kk),Bq,X,Ws,L[nd]);
              
                    break;
            }
        }

      
        MatrixXd dA_dq_cell_niXniXna = compute_dAdq(t,dtilde,tc_normINV,friction,Bq,nd,X,Ws,n,na,L[nd]);


        int count = 0;
        MatrixXd temp = MatrixXd::Zero(n[nd],n[nd]*na);
        for (int nn = 0; nn<n[nd]; nn++){
            for (int cc = nn; cc<na*n[nd]; cc=cc+n[nd]){
                temp.block(0,count,n[nd],1)=dA_dq_cell_niXniXna.block(0,cc,n[nd],1);
                count++;
            }
        }
        switch (nd) {
            case 0:
                K.block(0,0,n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(0,0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                // Fd_dq.block(0,0,n[0],n.sum()) = L[0]*compute_dFdq(qd_old,g.block(0,0,4*(g_pt[nd]+2),4),X,Ws,h,nd,n,L[nd],Mbar,g_pt);     
                Fd_dq.block(0,0,n[0],n.sum()) = L[0]*compute_dFdq(qd_old,g,X,Ws,h,nd,n,L[nd],Mbar,g_pt);     

                dA1_dq.block(0,0,n[0],na*n[0]) = temp;
           
                break;
            case 1:
                K.block(n[0],n[0],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(6*(g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                
                // Fd_dq.block(n[0],0,n[1],n.sum()) = L[1]*compute_dFdq(qd_old,g.block(4*(g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4),X,Ws,h,nd,n,L[nd],Mbar,g_pt);
                Fd_dq.block(n[0],0,n[1],n.sum()) = L[1]*compute_dFdq(qd_old,g,X,Ws,h,nd,n,L[nd],Mbar,g_pt);

                dA2_dq.block(0,na*n[0],n[1],na*n[1]) = temp;
                 break;
            case 2:
                K.block(n[0]+n[1],n[0]+n[1],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(6*(g_pt[nd-2]+2 + g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]); 
                // Fd_dq.block(n[0]+n[1],0,n[2],n.sum()) = L[2]*compute_dFdq(qd_old,g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4),X,Ws,h,nd,n,L[nd],Mbar,g_pt);
                Fd_dq.block(n[0]+n[1],0,n[2],n.sum()) = L[2]*compute_dFdq(qd_old,g,X,Ws,h,nd,n,L[nd],Mbar,g_pt);
                dA3_dq.block(0,na*(n[0]+n[1]),n[2],na*n[2]) = temp;

                break;
            case 3:
                K.block(n[0]+n[1]+n[2],n[0]+n[1]+n[2],n[nd],n[nd])=compute_K(Sigma,nd,X,Ws,n,Bq,L[nd]);
                Fd = Fd + compute_F(J.block(6*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2),0,6*(g_pt[nd]+2),n.sum()),Fed,X,Ws,n,L[nd]);
                // Fd_dq.block(n[0]+n[1]+n[2],0,n[3],n.sum()) = L[3]*compute_dFdq(qd_old,g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2),0,4*(g_pt[nd]+2),4),X,Ws,h,nd,n,L[nd],Mbar,g_pt); 
                Fd_dq.block(n[0]+n[1]+n[2],0,n[3],n.sum()) = L[3]*compute_dFdq(qd_old,g,X,Ws,h,nd,n,L[nd],Mbar,g_pt); 

                dA4_dq.block(0,na*(n[0]+n[1]+n[2]),n[3],na*n[3]) = temp;

               break;
        }
    }   
        MatrixXd Jhd = MatrixXd::Zero(A_T.rows() + n.sum() - na, A_T.cols());
        MatrixXd JhT = MatrixXd::Zero(n.sum(),n.sum());
        JhT.block(0,0,n.sum(), na) = A_T.transpose();
        JhT.block(na,na,n.sum() - na,n.sum()-na) = MatrixXd::Identity(n.sum()-na,n.sum()-na);
        
        
        Jhd.block(0,0,A_T.rows(),A_T.cols()) = A_T;
        Jhd.block(A_T.rows(),na,n.sum() - na,n.sum()-na) = MatrixXd::Identity(n.sum()-na,n.sum()-na);
        


    //     VectorXd F_load = VectorXd::Zero(n.sum());
    //     VectorXd F_cup = VectorXd::Zero(n.sum());
    //     VectorXd F_rings = VectorXd::Zero(n.sum());

       
    //     VectorXd grav_cup = VectorXd::Zero(6);
    // //   grav_cup << 0,0,0,0,0,-9.81*0.0016;
    //  grav_cup << 0,0,0,0,0,-9.81*0.0006;

    // //cout << "grav_cup" << grav_cup << endl;
   
    // VectorXd grav_ring1 = VectorXd::Zero(6);
    // //grav_ring1 << 0,0,0,0,0,-9.81*0.0013;
    // grav_ring1 << 0,0,0,0,0,-9.81*0.001;

    // VectorXd grav_ring2 = VectorXd::Zero(6);
    // //grav_ring2 << 0,0,0,0,0,-9.81*0.0013;
    // grav_ring2 << 0,0,0,0,0,-9.81*0.001;
   
    // VectorXd grav_ring3 = VectorXd::Zero(6);
    // //grav_ring3 << 0,0,0,0,0,-9.81*0.0013;
    // grav_ring3 << 0,0,0,0,0,-9.81*0.001;
    
    //     MatrixXd gnewd = MatrixXd::Zero(4,4);
    //     gnewd=g.bottomRows(4);
    //     gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
     
        
    //    F_cup=J.bottomRows(6).transpose()*coadjmap_fcn(gnewd.inverse())*grav_cup;    
    
    //     gnewd=g.block(4*(g_pt[0]+2 + 4),0,4,4);
    //     gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
    //     F_rings = F_rings + J.block(6*(g_pt[0]+2 + 4),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring1;    
     
    //     gnewd=g.block(4*(g_pt[0]+2 + 6),0,4,4);
    //     gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
    //     F_rings = F_rings + J.block(6*(g_pt[0]+2 + 6),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring2;    
     
    //     gnewd=g.block(4*(g_pt[0]+2 + 9),0,4,4);
    //     gnewd.block(0,3,3,1)=VectorXd::Zero(3,1);
    //     F_rings = F_rings + J.block(6*(g_pt[0]+2 + 9),0,6,n.sum()).transpose()*coadjmap_fcn(gnewd.inverse())*grav_ring3;    
     
        // Fd = Fd + F_load + F_rings + F_cup;
            

        VectorXd kthetad = (Jhd.transpose()).inverse() * K * qd_old;
        VectorXd gthetad = (Jhd.transpose()).inverse() * Fd;
       
        VectorXd G_theta = kthetad - gthetad;
        VectorXd G_u = G_theta.tail(n.sum()-na);

       
        MatrixXd dA_dq_SnixnaxSni = MatrixXd::Zero(n.sum(),na*n.sum());
        MatrixXd dJhT_dq = MatrixXd::Zero(n.sum(),n.sum()*n.sum()); 
        MatrixXd dJhmT_dq = MatrixXd::Zero(n.sum(),n.sum()*n.sum());
        MatrixXd dJhmT_dq_KqF = MatrixXd::Zero(n.sum(),n.sum());
        
        for (int nn=0; nn < n.sum(); nn++) {
            dA_dq_SnixnaxSni.block(0,na*(nn),n[0],na) = dA1_dq.block(0,na*(nn),n[0],na); 
            dA_dq_SnixnaxSni.block(n[0],na*(nn),n[1],na) = dA2_dq.block(0,na*(nn),n[1],na);
            dA_dq_SnixnaxSni.block(n[0]+n[1],na*(nn),n[2],na) = dA3_dq.block(0,na*(nn),n[2],na);
            dA_dq_SnixnaxSni.block(n[0]+n[1]+n[2],na*(nn),n[3],na) = dA4_dq.block(0,na*(nn),n[3],na);
        
            dJhT_dq.block(0,nn*n.sum(),n.sum(),na) = dA_dq_SnixnaxSni.block(0,nn*na,n.sum(),na);
            dJhT_dq.block(0,na+nn*n.sum(),n.sum(),n.sum()-na) = MatrixXd::Zero(n.sum(),n.sum()-na);

            dJhmT_dq.block(0,nn*n.sum(),n.sum(),n.sum()) = - JhT.inverse() * dJhT_dq.block(0,nn*n.sum(),n.sum(),n.sum()) * JhT.inverse();
            dJhmT_dq_KqF.block(0,nn,n.sum(),1) = dJhmT_dq.block(0,nn*n.sum(),n.sum(),n.sum()) * (K*qd_old - Fd);
        
        }

        MatrixXd H_Gthq = dJhmT_dq_KqF + JhT.inverse() * (K-Fd_dq);
        
        MatrixXd H_Guq = H_Gthq.block(na,0,n.sum()-na,n.sum());
        
        MatrixXd H_ux = H_Guq*Jhd.inverse();
        MatrixXd H_ua = H_ux.block(0,0,n.sum()-na,na);
        MatrixXd H_uu = H_ux.block(0,na,n.sum()-na,n.sum()-na);

        //tip_ctrl
        //MatrixXd J_v = J.bottomRows(3);
        //MatrixXd R = g.bottomRows(4).block(0,0,3,3);
        //MatrixXd J_v_gl = R*J_v;


        //tipMarker_ctrl
        MatrixXd g_tip = MatrixXd::Zero(4,4);
        g_tip.block(0,0,3,3) = MatrixXd::Identity(3,3);
        VectorXd r_tip = VectorXd::Zero(3);
        //r_tip << 0.01,0,0;
        r_tip << 0.008,0,0;

        g_tip.block(0,3,3,1) = r_tip;
        g_tip.block(3,0,1,3) = MatrixXd::Zero(1,3);
        g_tip(3,3) = 1;
        MatrixXd J_t = adjmap_fcn(ginv(g_tip)) *  J.bottomRows(6);
        MatrixXd J_v = J_t.bottomRows(3);
        MatrixXd R = MatrixXd::Identity(3,3)*g.bottomRows(4).block(0,0,3,3);
        MatrixXd J_v_gl = R*J_v;

        MatrixXd g_tip2 = MatrixXd::Zero(4,4);
        g_tip2.block(0,0,3,3) = MatrixXd::Identity(3,3);
        VectorXd r_tip2 = VectorXd::Zero(3);
        //r_tip2 << 0,0,0.018;
        r_tip2 << 0,0,0.02;
        
        g_tip2.block(0,3,3,1) = r_tip2;
        g_tip2.block(3,0,1,3) = MatrixXd::Zero(1,3);
        g_tip2(3,3) = 1;
        MatrixXd J_t2 = adjmap_fcn(ginv(g_tip2)) *  J.block(6*(g_pt[0]+2 +g_pt[1]+2)-6,0,6,n.sum());
        MatrixXd J_v2 = J_t2.bottomRows(3);
        MatrixXd R2 = MatrixXd::Identity(3,3)*g.block(4*(g_pt[0]+2 +g_pt[1]+2 )-4,0,3,3);
        MatrixXd J_v2_gl = R2*J_v2;

        MatrixXd J_v_gl_all = MatrixXd::Zero(na, n.sum());
    
        //na=4, to be changed according to the desired points
        J_v_gl_all.topRows(3) = J_v_gl;
        J_v_gl_all.bottomRows(1) = J_v2_gl.row(1);

 
        MatrixXd J_hoha = J_v_gl_all*Jhd.inverse();
        MatrixXd J_a = J_hoha.block(0,0,na,na);
        MatrixXd J_u = J_hoha.block(0,na,na,n.sum()-na);


        MatrixXd temp_Peps = MatrixXd::Zero(n.sum(),na);
        temp_Peps.topRows(na) = MatrixXd::Identity(na,na);
        temp_Peps.bottomRows(n.sum()-na) = -H_uu.inverse() * H_ua;

        MatrixXd J_aINV = MatrixXd::Zero(na,na); 
        double lambda = 0.0001;
        if (abs(J_a.determinant()) < 1e-10){
            J_aINV = J_a.transpose()*(J_a*J_a.transpose()+lambda*lambda*MatrixXd::Identity(na,na) ).inverse();
        }else{
            J_aINV = J_a.inverse(); 
        }

        MatrixXd temp_PZ = MatrixXd::Zero(n.sum(),n.sum()-na);
        temp_PZ.topRows(na) = -J_aINV* J_u;
        temp_PZ.bottomRows(n.sum()-na) = MatrixXd::Identity(n.sum()-na,n.sum()-na);

        MatrixXd P_eps = temp_Peps * (MatrixXd::Identity(na,na) - J_aINV*J_u*H_uu.inverse()*H_ua).inverse();
        MatrixXd P_Z = temp_PZ * (MatrixXd::Identity(n.sum()-na,n.sum()-na) - H_uu.inverse()*H_ua*J_aINV*J_u).inverse();

        //tipMarker_ctrl
        VectorXd x = xe;

        VectorXd thetadot = P_eps* J_aINV* Gamma_Z * (xd - x) + P_Z * (-H_uu.inverse()*Gamma_eps*G_u);

        VectorXd qdot = Jhd.inverse()*thetadot;
        qd = qd_old + tsample * qdot;
     

    return qd;

}

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
    double L){

    MatrixXd dA_dq = MatrixXd :: Zero(n[nd],n[nd]*na);
    MatrixXd dAfun(X.size()*n[nd],n[nd]*na);  
    dAfun.setZero();
    int numSections = X.size();

    for (int ii = 0; ii < numSections; ii++) {
        MatrixXd dBtau_dq_section = MatrixXd::Zero(6, n[nd]*na); 
        MatrixXd Bq_section = Bq.block(6 * ii, 0, 6, Bq.cols());  
        MatrixXd temp = MatrixXd::Zero(6,6);
        
        for (int kk = 0; kk < na; kk++){
            MatrixXd ttilde = skewMat_fcn(t.block(ii*3,kk,3,1));
            temp.block(0,0,3,3) = dtilde.block(3*ii,3*kk,3,3) * ttilde* ttilde *dtilde.block(3*ii,3*kk,3,3) ;
            temp.block(0,3,3,3) = -dtilde.block(3*ii,3*kk,3,3) *ttilde* ttilde;
            temp.block(3,0,3,3) =  ttilde* ttilde*dtilde.block(3*ii,3*kk,3,3);
            temp.block(3,3,3,3) = -ttilde* ttilde;
            dBtau_dq_section.block(0,kk*n[nd],6,n[nd]) = tc_normINV(ii,kk)*friction(ii,kk) * temp * Bq_section;
        }

        dAfun.block(n[nd]*ii, 0, n[nd],n[nd]*na) = Bq_section.transpose() * dBtau_dq_section;
    }

    dA_dq = compute_intGQ(n[nd], dAfun, X, Ws, L);
    return dA_dq;

    }

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
    const VectorXi& g_pt){

        MatrixXd dFdq = MatrixXd::Zero(n[nd],n.sum());
        VectorXd grav6(6);
        grav6 << 0,0,0,0,0,-9.81;

        double eps_min = 1e-4;
        double eps_max = 1e8;
        for (int aa = 0; aa < g_pt[nd]+2-1; aa++) {
            VectorXd FC = MatrixXd::Zero(6,1);
            for (int kk = aa+1;kk < g_pt[nd]+2; kk++){
                MatrixXd grel = MatrixXd::Zero(4,4);
                MatrixXd gtemp = MatrixXd::Zero(4,4);
                switch (nd){
                    case 0:
                        grel = ginv(g.block(4*aa,0,4,4))*g.block(4*kk,0,4,4);
                        gtemp = g.block(4*kk,0,4,4);
                    break;
                    case 1:
                        grel = ginv(g.block(4*(g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-1]+2 + kk),0,4,4);
                        gtemp = g.block(4*(g_pt[nd-1]+2 + kk),0,4,4);
                    break;
                    case 2:
                        grel = ginv(g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                        gtemp = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);

                    break;
                    case 3:
                        grel = ginv(g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                        gtemp = g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                    break;
                }
                MatrixXd Fcal = Ws[kk]*Mbar.block(6*kk,6*kk,6,6)*adjmap_fcn(ginv(gtemp))*grav6;
                FC = FC + Fcal;
            }
            MatrixXd SUM_dAddq_Fcal = MatrixXd::Zero(6, n.sum());
            MatrixXd SUM_Ad_dFcaldq = MatrixXd::Zero(6, n.sum());

            double X1 = X[aa] + h[aa]/2.0 - sqrt(3.0)*h[aa]/6.0;
            double X2 = X[aa] + h[aa]/2.0 + sqrt(3.0)*h[aa]/6.0;

            VectorXd xi1 = compute_xi(qd_old, X1, nd, n, L); 
            VectorXd xi2 = compute_xi(qd_old, X2, nd, n, L);
            MatrixXd adjop_xi1 = adjop_fcn(xi1);
            VectorXd omega = h[aa]/2.0 * (xi1 + xi2) + sqrt(3.0)*h[aa]*h[aa]/12.0 * (adjop_xi1 * xi2);

            double theta= omega.head(3).norm();
            MatrixXd adjop_omega = adjop_fcn(omega);

            MatrixXd T = compute_T(theta, adjop_omega);
            MatrixXd Z = compute_Z(h[aa], xi1, xi2, nd, X1, X2, L); 

            MatrixXd S = T*Z;

            MatrixXd SUM1 = MatrixXd::Zero(6, n[nd]);
            MatrixXd SUM2 = MatrixXd::Zero(6, n[nd]);    

            for (int rr = 1; rr < 4; rr++){
                double fth = 0.0;
                double fdth_fracth = 0.0;
                double th = theta;
                double c = cos(th);
                double s = sin(th);
                switch (rr) {
                    case 1:
                        if (theta < eps_min) {
                            fth = 0.5;
                            fdth_fracth = 0.0;
                        } else if (theta > eps_max) {
                            fth = 0.0;
                            fdth_fracth = 0.0;
                        } else {
                           fth = (4.0 - 4.0*c - th*s) / (2.0 * th * th);
                           fdth_fracth = (-8.0 + (8.0 - th*th)*c + 5.0*th*s) / (2.0 * pow(th,4));
                        }
                       
                        break;
                    case 2:
                        if (theta < eps_min) {
                            fth = 1.0/6.0;
                            fdth_fracth = 0.0;
                        } else if (theta > eps_max) {
                            fth = 0.0;
                            fdth_fracth = 0.0;
                        } else {
                        fth = (4.0*th - 5.0*s + th*c) / (2.0 * pow(th,3));
                        fdth_fracth = (-8.0*th + (15.0 - th*th)*s - 7.0*th*c) / (2.0 * pow(th,5));
                        }
                        break;
                    case 3:
                         if (theta < eps_min) {
                            fth = 1.0/24.0;
                            fdth_fracth = -1.0/180.0;
                        } else if (theta > eps_max) {
                            fth = 0.0;
                            fdth_fracth = 0.0;
                        } else {
                            fth = (2.0 - 2.0*c - th*s) / (2.0 * pow(th,4));
                            fdth_fracth = (-8.0 + (8.0 - th*th)*c + 5.0*th*s) / (2.0 * pow(th,6));
                         }
                        break;
                    case 4:
                         if (theta < eps_min) {
                            fth = 1.0/120.0;
                            fdth_fracth = -1.0/1260.0;
                        } else if (theta > eps_max) {
                            fth = 0.0;
                            fdth_fracth = 0.0;
                        } else {
                           fth = (2.0*th - 3.0*s + th*c) / (2.0 * pow(th,5));
                           fdth_fracth = (-8.0*th + (15.0 - th*th)*s - 7.0*th*c) / (2.0 * pow(th,7));
                        }
                        break;
                }

                MatrixXd Ith = MatrixXd::Zero(6,6);
                VectorXd Ith_comp = VectorXd::Zero(6);
                Ith_comp << 1,1,1,0,0,0;
                Ith = Ith_comp.asDiagonal();
                MatrixXd coadjop_pow = MatrixXd::Identity(6,6);
                for (int pp = 1 ; pp <= rr; pp++){
                    coadjop_pow = coadjop_pow * coadjop_fcn(omega);
                }
                SUM1 = SUM1 + pow(-1, rr) * fdth_fracth * coadjop_pow * FC * omega.transpose() * Ith * Z;
                MatrixXd SUM3 = MatrixXd::Zero(6, n[nd]);
                for (int uu = 1; uu <= rr; uu++){
                    MatrixXd coadjop_pow_u = MatrixXd::Identity(6,6);
                    MatrixXd coadjop_pow_ru = MatrixXd::Identity(6,6);
                    if (uu > 1){
                        for (int pp = 1 ; pp <= uu-1; pp++){
                            coadjop_pow_u = coadjop_pow_u * coadjop_fcn(omega);
                        }
                    }
                    if (rr - uu > 0){
                        for (int pp = 1 ; pp <= (rr - uu); pp++){
                            coadjop_pow_ru = coadjop_pow_ru * coadjop_fcn(omega);
                        }
                    }

                   SUM3 = SUM3 + coadjop_pow*coadjopbar_fcn(coadjop_pow_ru*FC)*Z;
                }
                SUM2 = SUM2 + pow(-1, rr) * fth * SUM3;
            }


            MatrixXd Bq_xi1 = MatrixXd::Zero(6,n[nd]);
            MatrixXd Bq_xi2 = MatrixXd::Zero(6,n[nd]);
            
            switch (nd) {
                case 0:
                    Bq_xi1 = compute_Bq1(X1,L);
                    Bq_xi2 = compute_Bq1(X2,L);
                  break;
                case 1:
                    Bq_xi1 = compute_Bq2(X1,L);
                    Bq_xi2 = compute_Bq2(X2,L);
                    break;
                case 2:
                    Bq_xi1 = compute_Bq3(X1,L);
                    Bq_xi2 = compute_Bq3(X2,L);
                 
                    break;
                case 3:
                    Bq_xi1 = compute_Bq4(X1,L);
                    Bq_xi2 = compute_Bq4(X2,L);    
                    break;
            }
            MatrixXd dSTdq_FC_T1 = sqrt(3.0) / 12.0 * pow(h[aa], 2) * (Bq_xi1.transpose() * coadjopbar_fcn(T.transpose() * FC) * Bq_xi2 - Bq_xi2.transpose() * coadjopbar_fcn(T.transpose() * FC) * Bq_xi1);
            MatrixXd dSTdq_FC = dSTdq_FC_T1 + Z.transpose() * SUM1 + Z.transpose() * SUM2;

            for (int kk = aa + 1; kk < g_pt[nd] + 2; kk++ ){
                MatrixXd grel = MatrixXd::Zero(4,4);
                MatrixXd gtemp = MatrixXd::Zero(4,4);
               
                switch (nd){
                    case 0:
                        grel = ginv(g.block(4*aa,0,4,4))*g.block(4*kk,0,4,4);
                        gtemp = g.block(4*kk,0,4,4);
                    break;
                    case 1:
                        grel = ginv(g.block(4*(g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-1]+2 + kk),0,4,4);
                        gtemp = g.block(4*(g_pt[nd-1]+2 + kk),0,4,4);
                    
                    break;
                    case 2:
                        grel = ginv(g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                        gtemp = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);

                    break;
                    case 3:
                        grel = ginv(g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                        gtemp = g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                    break;
                }
                MatrixXd Fcal = Ws[kk]*Mbar.block(6*kk,6*kk,6,6)*adjmap_fcn(ginv(gtemp))*grav6;
                MatrixXd SUM_Ad_ad_S = MatrixXd::Zero(6,n[nd]);

                for (int bb = 0; bb < kk - 1; bb++ ){
                    MatrixXd g_relb = MatrixXd::Zero(4,4);
                    MatrixXd gtemp = MatrixXd::Zero(4,4);
                    switch (nd){
                        case 0:
                            g_relb = ginv(g.block(4*bb,0,4,4))*g.block(4*kk,0,4,4);
                            gtemp = g.block(4*bb,0,4,4);
                        break;
                        case 1:
                            g_relb = ginv(g.block(4*(g_pt[nd-1]+2 + bb),0,4,4))*g.block(4*(g_pt[nd-1]+2 + kk),0,4,4);
                            gtemp = g.block(4*(g_pt[nd-1]+2 + bb),0,4,4);
                        break;
                        case 2:
                            g_relb = ginv(g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4))*g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                            gtemp = g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4);
                        break;
                        case 3:
                            g_relb = ginv(g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4))*g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                            gtemp = g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4);
                        break;  
                    }
                    double X1b = X[bb] + h[bb]/2.0 - sqrt(3.0)*h[bb]/6.0;
                    double X2b = X[bb] + h[bb]/2.0 + sqrt(3.0)*h[bb]/6.0;
                
                    VectorXd xi1b = compute_xi(qd_old, X1b, nd, n, L); 
                    VectorXd xi2b = compute_xi(qd_old, X2b, nd, n, L);
                    MatrixXd adjop_xi1b = adjop_fcn(xi1b);
                    VectorXd omega_b = h[bb]/2.0 * (xi1b + xi2b) + sqrt(3.0)*h[bb]*h[bb]/12.0 * (adjop_xi1b * xi2b);
                    double theta_b= omega_b.head(3).norm();
                    MatrixXd adjop_omega_b = adjop_fcn(omega_b);
                    MatrixXd T_b = compute_T(theta_b, adjop_omega_b);
                    MatrixXd Z_b = compute_Z(h[bb], xi1b, xi2b, nd, X1b, X2b, L); 
                    MatrixXd S_b = T_b*Z_b;

                    SUM_Ad_ad_S = SUM_Ad_ad_S + adjmap_fcn(ginv(g_relb))*adjop_fcn(adjmap_fcn(ginv(gtemp)) * grav6) * S_b;
                }



                MatrixXd dFcal_dq = Ws[kk]*Mbar.block(6*kk,6*kk,6,6)*SUM_Ad_ad_S;

                MatrixXd dAd_dq_Fcal = MatrixXd::Zero(6,n.sum());
              
                for (int bb = aa; bb < kk -1 ; bb ++){
                    MatrixXd g_relb = MatrixXd::Zero(4,4);
                    MatrixXd g_rel_ab = MatrixXd::Zero(4,4);
                    switch (nd){
                        case 0:
                            g_relb = ginv(g.block(4*bb,0,4,4))*g.block(4*kk,0,4,4);
                            g_rel_ab = ginv(g.block(4*aa,0,4,4))*g.block(4*bb,0,4,4);
                        break;
                        case 1:
                            g_relb = ginv(g.block(4*(g_pt[nd-1]+2 + bb),0,4,4))*g.block(4*(g_pt[nd-1]+2 + kk),0,4,4);
                            g_rel_ab = ginv(g.block(4*(g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-1]+2 + bb),0,4,4);
                        break;
                        case 2:
                            g_relb = ginv(g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4))*g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                            g_rel_ab = ginv(g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4);
                        break;
                        case 3:
                            g_relb = ginv(g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4))*g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + kk),0,4,4);
                            g_rel_ab = ginv(g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + aa),0,4,4))*g.block(4*(g_pt[nd-3]+2 + g_pt[nd-2]+2 + g_pt[nd-1]+2 + bb),0,4,4);
                        break;
                    }
                    double X1b = X[bb] + h[bb]/2.0 - sqrt(3.0)*h[bb]/6.0;
                    double X2b = X[bb] + h[bb]/2.0 + sqrt(3.0)*h[bb]/6.0;
                    VectorXd xi1b = compute_xi(qd_old, X1b, nd, n, L);
                    VectorXd xi2b = compute_xi(qd_old, X2b, nd, n, L);
                    MatrixXd adjop_xi1b = adjop_fcn(xi1b);
                    VectorXd omega_b = h[bb]/2.0 * (xi1b + xi2b) + sqrt(3.0)*h[bb]*h[bb]/12.0 * (adjop_xi1b * xi2b);
                    double theta_b= omega_b.head(3).norm();
                    MatrixXd adjop_omega_b = adjop_fcn(omega_b);
                    MatrixXd T_b = compute_T(theta_b, adjop_omega_b);
                    MatrixXd Z_b = compute_Z(h[bb], xi1b, xi2b, nd, X1b, X2b, L); 
                    MatrixXd S_b = T_b*Z_b;
                    MatrixXd Mcomp = MatrixXd::Zero(6,n.sum());
                    switch (nd) {
                        case 0:
                            Mcomp.block(0,0,6,n[0]) = coadjopbar_fcn(coadjmap_fcn(g_relb)*Fcal)*S_b;
                        break;
                        case 1:
                            Mcomp.block(0,n[0],6,n[1]) = coadjopbar_fcn(coadjmap_fcn(g_relb)*Fcal)*S_b;
                        break;
                        case 2:
                            Mcomp.block(0,n[0]+n[1],6,n[2]) = coadjopbar_fcn(coadjmap_fcn(g_relb)*Fcal)*S_b;
                        break;
                        case 3:  
                            Mcomp.block(0,n[0]+n[1]+n[2],6,n[3]) = coadjopbar_fcn(coadjmap_fcn(g_relb)*Fcal)*S_b;
                        break;
                    }
                    dAd_dq_Fcal = dAd_dq_Fcal + coadjmap_fcn(g_rel_ab)*Mcomp;
                }
                SUM_dAddq_Fcal = SUM_dAddq_Fcal + dAd_dq_Fcal;

                MatrixXd Mcomp = MatrixXd::Zero(6,n.sum());
                switch (nd) {
                    case 0:
                        Mcomp.block(0,0,6,n[0]) = dFcal_dq;
                    break;
                    case 1:
                        Mcomp.block(0,n[0],6,n[1]) = dFcal_dq;
                    break;
                    case 2:
                        Mcomp.block(0,n[0]+n[1],6,n[2]) = dFcal_dq;
                    break;
                    case 3:  
                        Mcomp.block(0,n[0]+n[1]+n[2],6,n[3]) = dFcal_dq;
                    break;
                }
                SUM_Ad_dFcaldq = SUM_Ad_dFcaldq + coadjmap_fcn(grel)*Mcomp;
                       
            }

            MatrixXd dFalpha_dq = MatrixXd::Zero(n[nd],n.sum());
            MatrixXd Mcomp = MatrixXd::Zero(n[nd],n.sum());  
            switch (nd) {
                case 0:
                    Mcomp.block(0,0,n[nd],n[nd]) = dSTdq_FC;
                break;
                case 1:
                    Mcomp.block(0,n[nd-1],n[nd],n[nd]) =dSTdq_FC ;
                break;
                case 2:
                    Mcomp.block(0,n[nd-1]+n[nd-2],n[nd],n[nd]) =dSTdq_FC ;
                break;
                case 3:  
                    Mcomp.block(0,n[nd-1]+n[nd-2]+n[nd-3],n[nd],n[nd]) =dSTdq_FC;
                break;
            }
            dFalpha_dq = Mcomp + S.transpose() * SUM_dAddq_Fcal + S.transpose() * SUM_Ad_dFcaldq;
            dFdq = dFdq + dFalpha_dq;
      
        }

        return dFdq;
    
    }