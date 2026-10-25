// ============================================================================
// WARNING: RECENTLY UPDATED, NOT YET VALIDATED ON THE ROBOT
//   - motor command saturation added (TAU_SAT, POS_SAT)
//   - tendon insets, friction coefficients and rest cable lengths (Lc0)
//     updated to robot3 values
//   - CMakeLists.txt cleaned: only the files in src/ are built
// Double check this file before running it on the hardware.
// ============================================================================
#include "md80_node.hpp"
#include "std_msgs/msg/float64.hpp"
#include "general_fcns.h"
#include "softrob_gvs.h"
#include "softrob_CCandCLIK.h"

#include "mocap_optitrack_interfaces/msg/marker_array.hpp"
#include "std_msgs/msg/float64_multi_array.hpp"
//#include <cmath>

#define TAU_SAT 0.6            // symmetric saturation on cmd_tau  [N·m at motor]
#define POS_SAT (M_PI / 2.0)   // saturation on cmd_pos            [rad]
// NOTE: if strictly necessary, POS_SAT can be extended up to M_PI.
// This is borderline for the arm: do it only if the task really requires it.

class MyControlNode : public rclcpp::Node
{
public:
    MyControlNode(const MatrixXd& xd_all)
    : Node("my_control"), base_id(5), marker3_id(349) ,tip_id(365), xd_all_(xd_all), wait_steps(0),intThadTha(VectorXd::Zero(4)), pos_m1(0.0), pos_m2(0.0),  pos_m3(0.0), pos_m4(0.0), vel_m1(0.0), vel_m2(0.0), vel_m3(0.0), vel_m4(0.0),qd(VectorXd::Zero(17)),traj_col(0),traj_pos(1), mbx(0.0), mby(0.0), mbz(0.0), mtx(0.0), mty(0.0), mtz(0.0), m3x(0.0), m3y(0.0), m3z(0.0),mbxlcl(0.0), mbylcl(0.0), mbzlcl(0.0), mtxlcl(0.0), mtylcl(0.0), mtzlcl(0.0), m3xlcl(0.0), m3ylcl(0.0), m3zlcl(0.0)
    {
        start_time_ = this->now();
        
        // Publisher
        xd1_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xd1", 1);
        xd2_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xd2", 1);
        xd3_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xd3", 1);
        xd4_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xd4", 1);
        xe1_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xe1", 1);
        xe2_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xe2", 1);
        xe3_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xe3", 1);
        xe4_pub_ = this->create_publisher<std_msgs::msg::Float64>("/xe4", 1);
    
        current_time_pub_ = this->create_publisher<std_msgs::msg::Float64>("/current_time", 1);
  

        torque1_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor1/torque", 1);
        ref_torque1_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor1/Ref_torque", 1); 
        torque2_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor2/torque", 1);
        ref_torque2_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor2/Ref_torque", 1);
        torque3_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor3/torque", 1);
        ref_torque3_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor3/Ref_torque", 1);     
		torque4_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor4/torque", 1);
        ref_torque4_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor4/Ref_torque", 1);
      
        Lc1_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor1/Lc", 1);
        Lc1d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor1/Lcd", 1);
        Lc2_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor2/Lc", 1);
        Lc2d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor2/Lcd", 1);
        Lc3_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor3/Lc", 1);
        Lc3d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor3/Lcd", 1);
        Lc4_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor4/Lc", 1);
        Lc4d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor4/Lcd", 1);
        ybar1_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor1/ybar", 1);
        ybar1d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor1/yd_bar", 1);
        ybar2_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor2/ybar", 1);
        ybar2d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor2/yd_bar", 1);
        ybar3_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor3/ybar", 1);
        ybar3d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor3/yd_bar", 1);
        ybar4_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor4/ybar", 1);
        ybar4d_pub_ = this->create_publisher<std_msgs::msg::Float64>("/motor4/yd_bar", 1);

        qd_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>("/qd_vec", 1);

        motion_command_pub_ = this->create_publisher<candle_ros2::msg::MotionCommand>("/md80/motion_command", 1);


        // Subscriber
        joint_states_sub_ = this->create_subscription<sensor_msgs::msg::JointState>(
            "/md80/joint_states", 1, std::bind(&MyControlNode::jointStatesCallback, this, std::placeholders::_1));

       marker_sub_ = this->create_subscription<mocap_optitrack_interfaces::msg::MarkerArray>(
            "/mocap_markers", 1, std::bind(&MyControlNode::markerCallback, this, std::placeholders::_1));
            

        timer_ = this->create_wall_timer(
            std::chrono::milliseconds(10),
            std::bind(&MyControlNode::timerCallback, this));
    }

private:
    rclcpp::Time start_time_;
    MatrixXd xd_all_;
    double pos_m1, pos_m2, pos_m3, pos_m4;
    double vel_m1, vel_m2, vel_m3, vel_m4;
    double mbx,mby,mbz;
    double mtx,mty,mtz;
    double m3x,m3y,m3z;
	VectorXd intThadTha,qd;
    int traj_col,traj_pos,wait_steps;

    double mbxlcl,mbylcl,mbzlcl;
    double mtxlcl,mtylcl,mtzlcl;
    double m3xlcl,m3ylcl,m3zlcl;
    int base_id, marker3_id, tip_id;

    void markerCallback(const mocap_optitrack_interfaces::msg::MarkerArray::SharedPtr msg)
    {
        //Debug logging
        //RCLCPP_INFO(this->get_logger(), "Received MarkerArray with %zu markers", msg->markers.size());
        //for (const auto &marker : msg->markers)
        //{
        //    RCLCPP_INFO(this->get_logger(), "Marker ID: %d, Type: %s, Position: [%.2f, %.2f, %.2f], Size: %.2f",
        //                marker.id, marker.type.c_str(), marker.position.x, marker.position.y, marker.position.z, marker.size);
        //}
    
        
//    //THE IDs NEED TO BE CHANGED EVERY TIME!!!!
        
//     int base_id = 52228;
//     //const int marker1_id = 50945;
//     //const int marker2_id = 50949;
//     int marker3_id = 52262;
//     int tip_id = 52519;
        double thresh_dist = 0.03; // Threshold distance to consider a marker "found"
        
        bool found_base = false, found_marker3 = false, found_tip = false;
        Vector3d base_pos_prev(mbxlcl, mbylcl, mbzlcl);
        Vector3d marker3_pos_prev(m3xlcl, m3ylcl, m3zlcl);
        Vector3d tip_pos_prev(mtxlcl, mtylcl, mtzlcl);
         for (const auto &marker : msg->markers) {
            if (marker.id == base_id) found_base = true;
            if (marker.id == marker3_id) found_marker3 = true;
            if (marker.id == tip_id) found_tip = true;
        }

        if (!found_base) {
            double min_dist = std::numeric_limits<double>::max();
            int new_id = -1;
            for (const auto &marker : msg->markers) {
                Vector3d pos(marker.position.x, marker.position.y, marker.position.z);
                double dist = (pos - base_pos_prev).norm();
                if (dist < min_dist) {
                    min_dist = dist;
                    new_id = marker.id;
                }
            }
            if (new_id != -1 && min_dist < thresh_dist) base_id = new_id;
        }
        if (!found_marker3) {
            double min_dist = std::numeric_limits<double>::max();
            int new_id = -1;
            for (const auto &marker : msg->markers) {
                Vector3d pos(marker.position.x, marker.position.y, marker.position.z);
                double dist = (pos - marker3_pos_prev).norm();
                if (dist < min_dist) {
                    min_dist = dist;
                    new_id = marker.id;
                }
            }
            if (new_id != -1 && min_dist < thresh_dist) marker3_id = new_id;
        }
        if (!found_tip) {
            double min_dist = std::numeric_limits<double>::max();
            int new_id = -1;
            for (const auto &marker : msg->markers) {
                Vector3d pos(marker.position.x, marker.position.y, marker.position.z);
                double dist = (pos - tip_pos_prev).norm();
                if (dist < min_dist) {
                    min_dist = dist;
                    new_id = marker.id;
                }
            }
            if (new_id != -1 && min_dist < thresh_dist) tip_id = new_id;
        }

         Vector3d offset;
         //offset << -0.021,- 0.008, 0.0;
        offset << -0.03,- 0.013, 0.0;
         
         Matrix3d R_alg;
        R_alg << 1, 0, 0,
                 0, 0, 1,
                 0, -1, 0;

        double temp_x = 0.0;
        double temp_y = 0.0;
        double temp_z = 0.0;
        double tempb_x = 0.0;
        double tempb_y = 0.0;
        double tempb_z = 0.0;
        
        Vector3d temp_pos,tempb_Rpos,pos_gl;
        temp_pos.setZero();
        tempb_Rpos.setZero();
        pos_gl.setZero();

        //Looking for the base marker
        for (const auto &marker : msg->markers)
        {

            if (marker.id == base_id)
            {
                tempb_x = marker.position.x;
                tempb_y = marker.position.y;
                tempb_z = marker.position.z;

                 mbxlcl = marker.position.x;
                 mbylcl = marker.position.y;
                 mbzlcl = marker.position.z;
                //cout << "base marker lcl:" << tempb_x << " " << tempb_y << " " << tempb_z << endl;
              
                temp_pos <<-offset[0], -offset[1], -offset[2];
                pos_gl = R_alg*temp_pos;

                mbx = pos_gl[0];
                mby = pos_gl[1];
                mbz = pos_gl[2];

            }
        }

        // Save the other markers
        for (const auto &marker : msg->markers)
        {  
             if(marker.id == marker3_id){

                temp_x = marker.position.x;
                temp_y = marker.position.y;
                temp_z = marker.position.z;

                 m3xlcl = marker.position.x;
                 m3ylcl = marker.position.y;
                 m3zlcl = marker.position.z;

                //cout<< "marker 3 lcl:" << temp_x << " " << temp_y << " " << temp_z << endl;
                temp_pos<< temp_x-(tempb_x+offset[0]),
                           temp_y-(tempb_y+offset[1]),
                           temp_z-(tempb_z+offset[2]);
                pos_gl = R_alg*temp_pos;
                m3x = pos_gl[0];
                m3y = pos_gl[1];
                m3z = pos_gl[2];

             }
           else if (marker.id == tip_id) {
                temp_x = marker.position.x;
                temp_y = marker.position.y;
                temp_z = marker.position.z;

                 mtxlcl = marker.position.x;
                 mtylcl = marker.position.y;
                 mtzlcl = marker.position.z;

                //cout<< "tip marker lcl:" << temp_x << " " << temp_y << " " << temp_z << endl;
                temp_pos<< temp_x-(tempb_x+offset[0]),
                           temp_y-(tempb_y+offset[1]),
                           temp_z-(tempb_z+offset[2]);

                pos_gl = R_alg*temp_pos;
                mtx = pos_gl[0];
                mty = pos_gl[1];
                mtz = pos_gl[2];
            }
        }
        //cout << "mbx:" << mbx << " mby:" << mby << " mbz:" << mbz << endl; 
        //cout << "m1x:" << m1x << " m1y:" << m1y << " m1z:" << m1z << endl;
        //cout << "m2x:" << m2x << " m2y:" << m2y << " m2z:" << m2z << endl;  
        //cout << "m3x:" << m3x << " m3y:" << m3y << " m3z:" << m3z << endl;
        //cout << "mtx:" << mtx << " mty:" << mty << " mtz:" << mtz << endl;

}

    // Callback for joint_states topc
    void jointStatesCallback(const sensor_msgs::msg::JointState::SharedPtr msg)
    {
        if (!msg->position.empty()) {
            std_msgs::msg::Float64 torque1_msg,torque2_msg,torque3_msg,torque4_msg;       
            
            pos_m1 = msg->position[0];
            pos_m2 = msg->position[1];
            pos_m3 = msg->position[2];
            pos_m4 = msg->position[3];
            
            vel_m1 = msg->velocity[0];
            vel_m2 = msg->velocity[1];
            vel_m3 = msg->velocity[2];
            vel_m4 = msg->velocity[3];
            
            torque1_msg.data = msg->effort[0];  
            torque2_msg.data = msg->effort[1];
            torque3_msg.data = msg->effort[2];  
            torque4_msg.data = msg->effort[3];
            
            torque1_pub_->publish(torque1_msg);
            torque2_pub_->publish(torque2_msg);
            torque3_pub_->publish(torque3_msg);
            torque4_pub_->publish(torque4_msg);
        }
    }

    // Control Callback 
    void timerCallback()
    {

    // Compute currente time
    auto now = this->now();
    double elapsed = (now - start_time_).seconds();
    std_msgs::msg::Float64 time_msg;
    time_msg.data = elapsed;
    current_time_pub_->publish(time_msg);
    
    double radius = 0.0325;
    int n_divisions = 4;

    VectorXi n(n_divisions);
	n << 3, 8, 4, 2;
	int na=4;
    
    VectorXd L(4);     // Vector for lengths
    L << 0.025, 0.255, 0.025, 0.055;
    
    double E = 302000;
    double v = 0.4;
    double ups = 5000;   // material damping [Pa*s] (see data/soft_arm_parameters.yaml)
    double rho = 1310;


    double G = E / (2 * (1 + v));
    VectorXd mu (na); 

    //old
	//mu << 0.11,0.24,0.12,0.125;  //friction mu_f
	//mu << 0.47,0.835,0.437,0.523; //friction mu_x
	
    // mu << 0.0651,0.145,0.0353,0.055;  //friction mu_f (thesis arm)
    mu << 0.003,0.021,0.005,0.023;  // friction_model_phi, robot3 (see data/soft_arm_parameters.yaml)
	//mu << 0.2862,0.5424,0.1483,0.2493; //friction mu_x
	


    VectorXd X,Ws;
    VectorXd h;
    MatrixXd g,J;
    VectorXd Lc(na),Lcdot(na),Lcd(na),ybar(na),yd_bar(na);
    yd_bar.setZero();
    ybar.setZero();
    Lc.setZero();
    Lcd.setZero(); 
    Lcdot.setZero(); 
    

    double L0=0.02;
    double Rb0 = 0.016;
    double Rt4 = 0.0048;
    double Rt0 = (Rb0 + (Rt4 - Rb0) * (L0 / (L0+L.sum())));
 
    // Assuming a continuous conical decrease
    double Rt1 = (Rb0 + (Rt4 - Rb0) * ((L0+L[0]) / (L0+L.sum())));
    double Rt2 = (Rb0 + (Rt4 - Rb0) * ((L0+L[0] + L[1]) / (L0+L.sum())));
    double Rt3 = (Rb0 + (Rt4 - Rb0) * ((L0+L[0] + L[1] + L[2]) / (L0+L.sum())));

    double Rb1 = Rt0, Rb2 = Rt1, Rb3 = Rt2, Rb4 = Rt3;


    VectorXd Rb(4);
    Rb << Rb1, Rb2, Rb3, Rb4;
    VectorXd Rt(4);
    Rt << Rt1, Rt2, Rt3, Rt4;

    // Tendon insets, robot3 (see data/soft_arm_parameters.yaml)
    // Thesis arm used [0.004, 0.004, 0.003, 0.003]
    double d1_1 = Rb[0] - 0.002;
    double d2_1 = Rb[1] - 0.004;
    double d3_1 = Rb[2] - 0.002;
    double d4_1 = Rb[3] - 0.001;

    double d1_23 = Rb[0] - 0.002;
    double d2_23 = Rb[1] - 0.004;
    double d3_23 = -(Rb[2] - 0.002);
    double d4_23 = -(Rb[3] - 0.001);
        
    MatrixXd dist(2, 4);
    dist << d1_1, d2_1, d3_1, 0,
            d1_23, d2_23, d3_23, 0;

    double m1_1 = (d2_1 - d1_1) / L[0];
    double m2_1 = (d3_1 - d2_1) / L[1];
    double m3_1 = (d4_1 - d3_1) / L[2];

    double m1_23 = (d2_23 - d1_23) / L[0];
    double m2_23 = (d3_23 - d2_23) / L[1];
    double m3_23 = (d4_23 - d3_23) / L[2];

    MatrixXd m(2, 4);
    m << m1_1, m2_1, m3_1, 0,
         m1_23, m2_23, m3_23, 0;

    VectorXd alpha(4);
    alpha << 30 * M_PI / 180, 60 * M_PI / 180, 120 * M_PI / 180, 150 * M_PI / 180;

    // Set time interval
    double tsample = 0.01; // 0.02;  
 
 
    VectorXd xd = VectorXd::Zero(na);   

    double threshold_p = 4e-3; 
    double threshold_gains = 2e-2; 

    VectorXd xe = VectorXd::Zero(na);
    xe<< mtx, mty,mtz,m3y;

    // triangle
    if (traj_col<=4000) { // i need about 40s for each point to converge with variable gains 2cm -> CON IL JUMP ERA PIU VELOCE -> PROVARE QUELLO?
        //xd<<0.12,0.1,-0.27,0.04;
        xd<<0.12,-0.1,-0.27,-0.04;
        traj_col++;
    }   
    else if(traj_col<=8000) {
        xd<<0.12,0,-0.318,-0.002;        
        traj_col++;
    }
    else if(traj_col<12000) {
        xd<<0.12,0.1,-0.27,0.04;    
        //xd<<0.12,-0.1,-0.27,-0.04;
        traj_col++;
    }
    else {
        xd<<0.12,0.1,-0.27,0.04;       
    }
    


    // bool change_g = false;
    // if ((xe - xd).norm() < threshold_gains){
    //     change_g = true;
    // }

    cout << "error norm:" << (xe - xd).norm() << endl;


    cout <<"traj_col:" << traj_col << endl;
    cout << "xd:" << endl << xd << endl;
    cout << "xe:" << endl << xe << endl;
    

    qd = clik_rt(n_divisions,n,na,L,qd,xd,xe,tsample,dist,m,alpha,Rb,Rt,E,G,rho,mu,0);
    // qd = clik_rt_analytic(n_divisions,n,na,L,qd,xd,xe,tsample,dist,m,alpha,Rb,Rt,E,G,rho,mu,change_g);
    cout << "qd:" << endl << qd << endl;

    VectorXd u(na);
    u.setZero();
    
    VectorXd delta_L(na),Lc0(na);
    Lc0.setZero();
    delta_L.setZero();


    //Lc0 << 0.305105,0.305192,0.305192,0.305105;    
    //Lc0 << 0.306721,0.306807,0.306807,0.306721;  
    
    // Lc0 << 0.306614,0.306701,0.306701,0.306614;   // thesis arm insets
    // Rest cable lengths, robot3 insets (from ConstrainedCLIK_PIDffdCC_3DOA_rt):
    // tendons 1 and 4 (straight) / tendons 2 and 3 (crossing)
    Lc0 << 0.306722,0.306829,0.306829,0.306722;

    delta_L << pos_m1*radius,pos_m2*radius,pos_m3*radius,pos_m4*radius;
	Lc = Lc0 + delta_L;
    
   
    //Lcdot << vel_m1*radius,vel_m2*radius,vel_m3*radius,vel_m4*radius; 
  
    input_collocatedCtrlPsIDFFd_MIXED(n_divisions,n,na,L,tsample,dist,m,alpha,qd,intThadTha,Rb,Rt,E,G,rho,mu,u,Lc,Lcdot,Lcd,ybar,yd_bar);
    float pos1d = (Lcd[0]-Lc0[0])/radius;
    float pos2d = (Lcd[1]-Lc0[1])/radius;
    float pos3d = (Lcd[2]-Lc0[2])/radius;
    float pos4d = (Lcd[3]-Lc0[3])/radius;


    std_msgs::msg::Float64MultiArray qd_msg;
    qd_msg.data.resize(qd.size());
    for (int i = 0; i < qd.size(); ++i) {
        qd_msg.data[i] = qd[i];
    }
    qd_pub_->publish(qd_msg);

    // saturation
    //for(int i =0; i<na; i++){
    //    if(u[i]>0){
    //         u[i] = 0.0;
    //     }
    //}
    
    // friction compensation at low velocity
    //float vel_threshold = 0.5f;
	//float tauf_pos = 0.025f;//23
	//float tauf_neg = -0.025f;//22    		
    //if(vel_m1>=-vel_threshold && vel_m1<=vel_threshold ){
    //	if(u[0]<0)
    //    u[0] += tauf_neg;
    //	else if(u[0]>=0)
    //    u[0] += tauf_pos;
    //}	
    //if(vel_m2>=-vel_threshold && vel_m2<=vel_threshold ){
    //	if(u[1]<0)
    //    u[1] += tauf_neg;
    //	else if(u[1]>=0)
    //    u[1] += tauf_pos;
    //}




    // --- Safety: saturate motor commands (same logic as ConstrainedCLIK_PIDffdCC_3DOA_rt) --- //
    // Position: clamp to [-POS_SAT, 0] (negative = pulling the tendon).
    // Torque:   clamp to [-TAU_SAT, +TAU_SAT]; if the position is saturated at its
    //           upper bound (0), the upper torque limit becomes 0 so the motor
    //           cannot push the tendon further in the releasing direction.
    VectorXd cmd_pos(na), cmd_tau(na);
    cmd_pos << pos1d, pos2d, pos3d, pos4d;
    cmd_tau << u[0]*radius, u[1]*radius, u[2]*radius, u[3]*radius;

    for (int i = 0; i < na; ++i) {
        cmd_pos(i) = std::max(-POS_SAT, std::min(cmd_pos(i), 0.0));

        const bool pos_at_upper = (cmd_pos(i) >= 0.0);   // upper bound is 0.0
        const double tau_upper = pos_at_upper ? 0.0 : TAU_SAT;
        cmd_tau(i) = std::max(-TAU_SAT, std::min(cmd_tau(i), tau_upper));
    }

    pos1d = cmd_pos(0);
    pos2d = cmd_pos(1);
    pos3d = cmd_pos(2);
    pos4d = cmd_pos(3);

    float ref_torque1 = cmd_tau(0);
    float ref_torque2 = cmd_tau(1);
    float ref_torque3 = cmd_tau(2);
    float ref_torque4 = cmd_tau(3);
    
        // Publications
  		candle_ros2::msg::MotionCommand motion_msg;
  	
	  	//ADAPT ACCORDING TO MOTORS ID...
        motion_msg.drive_ids.push_back(277); 
        motion_msg.drive_ids.push_back(278);
        motion_msg.drive_ids.push_back(280);
        motion_msg.drive_ids.push_back(279);
    

    if (wait_steps > 10) {
        motion_msg.target_position.push_back(pos1d);
        motion_msg.target_position.push_back(pos2d);
        motion_msg.target_position.push_back(pos3d);
        motion_msg.target_position.push_back(pos4d);
            
        motion_msg.target_velocity.push_back(0.0);
        motion_msg.target_velocity.push_back(0.0);       
        motion_msg.target_velocity.push_back(0.0);
        motion_msg.target_velocity.push_back(0.0);

        motion_msg.target_torque.push_back(ref_torque1);
        motion_msg.target_torque.push_back(ref_torque2);
        motion_msg.target_torque.push_back(ref_torque3);
        motion_msg.target_torque.push_back(ref_torque4);
    }
    else {
        wait_steps++;

        motion_msg.target_position.push_back(0);
        motion_msg.target_position.push_back(0);
        motion_msg.target_position.push_back(0);
        motion_msg.target_position.push_back(0);

        motion_msg.target_velocity.push_back(0.0);
        motion_msg.target_velocity.push_back(0.0);       
        motion_msg.target_velocity.push_back(0.0);
        motion_msg.target_velocity.push_back(0.0);
        
        motion_msg.target_torque.push_back(0);
        motion_msg.target_torque.push_back(0);
        motion_msg.target_torque.push_back(0);
        motion_msg.target_torque.push_back(0);
    }
    
        motion_command_pub_->publish(motion_msg);

        // Pubblica la posizione di riferimento
        std_msgs::msg::Float64 ref_torque1_msg,ref_torque2_msg,ref_torque3_msg,ref_torque4_msg;
        ref_torque1_msg.data = ref_torque1;
        ref_torque1_pub_->publish(ref_torque1_msg);
        ref_torque2_msg.data = ref_torque2;
        ref_torque2_pub_->publish(ref_torque2_msg);
        ref_torque3_msg.data = ref_torque3;
        ref_torque3_pub_->publish(ref_torque3_msg);
        ref_torque4_msg.data = ref_torque4;
        ref_torque4_pub_->publish(ref_torque4_msg);     
        
        std_msgs::msg::Float64 Lc1_msg,Lc2_msg,Lc3_msg,Lc4_msg,Lc1d_msg,Lc2d_msg,Lc3d_msg,Lc4d_msg;
        Lc1_msg.data = Lc[0];
        Lc1_pub_->publish(Lc1_msg);
        Lc1d_msg.data = Lcd[0];
        Lc1d_pub_->publish(Lc1d_msg);
        
        Lc2_msg.data = Lc[1];
        Lc2_pub_->publish(Lc2_msg);
        Lc2d_msg.data = Lcd[1];
        Lc2d_pub_->publish(Lc2d_msg);
        
        Lc3_msg.data = Lc[2];
        Lc3_pub_->publish(Lc3_msg);
        Lc3d_msg.data = Lcd[2];
        Lc3d_pub_->publish(Lc3d_msg);
        
        Lc4_msg.data = Lc[3];
        Lc4_pub_->publish(Lc4_msg);
        Lc4d_msg.data = Lcd[3];
        Lc4d_pub_->publish(Lc4d_msg);     
        
        std_msgs::msg::Float64 ybar1_msg,ybar2_msg,ybar3_msg,ybar4_msg,ybar1d_msg,ybar2d_msg,ybar3d_msg,ybar4d_msg;
        ybar1_msg.data = ybar[0];
        ybar1_pub_->publish(ybar1_msg);
        ybar1d_msg.data = yd_bar[0];
        ybar1d_pub_->publish(ybar1d_msg);
        
        ybar2_msg.data = ybar[1];
        ybar2_pub_->publish(ybar2_msg);
        ybar2d_msg.data = yd_bar[1];
        ybar2d_pub_->publish(ybar2d_msg);
        
        ybar3_msg.data = ybar[2];
        ybar3_pub_->publish(ybar3_msg);
        ybar3d_msg.data = yd_bar[2];
        ybar3d_pub_->publish(ybar3d_msg);
        
        ybar4_msg.data = ybar[3];
        ybar4_pub_->publish(ybar4_msg);
        ybar4d_msg.data = yd_bar[3];
        ybar4d_pub_->publish(ybar4d_msg);     
        

        std_msgs::msg::Float64 xd1_msg, xd2_msg, xd3_msg,xd4_msg, xe1_msg, xe2_msg,xe3_msg, xe4_msg;
        xd1_msg.data = xd[0];
        xd2_msg.data = xd[1];
        xd3_msg.data = xd[2];
        xd4_msg.data = xd[3];
        xe1_msg.data = xe[0];
        xe2_msg.data = xe[1];
        xe3_msg.data = xe[2];
        xe4_msg.data = xe[3];

        // Pubblica i messaggi
        xd1_pub_->publish(xd1_msg);
        xd2_pub_->publish(xd2_msg);
        xd3_pub_->publish(xd3_msg);
        xd4_pub_->publish(xd4_msg);
        xe1_pub_->publish(xe1_msg);
        xe2_pub_->publish(xe2_msg);
        xe3_pub_->publish(xe3_msg);
        xe4_pub_->publish(xe4_msg);
        
              
    }

    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xd1_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xd2_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xd3_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xd4_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xe1_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xe2_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xe3_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr xe4_pub_;


    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr torque1_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ref_torque1_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr torque2_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ref_torque2_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr torque3_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ref_torque3_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr torque4_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ref_torque4_pub_;
 
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc1_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc1d_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc2_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc2d_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc3_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc3d_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc4_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr Lc4d_pub_;    
    
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar1_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar1d_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar2_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar2d_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar3_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar3d_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar4_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr ybar4d_pub_;    
 
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr current_time_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr qd_pub_;

    rclcpp::Publisher<candle_ros2::msg::MotionCommand>::SharedPtr motion_command_pub_;
    rclcpp::Subscription<sensor_msgs::msg::JointState>::SharedPtr joint_states_sub_;
    rclcpp::Subscription<mocap_optitrack_interfaces::msg::MarkerArray>::SharedPtr marker_sub_; 
    rclcpp::TimerBase::SharedPtr timer_;
};

int main(int argc, char **argv)
{
    // Initialize ROS 2 node
    rclcpp::init(argc, argv);

    MatrixXd xd_all = MatrixXd::Zero(4, 1);

    // Load trajectory from file
    // MatrixXd xd_all;
    //try {
    //    xd_all = loadMatrixFromTxt("/home/vitodperfetta/Desktop/md_motor/code/ros2_ws/src/candle_ros2/src/xdm34_all.txt");
    //} catch (const std::exception& e) {
    //    std::cerr << "Error loading xd_all: " << e.what() << std::endl;
    //    return -1;
    //}

    // Control node spinng
    rclcpp::spin(std::make_shared<MyControlNode>(xd_all));

    // Closing ROS 2
    rclcpp::shutdown();
    return 0;
}