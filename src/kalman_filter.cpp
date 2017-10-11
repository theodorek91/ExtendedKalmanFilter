#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
	MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in; // object state
	P_ = P_in; // object covariance matrix
	F_ = F_in; // state transition matrix
	H_ = H_in; // measurement matrix
	R_ = R_in; // measurement covariance matrix
	Q_ = Q_in; // process covariance matrix
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	x_ = F_ * x_; 
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	// For laser
	VectorXd y = z - H_ * x_;

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
	
	// New estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	// for radar due to non linearity
	
	const float  PI = 3.14159265;
	
	// Required Calculations to check of consistency of phi to lie between pi and -pi
	VectorXd z_copy(3);
	z_copy << z(0), z(1), z(2);

	if (z_copy(1) > PI) { z_copy(1) = z_copy(1) - (2 * PI); }
	else if (z_copy(1) < -PI) { z_copy(1) = (2 * PI) + z_copy(1); }
	// calculate h(x')
	float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
	float phi = atan2(x_(1), x_(0));
	float rho_dot = 0.0;
	if (fabs(rho) > 1e-6) { rho_dot = (x_(0)*x_(2) + x_(1) * x_(3)) / rho; }

	VectorXd hx(3);
	hx << rho, phi, rho_dot;

	VectorXd y = z_copy - hx;

	// check to make sure phi lies between pi and -pi
	while (y(1) > PI) { y(1) = y(1) - (2 *PI) ; }
	while (y(1) < -PI) { y(1) = y(1) + (2 * PI); }


	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	// New estimate
	
	x_ = x_ + (K * y);
	
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I - K * H_) * P_;


}
