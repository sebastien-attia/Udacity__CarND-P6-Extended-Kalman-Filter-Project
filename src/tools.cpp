#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);

	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	//check division by zero
	double dist2 = pow(px, 2)+pow(py, 2);
	double dist = sqrt(dist2);
  double dist15 = dist2 * dist;
	if (isZero(dist2)) {
	    cout << "CalculateJacobian(): Division by 0" << endl;
      return Hj;
	}

	//compute the Jacobian matrix
  Hj << px/dist,  py/dist,  0, 0,
       -py/dist2, px/dist2, 0, 0,
        py*(vx*py-vy*px)/dist15, px*(vy*px - vx*py)/dist15, px/dist, py/dist;

	return Hj;
}

VectorXd Tools::ConvertToCartesian(const VectorXd& x_state) {
  double rho = x_state[0];
  double theta = x_state[1];
  double rho_dot = x_state[2];

  double cos_theta = cos(theta);
  double sin_theta = sin(theta);

  double x = cos_theta * rho;
  double y = sin_theta * rho;
  double vx = rho_dot * cos_theta;
  double vy = rho_dot * sin_theta;

  VectorXd result = VectorXd(4);
  result << x, y, vx, vy;
  return result;
}

VectorXd Tools::ConvertToPolar(const VectorXd& x_state) {
  double x = x_state[0];
  double y = x_state[1];
  double vx = x_state[2];
  double vy = x_state[3];

  double rho = sqrt(pow(x, 2) + pow(y, 2));
  double theta = atan2(y, x);
  double rho_dot = (x*vx + y*vy)/rho;

  VectorXd result = VectorXd(3);
  result << rho, theta, rho_dot;
  return result;
}

double Tools::rangeAngle(const double theta) {
  double theta_rem = std::remainder(theta + M_PI, 2 * M_PI);
  return (theta_rem >= 0 ? (theta_rem - M_PI) : (theta_rem + M_PI));
}

bool Tools::isZero(const double var2test) {
  return (abs(var2test) < zero_threshold_ ? true : false);
}
