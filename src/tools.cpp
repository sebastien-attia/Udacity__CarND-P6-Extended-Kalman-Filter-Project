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
  TODO:
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
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	float dist2 = pow(px, 2)+pow(py, 2);
	float dist = pow(dist2, 0.5);
	if (dist2 < 0.0001) {
	    cout << "CalculateJacobian(): Division by 0" << endl;
      return Hj;
	}

	//compute the Jacobian matrix
  Hj << px/dist,  py/dist,  0, 0,
       -py/dist2, px/dist2, 0, 0,
        py*(vx*py-vy*px)/pow(dist2, 3), px*(vy*px - vx*py)/pow(dist2, 3), px/dist, py/dist;


	return Hj;
}

VectorXd Tools::ConvertToCartesian(const VectorXd& x_state) {
  float rho = x_state[0];
  float theta = x_state[1];
  float rho_dot = x_state[2];

  float cos_theta = cos(theta);
  float sin_theta = sin(theta);

  float x = cos_theta * rho;
  float y = sin_theta * rho;
  float vx = rho_dot * cos_theta;
  float vy = rho_dot * sin_theta;

  VectorXd result = VectorXd(4);
  result << x, y, vx, vy;
  return result;
}

VectorXd Tools::ConvertToPolar(const VectorXd& x_state) {
  float x = x_state[0];
  float y = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];

  float rho = sqrt(pow(x, 2) + pow(y, 2));
  float theta = atan2(y, x);
  cout << "Theta: " << theta *180 / 3.14169 << endl;
  float rho_dot = (x*vx + y*vy)/rho;

  VectorXd result = VectorXd(3);
  result << rho, theta, rho_dot;
  return result;
}
