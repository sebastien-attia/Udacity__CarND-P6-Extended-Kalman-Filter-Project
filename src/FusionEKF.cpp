#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
			        0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd::Identity(4, 4);

    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_ = tools.ConvertToCartesian(measurement_pack.raw_measurements_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  {
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1,  0, dt, 0,
               0,  1, 0, dt,
               0,  0, 1, 0,
               0,  0, 0, 1;
  }
  {
    float noise_ax = 9, noise_ay = 9;
    float nx4 = pow(dt, 4) * noise_ax/4, nx3 = pow(dt, 3) * noise_ax/2, nx2 = pow(dt, 2) * noise_ax;
	  float ny4 = pow(dt, 4) * noise_ay/4, ny3 = pow(dt, 3) * noise_ay/2, ny2 = pow(dt, 2) * noise_ay;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << nx4,   0, nx3,   0,
		             0, ny4,   0, ny3,
	             nx3,   0, nx2,   0,
		             0, ny3,   0, ny2;
  }

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  VectorXd raw_measurements = measurement_pack.raw_measurements_;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    VectorXd cartesian = tools.ConvertToCartesian(raw_measurements);
    Hj_ = tools.CalculateJacobian(cartesian);

    ekf_.Init(Hj_, R_radar_);
    ekf_.UpdateEKF(raw_measurements);
  } else {
    // Laser updates
    ekf_.Init(H_laser_, R_laser_);
    ekf_.Update(raw_measurements);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
