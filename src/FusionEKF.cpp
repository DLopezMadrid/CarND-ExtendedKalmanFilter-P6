#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

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
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // H is the measurement matrix. H is a simple matrix in the case of the lidar
  // In the case of the radar the measurement noise is a function of x -> h(x)
  // since KF only work on linear systems, it needs to be approximated using a
  // taylor series expansion in the case of the radar (need to calculate the jacobian)

  H_laser_ << 1,0,0,0,
              0,1,0,0;

  // we can initialize the jacobian like this, it will be later updated for each actual position
  Hj_ << 1,1,0,0,
         1,1,0,0,
         1,1,1,1;

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
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Initialize the object covariance matrix P

    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1,0,   0,   0,
               0,1,   0,   0,
               0,0,1000,   0,
               0,0,   0,1000;

    // Initialize the state transition matrix F. This depends in dt and will be updated later for each measurement
    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1,0,1,0,
               0,1,0,1,
               0,0,1,0,
               0,0,0,1;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);

      //initialize state of KF
      ekf_.x_ << rho * cos(phi), rho * sin (phi), rho_dot * cos(phi), rho_dot * sin (phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1),0,0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // calculate dt in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0 ;
  previous_timestamp_ = measurement_pack.timestamp_;


  // update F for the current measurement
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, dt,  0,
             0, 1,  0, dt,
             0, 0,  1,  0,
             0, 0,  0,  1;

  noise_ax = 9;
  noise_ay = 9;

  float dt_2 = pow(dt, 2.0);
  float dt_3 = pow(dt, 3.0);
  float dt_4 = pow(dt, 4.0);

  // set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << dt_4 / 4 * noise_ax,                   0, dt_3 / 2 * noise_ax,                   0,
                               0, dt_4 / 4 * noise_ay,                   0, dt_3 / 2 * noise_ay,
             dt_3 / 2 * noise_ax,                   0,       dt_2*noise_ax,                   0,
                               0, dt_3 / 2 * noise_ay,                   0,       dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
