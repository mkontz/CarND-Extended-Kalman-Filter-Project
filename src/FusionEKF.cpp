#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
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

  // measurement update matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  /**
   * Finish initializing the FusionEKF.
   *  1) Set the process
   *  2) State transition matrix
   *  3) Set size of preocess noises
   *  4) Set initial state covariance matrix
   */

  // the initial transition matrix F
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);

  // Initialize size of Process Noise Matrix
  ekf_.Q_ = MatrixXd(4, 4);

  // state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     */

    // first measurement
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       * Convert radar from polar to cartesian coordinates.
       * Initialize states.
      **/
      float rho = measurement_pack.raw_measurements_[1];
      float phi = measurement_pack.raw_measurements_[2];
      float px = rho * cos(phi);
      float py = rho * sin(phi);

      ekf_.x_ << px, py, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
       * Initialize states.
      **/
      float px = measurement_pack.raw_measurements_[1];
      float py = measurement_pack.raw_measurements_[2];
      ekf_.x_ << px, py, 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Calculate elapse time
   * Time is measured in seconds.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   */
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  /**
   * Update the process noise covariance matrix.
  // Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  **/
  {
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    float noise_ax = 9;
    float noise_ay = 9;

    ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  }

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
     * Radar update:
     *  1. Set R to R_radar
     *  2. Set H to Hj
     *  3. Construct z
     *  4. Cal UdpateEKF
    **/

    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    VectorXd z;
    z << measurement_pack.raw_measurements_[1], // p (rho)
         measurement_pack.raw_measurements_[2], // phi
         measurement_pack.raw_measurements_[3]; // d_rho/dt

    ekf_.UpdateEKF(z);
  }
  else
  {
    /**
     * Lase update:
     *  1. Set R to R_laser
     *  2. Set H to H_laser
     *  3. Construct z
     *  4. Cal Udpate
    **/

    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;

    VectorXd z;
    z << measurement_pack.raw_measurements_[1], // px
         measurement_pack.raw_measurements_[2]; // py

    ekf_.Update(z);
  }
}
