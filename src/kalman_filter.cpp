#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * Predict the state
   */

    // KF Prediction step
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Update the state by using Kalman Filter equations
   */

  // Temporary vectors and matrices
  VectorXd e;
  MatrixXd S;
  MatrixXd K;
  MatrixXd Ht;

  // KF Measurement update step
  Ht = H_.transpose();
  e = z - H_ * x_;
  S = H_ * P_ * Ht + R_;
  K = P_ * Ht * S.inverse();
  x_ += K * e;
  P_ -= K * H_ * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * Update the state by using Extended Kalman Filter equations
   *
   * The EFK Update will have two difference:
   *    1) H(x) will be used instead of H*x to calculate the error
   *    2) The H matrix used in the measurement update will be the Jacobian (Hj)
  **/

  // Temporary vectors and matrices
  VectorXd h_x(3);
  VectorXd e;
  MatrixXd S;
  MatrixXd K;
  MatrixXd Ht;

  // Calculate estimate of measurement from states
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  float p_hat = sqrt(px*px + py*py); // (p = rho)
  float phi_hat = atan2(py,px);
  float dp_hat = (px*vx + py*vy) / p_hat; // Note: dp_hat or p-dot is the dot produciton of x-y velocity and the unit vector in the direction of px-py

  h_x << p_hat, phi_hat, dp_hat;
  e = z - h_x;

  // Check that rho (angle) error has not rolled over -pi or pi
  e[1] = CheckAngleRollover(e[1]);

  // KF Measurement update step
  Ht = H_.transpose();

  S = H_ * P_ * Ht + R_;
  K = P_ * Ht * S.inverse();
  x_ += K * e;
  P_ -= K * H_ * P_;
}

float KalmanFilter::CheckAngleRollover(float angle)
{
  // Check that rho erro has not rolled over
  if (M_PI < angle)
  {
    while (M_PI < angle)
      angle -= 2*M_PI;
  }
  else if (angle < -M_PI)
  {
    while (angle < -M_PI)
      angle += 2*M_PI;
  }

  return angle;
}