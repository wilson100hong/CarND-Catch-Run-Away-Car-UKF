#include "ukf.h"
#include "Eigen/Dense"

#include <iostream>

#define EPS 1e-3

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

double NormalizeAngle(double rad) {
  while (rad > M_PI) rad -= 2*M_PI;
  while (rad < -M_PI) rad += 2*M_PI;
  return rad;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // If this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // If this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;

  // Initialize state vector
  x_ = VectorXd::Zero(n_x_);
  // Initalize measurement noise covariance matrix
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  R_radar_ = MatrixXd(3, 3); 
  R_radar_ <<  std_radr_ * std_radr_, 0, 0,
               0, std_radphi_ * std_radphi_, 0,
               0, 0, std_radrd_ * std_radrd_;

  // Initialize covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 
    std_laspx_*std_laspx_, 0, 0, 0, 0,
    0, std_laspy_*std_laspy_, 0, 0, 0,
    0, 0, 1, 0, 0, 
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1;

  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_+1);

  // Initialize sigma weights
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1;i<2*n_aug_+1;++i) {
    weights_(i) = 0.5/(lambda_+n_aug_);
  }

  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  auto sensor_type = meas_package.sensor_type_;
  if (!is_initialized_) {
    if (sensor_type == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      // and initialize state.
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx*vx + vy*vy);
      x_ << px, py, v, 0, 0;
    } else if (sensor_type == MeasurementPackage::LASER) {
      // Initialize state from lidar. 
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // Predict
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  // Update
  if (sensor_type == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (sensor_type == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  /*
   Estimate the object's location. 
   Modify the state vector, x_. Predict sigma points, the state, 
   and the state covariance matrix.
   */
   // Generate augmented sigma points: Xsig_aug
   VectorXd x_aug = VectorXd(n_aug_);  
   x_aug.head(n_x_) = x_;
   x_aug(5) = 0;
   x_aug(6) = 0;

   // Create augmented state covariance.
   MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
   P_aug.topLeftCorner(n_x_, n_x_) = P_;
   P_aug(5, 5) = std_a_ * std_a_;
   P_aug(6, 6) = std_yawdd_* std_yawdd_;

   // Create augmented sigma points.
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
   MatrixXd A = P_aug.llt().matrixL();

   Xsig_aug.col(0) = x_aug;
   double lns = sqrt(lambda_ + n_aug_);
   for (int i=0;i<n_aug_;++i) {
     Xsig_aug.col(i+1) = x_aug + lns * A.col(i);
     Xsig_aug.col(i+1+n_aug_) = x_aug - lns * A.col(i);
   }

   // Predict augmented sigma points
   Xsig_pred_.fill(0.0);
   for (int i=0;i<2*n_aug_+1;++i) {
     // Extract values for better readability
     double p_x = Xsig_aug(0,i);
     double p_y = Xsig_aug(1,i);
     double v = Xsig_aug(2,i);
     double yaw = Xsig_aug(3,i);
     double yawd = Xsig_aug(4,i);
     double nu_a = Xsig_aug(5,i);
     double nu_yawdd = Xsig_aug(6,i);

     // Predicted state values
     double px_p, py_p;

     // Avoid division by zero
     if (fabs(yawd) > EPS) {
       px_p = p_x + v/yawd * (sin (yaw + yawd*delta_t) - sin(yaw));
       py_p = p_y + v/yawd * (cos(yaw) - cos(yaw+yawd*delta_t) );
     } else {
       px_p = p_x + v*delta_t*cos(yaw);
       py_p = p_y + v*delta_t*sin(yaw);
     }

     double v_p = v;
     double yaw_p = yaw + yawd*delta_t;
     double yawd_p = yawd;

     // Add noise
     px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
     py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
     v_p = v_p + nu_a*delta_t;
     yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
     yawd_p = yawd_p + nu_yawdd*delta_t;

     // Write predicted sigma point into right column
     Xsig_pred_(0,i) = px_p;
     Xsig_pred_(1,i) = py_p;
     Xsig_pred_(2,i) = v_p;
     Xsig_pred_(3,i) = yaw_p;
     Xsig_pred_(4,i) = yawd_p;
   }

   // Predict state mean
   x_ = Xsig_pred_ * weights_;

   // Predict state covariance
   P_.fill(0);
   for (int i=0;i<2*n_aug_+1;++i) {
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     // Need to normalize angle
     x_diff(3) = NormalizeAngle(x_diff(3));
     P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
   }
}

// TODO: refactor UpdateLiar and UpdateRadar to share code.
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /*
   Use lidar data to update the belief 
   about the object's position. Modify the state vector, x_, and 
   covariance, P_.
   You can also calculate the lidar NIS, if desired.
   */
  // TODO: use XSig_pred
  int n_z = 2;

  MatrixXd H(n_z, n_x_);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  // Measurement residual
  VectorXd y = meas_package.raw_measurements_ - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R_laser_;
  MatrixXd K = P_* H.transpose() * S.inverse();

  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /*
   radar data to update the belief 
   about the object's position. Modify the state vector, x_, and 
   covariance, P_.
   You can also calculate the radar NIS, if desired.
   */
  int n_z = 3;
  // Sigma points in measuresument space.
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  // Transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // Extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v_x = cos(yaw)*v;
    double v_y = sin(yaw)*v;

    // Measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        // r
    Zsig(1,i) = atan2(p_y, p_x);                                // phi
    Zsig(2,i) = (p_x*v_x + p_y*v_y) / sqrt(p_x*p_x + p_y*p_y);  // r_dot
  }

  // Predicted measurement mean
  VectorXd z_pred = Zsig * weights_;

  // Innovation covariance
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Adds measurement noise
  S = S + R_radar_;

  // Cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // Measurement residual
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  y(1) = NormalizeAngle(y(1));

  // Update state mean and covariance matrix
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();
}
