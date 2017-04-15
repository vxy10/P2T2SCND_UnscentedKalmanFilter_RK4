#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// function for state Prediction
VectorXd StatePredict(VectorXd x, double dt){

    int n_x = 5;

    double px =  x(0);
    double py =  x(1);
    double  v =  x(2);
    double psi = x(3);
    double dpsi =x(4);
    double va =  x(5);
    double vp =  x(6);

    double dt2 = dt*dt;

    VectorXd x_pl = VectorXd(n_x);
    VectorXd fx = VectorXd(n_x);
    fx.setZero();

    if (fabs(dpsi)<0.001){

        fx << v*cos(psi)*dt + 1.0/2.0*dt2*cos(psi)*va,
              v*sin(psi)*dt + 1.0/2.0*dt2*sin(psi)*va,
              dt*va,
              1/2*dt2*vp,
              dt*vp;

    }else{
        fx << v/dpsi*(sin(psi+dpsi*dt)-sin(psi)) + 1.0/2.0*dt2*cos(psi)*va,
              v/dpsi*(-cos(psi+dpsi*dt)+cos(psi)) + 1.0/2.0*dt2*sin(psi)*va,
              dt*va,
              dpsi*dt+1/2*dt2*vp,
              dt*vp;

    }
    x_pl = x.head(n_x) + fx;
    return x_pl;
}


// Radar measurement function
VectorXd MeasurementFRadar(VectorXd x){
    double  px = x(0);
    double  py = x(1);
    double   v = x(2);
    double  psi = x(3);
    double dpsi = x(4);

    VectorXd Z = VectorXd(3);

    if (sqrt(pow(px*px + py*py,2))>0.001){

      Z(0) = sqrt(px*px+py*py);
      Z(1) = atan2(py,px);
      Z(2) = (px*v*cos(psi) + py*v*sin(psi))/Z(0);
      return Z;
  } else {
    Z(0) = sqrt(0.001);
    Z(1) = atan2(0.001,0.001);
    Z(2) = (px*v*cos(psi) + py*v*sin(psi))/Z(0);
    return Z;
  }
}

// LIDAR measurement function
VectorXd MeasurementFLidar(VectorXd x){
    double  px = x(0);
    double  py = x(1);
    VectorXd Z = VectorXd(2);
    Z(0) = px;
    Z(1) = py;
    return Z;
}







/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // define dimensions
  n_x_ = 5;
  n_aug_ = n_x_ + 2;

  lambda_ = 3-n_aug_;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  P_ << 1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1;


  Xsig_pred_ = MatrixXd(5, 2*7+1);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "UKF: " << endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_(0);
      float psi =  meas_package.raw_measurements_(1);
      float drho =  meas_package.raw_measurements_(2);

      x_(0) = rho*cos(psi);
      x_(1) = rho*sin(psi);
      x_(2) = drho; // Approximate value as dth is not known
      x_(3) = psi; // Approximate value as dth is not known
      x_(4) = 0; // Approximate value as dth is not known
      if (fabs(rho)>0.001){
        is_initialized_ = true;
      }
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) =meas_package.raw_measurements_(0);
      x_(1) =meas_package.raw_measurements_(1);
      x_(2) = 0; // Approximate value of 0
      x_(3) = 0; // Approximate value of 0
      x_(4) = 0; // Approximate value of 0

      if (sqrt(pow(x_(0),2)+pow(x_(1),2))>0.001){
        is_initialized_ = true;
      }

    }

    // done initializing, no need to predict or update

    time_us_ = meas_package.timestamp_;
    return;

  }

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;


  //cout << "dt = " << dt <<endl<<endl;

  //Use small dt to allow for turn effect
  const double diff_t = 0.1;

  while (dt > diff_t){
		    Prediction(diff_t);
		    dt -= diff_t;
	}
  Prediction(dt); // update states only if dt is above 0.001



  //cout << "x_ = " << x_<<endl<<endl;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    if (fabs(meas_package.raw_measurements_(0))>0.001){
      UpdateRadar(meas_package);
    }

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      double l_px = meas_package.raw_measurements_(0);
      double l_py = meas_package.raw_measurements_(1);
      if (sqrt(pow(l_px,2)+pow(l_py,2))>0.001){
        UpdateLidar(meas_package);
      }
  }


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

// Generate sigma points
  //initialization of matrices for sigma point calculations
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  MatrixXd A_aug = MatrixXd(n_aug_,n_aug_);
  MatrixXd Ones_nAug = MatrixXd(1, n_aug_);
  Ones_nAug.setOnes();
  MatrixXd Ones_nA = MatrixXd(1, 2*n_aug_+1);
  Ones_nA.setOnes();
  //calculate AUGMENTED points ...
  P_aug.setZero();
  P_aug.topLeftCorner( n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_; // 0 based indexing
  P_aug(6,6) = std_yawdd_*std_yawdd_; // 0 based indexing

  A_aug = P_aug.llt().matrixL();
  x_aug.setZero();
  x_aug.head(n_x_) << x_;

  // augmented sigma points

  Xsig_aug << x_aug,
              x_aug*Ones_nAug+sqrt(lambda_+n_aug_)*A_aug,
              x_aug*Ones_nAug-sqrt(lambda_+n_aug_)*A_aug; // generate sigma pts




  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred.col(0) = StatePredict(Xsig_aug.col(0),delta_t);
  for (int i=1;i<=n_aug_;i++){
    Xsig_pred.col(i) = StatePredict(Xsig_aug.col(i),delta_t);
    Xsig_pred.col(i+n_aug_) = StatePredict(Xsig_aug.col(i+n_aug_),delta_t);
  }

  VectorXd weights = VectorXd(2*n_aug_+1);
  MatrixXd Wts_diag = MatrixXd(2*n_aug_+1,2*n_aug_+1);
  //set weights
  weights(0) = lambda_/(lambda_+n_aug_);
  weights.tail(2*n_aug_).fill(1/2./(lambda_+n_aug_));
  Wts_diag  = MatrixXd(weights.asDiagonal());
  //predict state mean
  x_ = Xsig_pred*weights;
  // predict covariance
  P_ =(Xsig_pred-x_*Ones_nA)*Wts_diag*(Xsig_pred-x_*Ones_nA).transpose();
  Xsig_pred_ = Xsig_pred;



}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z_ = 2;

  VectorXd weights = VectorXd(2*n_aug_+1);
  MatrixXd Wts_diag = MatrixXd(2*n_aug_+1,2*n_aug_+1);
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  MatrixXd S = MatrixXd(n_z_,n_z_);
  VectorXd z_pred = VectorXd(n_z_);
  VectorXd z_true = VectorXd(n_z_);
  MatrixXd Ones_A = MatrixXd(1,2*n_aug_+1);
  Ones_A.setOnes();

  for (int i=0;i<2*n_aug_+1;i++){
      Zsig.col(i) = MeasurementFLidar(Xsig_pred_.col(i));
  }
  weights(0) = lambda_/(lambda_+n_aug_);
  weights.tail(2*n_aug_).fill(1/2./(lambda_+n_aug_));
  //predict state mean
  z_pred = Zsig*weights;
  //calculate measurement covariance matrix S
  S =(Zsig-z_pred*Ones_A)*MatrixXd(weights.asDiagonal())*(Zsig-z_pred*Ones_A).transpose();
  VectorXd R_d = VectorXd(n_z_);
  R_d << std_laspx_*std_laspx_,
        std_laspy_*std_laspy_;
  S = S + MatrixXd(R_d.asDiagonal());

  z_true = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  MatrixXd Z_diff = MatrixXd(n_z_,2*n_aug_+1);
  //calculate cross correlation matrix
  MatrixXd K = MatrixXd(n_x_,n_z_);
  Z_diff = (Zsig-z_pred*Ones_A);

  Tc = (Xsig_pred_-x_*Ones_A)*MatrixXd(weights.asDiagonal())*Z_diff.transpose();
  //calculate Kalman gain K;
  K = Tc*S.inverse();
  //update state mean and covariance matrix
  x_ = x_ + K *(z_true - z_pred);
  P_ = P_ - K*S*K.transpose();
  NIS_laser_ = (z_true - z_pred).transpose()*S*(z_true - z_pred);




}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z_ = 3;

  VectorXd weights = VectorXd(2*n_aug_+1);
  MatrixXd Wts_diag = MatrixXd(2*n_aug_+1,2*n_aug_+1);
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  MatrixXd S = MatrixXd(n_z_,n_z_);
  VectorXd z_pred = VectorXd(n_z_);
  VectorXd z_true = VectorXd(n_z_);
  MatrixXd Ones_A = MatrixXd(1,2*n_aug_+1);
  Ones_A.setOnes();

  for (int i=0;i<2*n_aug_+1;i++){
      Zsig.col(i) = MeasurementFRadar(Xsig_pred_.col(i));
  }
  weights(0) = lambda_/(lambda_+n_aug_);
  weights.tail(2*n_aug_).fill(1/2./(lambda_+n_aug_));
  //predict state mean
  z_pred = Zsig*weights;
  //calculate measurement covariance matrix S
  S =(Zsig-z_pred*Ones_A)*MatrixXd(weights.asDiagonal())*(Zsig-z_pred*Ones_A).transpose();
  VectorXd R_d = VectorXd(n_z_);
  R_d << std_radr_*std_radr_,
        std_radphi_*std_radphi_,
        std_radrd_*std_radrd_;
  S = S + MatrixXd(R_d.asDiagonal());
  //cout<<"S"<<S<<endl;
  z_true = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  MatrixXd Z_diff = MatrixXd(n_z_,2*n_aug_+1);
  //calculate cross correlation matrix
  MatrixXd K = MatrixXd(n_x_,n_z_);
  Z_diff = (Zsig-z_pred*Ones_A);

  for (int i=0;i<n_z_;i++){
      Z_diff(i,2) = atan2(sin(Z_diff(i,2)),cos(Z_diff(i,2)));

  }

  Tc = (Xsig_pred_-x_*Ones_A)*MatrixXd(weights.asDiagonal())*Z_diff.transpose();
  //calculate Kalman gain K;
  K = Tc*S.inverse();
  //update state mean and covariance matrix
  x_ = x_ + K *(z_true - z_pred);

  P_ = P_ - K*S*K.transpose();
  NIS_radar_ = (z_true - z_pred).transpose()*S*(z_true - z_pred);

}
