#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath> 

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

	H_laser_ << 1, 0, 0, 0,
		0, 1, 0, 0;

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

		// first measurement - initialize state
		ekf_.x_ = VectorXd(4);

		// if first measurement is radar - initialize from radar data by converting from polar coordinates 
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			double p = measurement_pack.raw_measurements_[0];
			double phi = measurement_pack.raw_measurements_[1];
			ekf_.x_ << p * cos(phi), p * sin(phi), 0, 0;

		}
		// if first measurement is laser - initialize from laser data taking  direct data
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
		}
		// just to ignore possible sensors we might dont know about 
		else {
			return;
		}

		// set initialization timestamp as a starting point 
		previous_timestamp_ = measurement_pack.timestamp_;

		// initialize transition matrix
		ekf_.F_ = MatrixXd(4, 4);

		ekf_.F_ << 1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;

		// initialize object covariance matrix
		ekf_.P_ = MatrixXd(4, 4);
		ekf_.P_ << 1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 100, 0,
			0, 0, 0, 100;

		// initalize process covariance matrix
		ekf_.Q_ = MatrixXd(4, 4);
		ekf_.Q_ << 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0;

		// done initializing, no need to predict or update
		is_initialized_ = true;
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


	  // Just for debug - to be able to switch one sensor of
	  //if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	  //{
	  //	  return;
	  // }


	  //compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;


	// Update F according to timestamp
	ekf_.F_ << 1, 0, dt, 0,
		0, 1, 0, dt,
		0, 0, 1, 0,
		0, 0, 0, 1;

	// Calculate process covariance matrix Q
	float dt2 = dt * dt;
	float dt3 = dt2 * dt;
	float dt4 = dt3 * dt;
	float qax2 = 9;
	float qay2 = 9;

	ekf_.Q_ << dt4 / 4 * qax2, 0, dt3 / 2 * qax2, 0,
		0, dt4 / 4 * qay2, 0, dt3 / 2 * qay2,
		dt3 / 2 * qax2, 0, dt2 * qax2, 0,
		0, dt3 / 2 * qay2, 0, dt2 * qay2;


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
		UpdateRadar(measurement_pack);

	}
	else {
		// Laser updates
		UpdateLaser(measurement_pack);
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}

void FusionEKF::UpdateRadar(const MeasurementPackage &measurement_pack) {
	float px = ekf_.x_(0);
	float py = ekf_.x_(1);
	float  radar_c1_cache = px * px + py * py;

	if (radar_c1_cache < 0.0001) {
		//set c1 to larger value if there is possibility of zero
		radar_c1_cache = 0.0001;
	}

	// Update measurement matrices
	// Calculate Hj jacobian
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_, &radar_c1_cache);
	ekf_.R_ = R_radar_;
	
	// calculate h(x) function
	ekf_.Hx_ = tools.Hx(ekf_.x_, &radar_c1_cache);
	
	// set angle normalization function as post-processing - this is needed to avoid radar specific code inside kalman filter class. 
	ekf_.EkfYPostProcessingFunc = Tools::NormalizeAngle;

	VectorXd z = VectorXd(3);
	z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
	ekf_.UpdateEKF(z);
}

void FusionEKF::UpdateLaser(const MeasurementPackage &measurement_pack) {
	//Apply lazer matrices
	ekf_.R_ = R_laser_;
	ekf_.H_ = H_laser_;

	//Get measurement vector
	VectorXd z = VectorXd(2);
	z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];

	ekf_.Update(z);
}