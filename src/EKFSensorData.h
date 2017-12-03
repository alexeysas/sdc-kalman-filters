#ifndef EKFSensorData_H_
#define EKFSensorData_H_

#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

class EKFSensorData {
public:
	/**
	* Constructor.
	*/
	EKFSensorData();

	/**
	* Destructor.
	*/
	virtual ~EKFSensorData();

	Eigen::MatrixXd Hj_;
	Eigen::VectorXd x_;
private:
	float* radar_c1_cache ;
	Eigen::VectorXd RadarHFunc(VectorXd &x);
};

#endif /* FusionEKF_H_ */
