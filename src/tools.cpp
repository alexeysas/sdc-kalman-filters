#include <iostream>
#include "tools.h"

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
	rmse << 0, 0, 0, 0;

	VectorXd residual;

	for (int i = 0; i < estimations.size(); ++i)
	{
		if (estimations[i].size() == 0)
		{
			throw std::invalid_argument("Estimation vector should not be zero size");
		}

		if (estimations[i].size() != ground_truth[i].size())
		{
			throw std::invalid_argument("Estimation vector should be same size as a ground truth");
		}
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		residual = (estimations[i] - ground_truth[i]);
		residual = residual.array() * residual.array();
		rmse = rmse + residual;
	}

	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state, float* c1_cache) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3, 4);

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1;

	// skip calculation if we already have them
	if (c1_cache == 0)
	{
		c1 = px * px + py * py;
	}
	else
	{
		c1 = *c1_cache;
	}

	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero - in kalman filter algorithm we should check this earlier to avoid exception
	if (fabs(c1) < 0.0001) {
		throw "Division by zero condition!";
	}

	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy*px) / c3, px*(px*vy - py*vx) / c3, px / c2, py / c2;

	return Hj;
}

VectorXd Tools::Hx(const VectorXd &x, float* c1_cache) {
	
	float c1;

	// skip calculation if we already have them
	if (c1_cache == 0)
	{
		c1 = x[0] * x[0] + x[1] * x[1];
	}
	else
	{
		c1 = *c1_cache;
	}

	//check division by zero - in kalman filter algorithm we should check this earlier to avoid exception
	if (fabs(c1) < 0.0001) {
		throw "Division by zero condition!";
	}
	
	VectorXd res = VectorXd(3);
	double p = sqrt(c1);

	double phi = atan2(x[1], x[0]);
	double v = (x[0] * x[2] + x[1] * x[3]) / p;
	res << p, phi, v;
	return res;
}


void Tools::NormalizeAngle(Eigen::VectorXd &input) {
	double res = input[1];

	while (res >= M_PI) {
		res -= 2 * M_PI;
	}

	while (res < -M_PI) {
		res += 2 * M_PI;
	}

	input[1] = res;
}


