#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::endl;
using std::cout;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  // Init rmse vector
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  if (estimations.size() != ground_truth.size())
  {
    cout << "Lengths estimations and ground truth don't match." << endl;
  }
  else if (0 == estimations.size())
  {
     cout << "Lengths estimations and ground truth are zero." << endl;
  }
  else
  {
	// accumulate squared residuals
	VectorXd res;

	for (size_t i=0; i < estimations.size(); i++)
	{
	  res = estimations.at(i) - ground_truth.at(i);
	  rmse += res*res;
	}

	// Calc mean & sqrt
	rmse = rmse / estimations.size();
	rmse.array().sqrt();
  }

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
  /**
   * Recover state parameters
  **/
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  /*
   * Construct Hj to return
  **/
  MatrixXd Hj(3,4);

  float p2 = px*px + py*py;
  if (0.0001 < p2)
  {
    float p = sqrt(p2);
    float p3 = p2*p;
    Hj << px/p, py/p, 0, 0,
         -py/p2, px/p2, 0, 0,
          py*(vx*py - vy*px)/p3, px*(vy*px - vx*py)/p3, px/p, py/p;
  }
  else
  {
      cout << "Need to divide by zero" << endl;
  }

  return Hj;
}
