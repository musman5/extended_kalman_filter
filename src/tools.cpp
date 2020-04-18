#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // Check if inputs are valid
  // The estimation vector size should not be zero
  // The estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "invalid estimation or ground_truth data" << endl;
    return rmse;
}

  //  Squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    
    // multiply coefficent-wise
    residual = residual.array() * residual.array();
    rmse = rmse + residual;
  }
  
  // Mean calculation
  rmse = rmse/ estimations.size();
  
  // Square root calculation
  rmse = rmse.array().sqrt();
  
  return rmse;
}
    
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  
  double px = x_state[0];
  double py = x_state[1];
  double vx = x_state[2];
  double vy = x_state[3];
  
  // pre-compute a set of terms to avoid repeated calculation
  double sqAndAdd = (px * px) + (py * py);
  double sqAndAddRoot = sqrt(sqAndAdd);
  double c = sqAndAddRoot * sqAndAdd;
  MatrixXd Jacobian = MatrixXd(3, 4);
  // check division by zero
  if (fabs(sqAndAdd) < 0.0001)
    return Jacobian;
  
  // compute the Jacobian matrix
  Jacobian << px / sqAndAddRoot , py / sqAndAddRoot, 0, 0,
  		-py / sqAndAdd, px/sqAndAdd, 0, 0,
  		py * (vx * py - vy*px) / c, px * (vy * px - 	vx * py) / c, px / sqAndAddRoot, py / sqAndAddRoot;
  return Jacobian;
}
