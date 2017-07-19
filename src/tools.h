#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  VectorXd ConvertToCartesian(const VectorXd& x_state);
  VectorXd ConvertToPolar(const VectorXd& x_state);

  double rangeAngle(const double theta);

  bool isZero(const double var2test);

private:
  const double zero_threshold_ = 0.000001;
};

#endif /* TOOLS_H_ */
