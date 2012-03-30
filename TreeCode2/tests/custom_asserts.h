/*
 * custom_asserts.h
 *
 *  Created on: 30 Mar 2012
 *      Author: stefans
 */

#ifndef CUSTOM_ASSERTS_H_
#define CUSTOM_ASSERTS_H_

#include <Eigen/Dense>
void EIGEN_REQUIRE_CLOSE(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, double tol);

#endif /* CUSTOM_ASSERTS_H_ */
