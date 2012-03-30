#include "custom_asserts.h"
#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>

void EIGEN_REQUIRE_CLOSE(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, double tol){
	BOOST_REQUIRE_EQUAL(a.rows(), b.rows());
	BOOST_REQUIRE_EQUAL(a.cols(), b.cols());

	for(int i=0;i<a.rows();i++){
		for(int j=0;j<a.cols();j++){
			BOOST_REQUIRE_CLOSE(a(i,j), b(i,j), tol);
		}
	}
}
