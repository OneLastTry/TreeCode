#include <boost/test/unit_test.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "../distributions/UniformDistribution.h"

using namespace std;
using namespace treecode;
using namespace Eigen;

#define D 2
#define N 1000

BOOST_AUTO_TEST_SUITE(Distributions)

BOOST_AUTO_TEST_CASE(TwoDUniformTest){
	boost::random::mt19937 rng;
	Configuration c(2, 0.0, 0.1);
	Vector2d min(-1,-1);
	Vector2d max(1,1);
	distribution::UniformDistribution<boost::random::mt19937> udist(c, min, max);
	for(int i=0;i<N;i++){
		Vector2d v = udist.getVector(rng);
		cerr << v[0] << "\t" << v[1] << endl;
	}
}


BOOST_AUTO_TEST_SUITE_END()
