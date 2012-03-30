/*
 * TreeTest.cpp
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */


#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <cmath>
#include <Eigen/Dense>

#include <distributions/UniformDistribution.h>
#include <distributions/MaxwellDistribution.h>
#include <distributions/SphericalDistribution.h>
#include <distributions/ConstDistribution.h>
#include <distributions/ConstantChargeDistribution.h>

#include "custom_asserts.h"

#include <iostream>

#ifndef STAT_TOLERANCE
#define STAT_TOLERANCE 1.5
#endif
#ifndef TOLERANCE
#define TOLERANCE 1e-10
#endif

//Just to stop the CDT parser deciding EVERYTHING needs a beautiful red underline
#ifdef __CDT_PARSER__
#define BOOST_REQUIRE_GE(a,b)
#define BOOST_REQUIRE_LE(a,b)
#define BOOST_REQUIRE_CLOSE(a,b,c)
#define BOOST_REQUIRE_SMALL(a,b)
#endif

using namespace boost::accumulators;
using namespace treecode;
using boost::random::mt19937;
typedef Eigen::Vector2d V2d;
typedef Eigen::Vector3d V3d;

//Record min, max, mean, median, variance, skew, kurtosis,
typedef accumulator_set<double, features<
			tag::min,
			tag::max,
			tag::mean,
			tag::median,
			tag::variance,
			tag::skewness,
			tag::kurtosis> > acc_set;

BOOST_AUTO_TEST_SUITE(Distributions)

BOOST_AUTO_TEST_CASE(UDist){
	acc_set acc_x;
	acc_set acc_y;

	V2d min_vec(0.5,0.5);
	V2d max_vec(1,1);
	mt19937 rng;
	distribution::UniformDistribution<mt19937> udist(rng, min_vec, max_vec);
	//Can't be bothered defining it to work on Vec, just use each component.
	for(int i=0;i<50000;i++){
		V2d v = udist.getVector();
		acc_x(v[0]);
		acc_y(v[1]);
	}

	//Check min/max are within range
	BOOST_REQUIRE_GE( min(acc_x), min_vec[0] );
	BOOST_REQUIRE_GE( min(acc_y), min_vec[1] );
	BOOST_REQUIRE_LE( max(acc_x), max_vec[0] );
	BOOST_REQUIRE_LE( max(acc_x), max_vec[0] );

	//Check mean
	BOOST_REQUIRE_CLOSE( mean(acc_x), 0.5 * (min_vec[0] + max_vec[0]), STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( mean(acc_y), 0.5 * (min_vec[1] + max_vec[1]), STAT_TOLERANCE);
	//Median
	BOOST_REQUIRE_CLOSE( median(acc_x), 0.5 * (min_vec[0] + max_vec[0]), STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( median(acc_y), 0.5 * (min_vec[1] + max_vec[1]), STAT_TOLERANCE);
	//Variance
	BOOST_REQUIRE_CLOSE( variance(acc_x), 1.0/12 * (max_vec[0]-min_vec[0])*(max_vec[0]-min_vec[0]), STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( variance(acc_y), 1.0/12 * (max_vec[1]-min_vec[1])*(max_vec[1]-min_vec[1]), STAT_TOLERANCE);
	//Skew
	BOOST_REQUIRE_SMALL( skewness(acc_x), STAT_TOLERANCE);
	BOOST_REQUIRE_SMALL( skewness(acc_y), STAT_TOLERANCE);
	//Kurtosis
	BOOST_REQUIRE_CLOSE( kurtosis(acc_x), -6.0/5, STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( kurtosis(acc_y), -6.0/5, STAT_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(SphDist){
	acc_set acc_r;
	acc_set acc_x, acc_y, acc_z;

	V3d origin(0.5,0.5,0.5);
	double radius = 1;

	mt19937 rng;
	distribution::SphericalDistribution<mt19937> sph_dist(rng, 3, origin, radius);
	//Can't be bothered defining it to work on Vec, just use each component.
	for(int i=0;i<50000;i++){
		V3d v = sph_dist.getVector();
		acc_r( (v-origin).norm() );
		acc_x(v[0]);
		acc_y(v[1]);
		acc_z(v[1]);
	}

	//Check min/max are within range
	BOOST_REQUIRE_SMALL( min(acc_r),  STAT_TOLERANCE);
	BOOST_REQUIRE_LE( max(acc_r), radius );

	//Check mean
	BOOST_REQUIRE_CLOSE( mean(acc_x), origin[0], STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( mean(acc_y), origin[1], STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( mean(acc_z), origin[2], STAT_TOLERANCE);
	//Median
	BOOST_REQUIRE_CLOSE( median(acc_x), origin[0], STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( median(acc_y), origin[1], STAT_TOLERANCE);
	BOOST_REQUIRE_CLOSE( median(acc_z), origin[2], STAT_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(MaxwellDist){
	acc_set acc_v;

	V2d min_vec(0.5,0.5);
	V2d max_vec(1,1);
	mt19937 rng;
	distribution::MaxwellDistribution<mt19937> mdist(rng, 1, 1, 3);
	//Can't be bothered defining it to work on Vec, just use each component.
	for(int i=0;i<1000000;i++){
		V3d v = mdist.getVector();
		acc_v(v.norm());
	}

	//Check greater than zero
	BOOST_REQUIRE_GE(min(acc_v), 0);

	//Check mean
	BOOST_REQUIRE_CLOSE( mean(acc_v), 2 * sqrt(2.0 / M_PI), STAT_TOLERANCE);
	//Variance
	BOOST_REQUIRE_CLOSE( variance(acc_v), (3.0*M_PI-8.0)/M_PI, STAT_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(ConstDist){
	V2d const_v(M_PI, M_PI);
	distribution::ConstDistribution cdist(const_v);
	EIGEN_REQUIRE_CLOSE(const_v, cdist.getVector(), TOLERANCE);
}

BOOST_AUTO_TEST_CASE(ConstChargeDist){
	int charge = 1;
	distribution::ConstantChargeDistribution cdist(charge);
	BOOST_REQUIRE_EQUAL(charge, cdist.getCharge());
}

BOOST_AUTO_TEST_CASE(SinDist){
	//TODO
}

BOOST_AUTO_TEST_SUITE_END();
