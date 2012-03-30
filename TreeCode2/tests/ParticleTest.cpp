#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <Eigen/Dense>

#include <Particle.h>
#include "custom_asserts.h"

#ifndef TOLERANCE
#define TOLERANCE 1e-10
#endif

BOOST_AUTO_TEST_SUITE(Particle);

BOOST_AUTO_TEST_CASE(InterfaceTest){
	using treecode::Particle;
	Eigen::Vector3d position(1,2,3);
	Eigen::Vector3d velocity(1,2,3);

	Particle<3> p(1, 1, position, velocity);
	//Test getters
	EIGEN_REQUIRE_CLOSE(p.getPosition(), position, TOLERANCE);
	EIGEN_REQUIRE_CLOSE(p.getVelocity(), velocity, TOLERANCE);
	BOOST_REQUIRE_EQUAL(1, p.getCharge());
	BOOST_REQUIRE_CLOSE(1.0, p.getMass(), TOLERANCE);
	//Test setters
	position = Eigen::Vector3d(3,2,1);
	velocity = Eigen::Vector3d(3,2,1);
	p.setPosition(position);
	p.setVelocity(velocity);
	EIGEN_REQUIRE_CLOSE(p.getPosition(), position, TOLERANCE);
	EIGEN_REQUIRE_CLOSE(p.getVelocity(), velocity, TOLERANCE);
	BOOST_REQUIRE_EQUAL(1, p.getCharge());
	BOOST_REQUIRE_CLOSE(1.0, p.getMass(), TOLERANCE);
	//Check update commands
	p.setPosition(Eigen::Vector3d::Zero());
	p.setVelocity(Eigen::Vector3d::Zero());
	p.updatePosition(Eigen::Vector3d(1,1,1));
	p.updateVelocity(Eigen::Vector3d(1,1,1));
	EIGEN_REQUIRE_CLOSE(p.getPosition(), Eigen::Vector3d(1,1,1), TOLERANCE);
	EIGEN_REQUIRE_CLOSE(p.getVelocity(), Eigen::Vector3d(1,1,1), TOLERANCE);
}

BOOST_AUTO_TEST_SUITE_END();
