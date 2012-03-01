/*
 * InterpolatedEwaldSum.h
 *
 *  Created on: 5 Feb 2012
 *      Author: stefans
 */

#ifndef INTERPOLATEDEWALDSUM_H_
#define INTERPOLATEDEWALDSUM_H_

#include "EwaldForce.h"
#include "CoulombForce.h"
#include <boost/multi_array.hpp>

#include "EwaldForce.h"
#include "../bounds/BoundaryConditions.h"
#include "../Node.h"
#include "../Particle.h"
#include "../Configuration.h"
#include <cmath>

#define SQRT_PI 1.7724538509055159

namespace treecode {
namespace potentials {


template <class Vec, class Mat>
class EwaldNode{
public:
	EwaldNode(){
		t0 = 0;
		t1 = Vec::Zero();
		t2 = Mat::Zero();
	}

	EwaldNode operator+(const EwaldNode& n2){
		EwaldNode res;
		res.t0 = this->t0 + n2.t0;
		res.t1 = this->t1 + n2.t1;
		res.t2 = this->t2 + n2.t2;
		return res;
	}

	EwaldNode operator*(double x){
		EwaldNode res;
		res.t0 = this->t0*x;
		res.t1 = this->t1*x;
		res.t2 = this->t2*x;
		return res;
	}

	double t0;
	Vec t1;
	Mat t2;

};

template <class Vec, class Mat>
class InterpolatedEwaldSum : public Potential<Vec, Mat>{
	typedef boost::multi_array<EwaldNode<Vec,Mat>, 3> multi_arr;
	typedef typename multi_arr::index index;

public:
	InterpolatedEwaldSum(
			const Configuration<Vec>& conf,
			const BoundaryConditions<Vec>& bounds,
			unsigned int divisions,
			const EwaldForce<Vec, Mat>& ewald_pot,
			const CoulombForceThreeD<Vec, Mat>& coulomb_pot) :
			conf_(conf), bounds_(bounds), divisions_(divisions),
			ewald_pot_(ewald_pot), coulomb_pot_(coulomb_pot){

		length_per_div_ = bounds.getSize() / divisions;

		// Create a 3D array that is divisions x divisions x divisions
		field_ = new multi_arr(boost::extents[divisions + 1][divisions + 1][divisions + 1]);
	}

	void init() {
		/*
		 * Place a particle in the centre of the simulation region.
		 * Then, calculate the correction to the potential from
		 * the images, and store at that grid point.
		 */
		Vec centre_point = bounds_.getOrigin().array() + bounds_.getSize() / 2;

		#pragma omp parallel for
		for (index i = 0; i <= divisions_; i++) {
			for (index j = 0; j <= divisions_; j++) {
				for (index k = 0; k <= divisions_; k++) {
					Vec current_point = bounds_.getOrigin()
							+ Vec(i * length_per_div_, j * length_per_div_, k * length_per_div_);
					Vec disp_vec = current_point - centre_point;

	#ifndef __CDT_PARSER__		//Stop eclipse crapping out
					calculateNode(disp_vec, (*field_)[i][j][k]);
	#endif
				}
			}
		}
	}

	void outputField() {
		std::ofstream out("field.csv");
		for (int x = 0; x <= divisions_; x++) {
			for (int y = 0; y <= divisions_; y++) {
				out << x << "\t" << y << "\t" << (*field_)[x][y][0].t0 << std::endl;
			}
			out << std::endl;
		}
	}

	double getPotential(const Particle<Vec>& part, const Node<Vec, Mat>& node, Precision precision) const {
		Vec disp_vec = bounds_.getDisplacementVector(part.getPosition(), node.getCentreOfCharge());
		Vec centre_point = bounds_.getOrigin().array() + bounds_.getSize() / 2;
		Vec disp_vec_to_centre = disp_vec + centre_point;

		double potential = 0;
		EwaldNode<Vec,Mat> interpolated = interpolate(disp_vec_to_centre);
		potential -= node.getCharge() * interpolated.t0;
		potential -= node.getDipoleMoments().dot(interpolated.t1);
		potential -= 0.5 * (node.getQuadrupoleMoments().transpose() * interpolated.t2).trace();

		//Apply self-energy correction
		potential -= 2 * ewald_pot_.alpha_ / SQRT_PI * node.getCharge();
		//Apply scaling from unit system
		potential /= (3 * conf_.getDensity());
		//Add contribution from main cell
		potential += coulomb_pot_.getPotential(part, node, precision);
		return potential;
	}

	Vec getForce(const Particle<Vec>& part, const Node<Vec,Mat>& node, Precision precision) const {
		Vec disp_vec = bounds_.getDisplacementVector(part.getPosition(), node.getCentreOfCharge());
		Vec centre_point = bounds_.getOrigin().array() + bounds_.getSize() / 2;
		Vec disp_vec_to_centre = disp_vec + centre_point;

		EwaldNode<Vec,Mat> interpolated = interpolate(disp_vec_to_centre);
		Vec force = Vec::Zero();
		force += node.getCharge() * interpolated.t1;
		force += interpolated.t2 * node.getDipoleMoments();

		//Apply charge and scaling from unit system
		force = part.getCharge() * force / (3 * conf_.getDensity());
		//Add direct contribution
		force += coulomb_pot_.getForce(part, node, precision);
		return force;
	}

	EwaldNode<Vec,Mat> interpolate(const Vec& r) const {
		Vec dv = r - bounds_.getOrigin();

		double x = dv[0];
		double y = dv[1];
		double z = dv[2];

		int x0 = floor(x / length_per_div_);
		int y0 = floor(y / length_per_div_);
		int z0 = floor(z / length_per_div_);

		int x1 = ceil(x / length_per_div_);
		int y1 = ceil(y / length_per_div_);
		int z1 = ceil(z / length_per_div_);

		double xd = x/length_per_div_ - x0;
		double yd = y/length_per_div_ - y0;
		double zd = z/length_per_div_ - z0;

		EwaldNode<Vec,Mat> n000 = (*field_)[x0][y0][z0];
		EwaldNode<Vec,Mat> n001 = (*field_)[x0][y0][z1];
		EwaldNode<Vec,Mat> n010 = (*field_)[x0][y1][z0];
		EwaldNode<Vec,Mat> n011 = (*field_)[x0][y1][z1];
		EwaldNode<Vec,Mat> n100 = (*field_)[x1][y0][z0];
		EwaldNode<Vec,Mat> n101 = (*field_)[x1][y0][z1];
		EwaldNode<Vec,Mat> n110 = (*field_)[x1][y1][z0];
		EwaldNode<Vec,Mat> n111 = (*field_)[x1][y1][z1];

		EwaldNode<Vec,Mat> i1 = n000 * (1 - zd) + n001 * zd;
		EwaldNode<Vec,Mat> i2 = n010 * (1 - zd) + n011 * zd;
		EwaldNode<Vec,Mat> j1 = n100 * (1 - zd) + n101 * zd;
		EwaldNode<Vec,Mat> j2 = n110 * (1 - zd) + n111 * zd;

		EwaldNode<Vec,Mat> w1 = i1 * (1 - yd) + i2 * yd;
		EwaldNode<Vec,Mat> w2 = j1 * (1 - yd) + j2 * yd;

		EwaldNode<Vec,Mat> interpolated = w1 * (1 - xd) + w2 * xd;
		return interpolated;
	}

	void calculateNode(const Vec& disp_vec, EwaldNode<Vec,Mat>& node) {
		node.t0 = 0;
		node.t1 = Vec::Zero();
		node.t2 = Mat::Zero(3, 3);

		int real_space_iterations_ = ewald_pot_.real_space_iterations_;
		int fourier_space_iterations_ = ewald_pot_.fourier_space_iterations_;

		for (int i = -real_space_iterations_; i <= real_space_iterations_; i++) {
			for (int j = -real_space_iterations_; j <= real_space_iterations_; j++) {
				for (int k = -real_space_iterations_; k <= real_space_iterations_; k++) {
					//Displacement vector, plus periodic images
					Vec r_n(i * bounds_.getSize(), j * bounds_.getSize(), k * bounds_.getSize());
					r_n += disp_vec;
					double r_n_norm = r_n.norm();

	//				if(i >= -1 && i <= 1 && j >= -1 && j <= 1 && k >= -1 && k <= 1)
					if(i == 0 && j == 0 && k == 0)
						r_n_norm = sqrt(r_n_norm*r_n_norm + conf_.getForceSoftening()*conf_.getForceSoftening());

					//Then actually add contributions.
					node.t0 += ewald_pot_.a0_real(r_n, r_n_norm);
					node.t1 += ewald_pot_.a1_real(r_n, r_n_norm);
					node.t2 += ewald_pot_.a2_real(r_n, r_n_norm);
				}
			}
		}

		//This is exactly the same as above, but uses thje fourier space sum
		for (int i = -fourier_space_iterations_; i <= fourier_space_iterations_; i++) {
			for (int j = -fourier_space_iterations_; j <= fourier_space_iterations_; j++) {
				for (int k = -fourier_space_iterations_; k <= fourier_space_iterations_; k++) {
					if (i == 0 && j == 0 && k == 0)
						continue;
					Vec h(i, j, k);
					double h_norm = h.norm();

					node.t0 += ewald_pot_.a0_fourier(h, disp_vec, h_norm);
					node.t1 += ewald_pot_.a1_fourier(h, disp_vec, h_norm);
					node.t2 += ewald_pot_.a2_fourier(h, disp_vec, h_norm);
				}
			}
		}

		//Now, remove the contribution from the main cell
		double r = 1.0 / sqrt(disp_vec.squaredNorm() + conf_.getForceSoftening()*conf_.getForceSoftening());
		double r_3 = r * r * r;
		double r_5 = r_3 * r * r;
		node.t0 -= r;
		node.t1 -= disp_vec * r_3;
		node.t2 -= disp_vec * disp_vec.transpose() * 3 * r_5;
		node.t2 += Mat::Identity(3,3) * r_3;
	}

	~InterpolatedEwaldSum() {
		delete field_;
	}

private:
	multi_arr* field_;

	const Configuration<Vec>& conf_;
	const BoundaryConditions<Vec>& bounds_;
	unsigned int divisions_;
	double length_per_div_;

	const EwaldForce<Vec,Mat>& ewald_pot_;
	const CoulombForceThreeD<Vec,Mat>& coulomb_pot_;


};

} /* namespace output */
} /* namespace treecode */
#endif /* INTERPOLATEDEWALDSUM_H_ */
