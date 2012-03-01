/*
 * EwaldForce.h
 *
 *  Created on: 20 Jan 2012
 *      Author: stefans
 */

#ifndef EWALDFORCE_H_
#define EWALDFORCE_H_

#include "Potential.h"
#include "../bounds/BoundaryConditions.h"
#include "../Particle.h"
#include "../Configuration.h"
#include "../Node.h"
//#include "InterpolatedEwaldSum.h"

#define SQRT_PI 1.7724538509055159
#define PI_SQUARED 9.869604401089358

namespace treecode {
namespace potentials {

template <class Vec, class Mat> class InterpolatedEwaldSum;

template <class Vec, class Mat>
class EwaldForce : public Potential<Vec, Mat>{
	friend class InterpolatedEwaldSum<Vec,Mat>;

public:

	/**
	 * @class EwaldForce "potentials/EwaldForce.h
	 * @brief A potential using Ewald summations to approximate an infinitely periodic system.
	 *
	 * This class defines a potential used to approximate an infinitely periodic system.
	 * For information, see <em>Many Body Tree Methods in Physics</em>. Make sure to
	 * <em>check the errata</em>!
	 */

	/**
	 * @brief Construct a new EwaldForce object.
	 * @param conf		Global configuration. Mainly used for force softening param.
	 * @param bounds	Boundary conditions. These should be <em>periodic</em> boundaries.
	 * @param alpha		Coupling parameter. Usually @f$ \alpha \approx 2 / L @f$, where @f$ L @f$ is system size.
	 * @param real_space_iterations		Number of iterations to perform in real space.
	 * @param fourier_space_iterations	Number of iterations to perform in Fourier space.
	 */
	EwaldForce(
			const Configuration<Vec>& conf,
			const BoundaryConditions<Vec>& bounds,
			double alpha,
			int real_space_iterations, int fourier_space_iterations) :
			conf_(conf),
			bounds_(bounds), alpha_(alpha),
			real_space_iterations_(real_space_iterations),
			fourier_space_iterations_(fourier_space_iterations) {

	}

	~EwaldForce() {
		// TODO Auto-generated destructor stub
	}

	/**
	 * @brief Get the potential acting on a particle.
	 *
	 * This works up to quadrupole moment, but the force calculation
	 * only works up to dipole moment, so be careful.
	 *
	 * The potential, @f$ \Phi_p @f$, is given by the following expression:
	 * @f[
	 * \Phi_p = M a_0 + \vec{D} \cdot \vec{a}_1 + 0.5 \operatorname{Tr}(\mathbf{Q}^T \mathbf{a}_2),
	 * @f]
	 * where @f$ M @f$ is the monopole moment of a node,
	 * @f$ \vec{D} @f$ is the dipole moment of a node and
	 * @f$ \mathbf{Q} @f$ is the quadrupole moment of a node.
	 *
	 * The various $a$s in the above equation are given by other functions in this.
	 *
	 * @see a2_real(), a2_fourier(), a1_real(), a1_fourier(), a0_real(), a0_fourier()
	 *
	 * @param part			Particle to calculate the potential at.
	 * @param node			Node acting on the particle.
	 * @param precision		Precision to calculate to.
	 * @return	Potential on particle due to node.
	 *
	 * TODO: Pay attention to precision paremeter
	 */
	double getPotential(const Particle<Vec>& part, const Node<Vec, Mat>& node, Precision precision) const{
		double potential = 0;

		//We'll need this later.
		Mat Q_trans = node.getQuadrupoleMoments().transpose();

		//Real space sums. Everything is linear, so we can do fourier and real
		//sums separately.
		for (int i = -real_space_iterations_; i <= real_space_iterations_; i++) {
			for (int j = -real_space_iterations_; j <= real_space_iterations_; j++) {
				for (int k = -real_space_iterations_; k <= real_space_iterations_; k++) {
					//Displacement vector, plus periodic images
					Vec r_n(i*bounds_.getSize(), j*bounds_.getSize(), k*bounds_.getSize());
					r_n += (part.getPosition() - node.getCentreOfCharge());
					double r_n_norm = r_n.norm();

					//Apply force softening to own cell and nearest neighbours.
					//This is because we can be close to particles in the neighbours too
					//I'm not entirely sure this is right, but it seems so.
					if(i >= -1 && i <= 1 && j >= -1 && j <= 1 && k >= -1 && k <= 1)
						r_n_norm = sqrt(r_n_norm*r_n_norm + conf_.getForceSoftening()*conf_.getForceSoftening());

					//Then actually add contributions.
					potential += node.getCharge() * a0_real(r_n, r_n_norm);
					potential += node.getDipoleMoments().dot(a1_real(r_n, r_n_norm));
					potential += 0.5*(Q_trans * a2_real(r_n, r_n_norm)).trace();
				}
			}
		}

		//This is exactly the same as above, but uses thje fourier space sum
		for (int i = -fourier_space_iterations_; i <= fourier_space_iterations_; i++) {
			for (int j = -fourier_space_iterations_; j <= fourier_space_iterations_; j++) {
				for (int k = -fourier_space_iterations_; k <= fourier_space_iterations_; k++) {
					if(i == 0 && j == 0 && k == 0)
						continue;
					Vec h(i, j, k);
					Vec r0 = part.getPosition() - node.getCentreOfCharge();
					double h_norm = h.norm();

					potential += node.getCharge() * a0_fourier(h, r0, h_norm);
					potential += node.getDipoleMoments().dot(a1_fourier(h, r0, h_norm));
					potential += 0.5*(Q_trans * a2_fourier(h, r0, h_norm)).trace();
				}
			}
		}
		//Apply self energy correction
		potential -= 2 * alpha_ / SQRT_PI * node.getCharge();
		return potential / (3 * conf_.getDensity());
	}

	/**
	 * @brief Get force on a particle due to a node, in an infinitely periodic system.
	 *
	 * This only expands the force to dipole order. It becomes a bit of a nightmare
	 * to expand to quadrupole order, and <em>Many Body Tree Methods</em> says there
	 * is no real need.
	 *
	 * The force is given by:
	 * @f[
	 *   \vec{F} = M \vec{a}_1 + \mathbf{a}_2 \vec{D},
	 * @f]
	 * where definitions are as for getPotential().
	 *
	 * @param part		Particle to find the force upon.
	 * @param node		Node acting on the particle.
	 * @param precision	Precision to work to.
	 * @return	Force on particle due to node.
	 *
	 * TODO: Actually pay attention to precision.
	 */
	Vec getForce(const Particle<Vec>& part, const Node<Vec, Mat>& node, Precision precision) const{
		Vec force = Vec::Zero();

		//This  is pretty much the same as getPotential(), but it obeys the formula
		//given in the docs above.
		for (int i = -real_space_iterations_; i <= real_space_iterations_; i++) {
			for (int j = -real_space_iterations_; j <= real_space_iterations_; j++) {
				for (int k = -real_space_iterations_; k <= real_space_iterations_; k++) {
					Vec r_n(i*bounds_.getSize(), j*bounds_.getSize(), k*bounds_.getSize());
					r_n += (part.getPosition() - node.getCentreOfCharge());
					double r_n_norm = r_n.norm();
					if(i >= -1 && i <= 1 && j >= -1 && j <= 1 && k >= -1 && k <= 1)
						r_n_norm = sqrt(r_n_norm*r_n_norm + conf_.getForceSoftening()*conf_.getForceSoftening());

					force += node.getCharge() * a1_real(r_n, r_n_norm);
					force += a2_real(r_n, r_n_norm) * node.getDipoleMoments();
				}
			}
		}
		for (int i = -fourier_space_iterations_; i <= fourier_space_iterations_; i++) {
			for (int j = -fourier_space_iterations_; j <= fourier_space_iterations_; j++) {
				for (int k = -fourier_space_iterations_; k <= fourier_space_iterations_; k++) {
					if(i == 0 && j == 0 && k == 0)
						continue;
					Vec h(i, j, k);
					Vec r0 = part.getPosition() - node.getCentreOfCharge();
					double h_norm = h.norm();
					force += node.getCharge() * a1_fourier(h, r0, h_norm);
					force += a2_fourier(h, r0, h_norm) * node.getDipoleMoments();
				}
			}
		}

		return part.getCharge() *  force / (3 * conf_.getDensity());
	}

	//Not used
	Vec real_space_force(const Particle<Vec>& p1, const Particle<Vec>& p2){
		Vec force = Vec::Zero();
		for (int i = -real_space_iterations_; i < real_space_iterations_; i++) {
			for (int j = -real_space_iterations_; j < real_space_iterations_; j++) {
				for (int k = -real_space_iterations_; k < real_space_iterations_; k++) {
					Vec r_n(i*bounds_.getSize(), j*bounds_.getSize(), k*bounds_.getSize());
					Vec r_ni = r_n - (p2.getPosition() - p1.getPosition());
					double r = r_ni.norm();
					double r_2 = r_ni.squaredNorm();
					double r_3 = r * r_2;
					force += p1.getCharge() * p2.getCharge() *
							(r_ni / r_3) * (erfc(alpha_ * r) +
									2.0 * alpha_ * r / SQRT_PI *
									exp(-(alpha_*alpha_*r_2)));
				}
			}
		}
		return force;
	}

	//Not used
	Vec fourier_space_force(const Particle<Vec>& p1, const Particle<Vec>& p2){
		Vec force = Vec::Zero();

		double L = bounds_.getSize();

		for (int i = -fourier_space_iterations_; i < fourier_space_iterations_; i++) {
			for (int j = -fourier_space_iterations_; j < fourier_space_iterations_; j++) {
				for (int k = -fourier_space_iterations_; k < fourier_space_iterations_; k++) {
					if(i == 0 && j == 0 && k == 0)
						continue;
					Vec h(i, j, k);
					Vec k = h * 2.0 * M_PI / L;
					double r_2 = h.squaredNorm();
					force -= 2.0 / (L*L) * p1.getCharge() * p2.getCharge() * h / r_2 *
							exp(-(PI_SQUARED*r_2/(alpha_*alpha_*L*L))) *
							sin(k.dot(p2.getPosition() - p1.getPosition()));
				}
			}
		}
		return force;
	}

	//Not used
	double real_space_pot(const Particle<Vec>& p1, const Particle<Vec>& p2) {
		double potential = 0;
		for (int i = -real_space_iterations_; i < real_space_iterations_; i++) {
			for (int j = -real_space_iterations_; j < real_space_iterations_; j++) {
				for (int k = -real_space_iterations_; k < real_space_iterations_; k++) {
					Vec r_n(i*bounds_.getSize(), j*bounds_.getSize(), k*bounds_.getSize());
					Vec r_ni = r_n - (p2.getPosition() - p1.getPosition());
					potential += p2.getCharge() * erfc(alpha_ * r_ni.norm()) / r_ni.norm();
				}
			}
		}
		return potential;
	}

	//Not used
	double fourier_space_pot(const Particle<Vec>& p1, const Particle<Vec>& p2) {
		double potential = 0;

		double L = bounds_.getSize();

		for (int i = -fourier_space_iterations_; i < fourier_space_iterations_; i++) {
			for (int j = -fourier_space_iterations_; j < fourier_space_iterations_; j++) {
				for (int k = -fourier_space_iterations_; k < fourier_space_iterations_; k++) {
					if(i == 0 && j == 0 && k == 0)
						continue;
					Vec h(i, j, k);
					Vec k = h * 2.0 * M_PI / L;
					potential += p2.getCharge() / (M_PI*L) *
							(1.0 / h.squaredNorm()) *
							exp(- (M_PI*M_PI*h.squaredNorm() / (alpha_*alpha_*L*L))) *
							cos(k.dot(p2.getPosition() - p1.getPosition()));
				}
			}
		}
		return potential;
	}
private:
	/**
	 * @brief Used in other equations.
	 *
	 * Given by @f[ \frac{1}{h^2} \exp{\left[ \left( \frac{\pi h}{\alpha L} \right)^2 \right]} @f]
	 * @param h Norm of h vector.
	 * @return	Useful parameter for other equations.
	 */
	double A(double h) const{
		double L = bounds_.getSize();
		return 1.0 / (h*h) * exp(-PI_SQUARED * h*h / (alpha_*alpha_ * L*L));
	}

	double B1(double r) const{
		return erfc(alpha_ * r) + 2*alpha_*r / SQRT_PI * exp(- alpha_*alpha_*r*r);
	}

	double B2(double r) const{
		return alpha_*alpha_*alpha_ / SQRT_PI * exp(- alpha_*alpha_*r*r);
	}

	double a0_real(const Vec& r_n, double r) const{
		return erfc(alpha_ * r) / r;
	}

	double a0_fourier(const Vec& h, const Vec& r0, double h_norm) const{
		double L = bounds_.getSize();
		return 1.0 / (M_PI * L) * A(h_norm) * cos(2*M_PI / L * h.dot(r0));
	}

	Vec a1_real(const Vec& r_n, double r) const{
		return r_n * B1(r) / (r*r*r);
	}

	Vec a1_fourier(const Vec& h, const Vec& r0, double h_norm) const{
		double L = bounds_.getSize();
		return 2.0 / (L*L) * h * A(h_norm) * sin(2*M_PI / L * h.dot(r0));
	}

	Mat a2_real(const Vec& r_n, double r) const{
		double r_2 = r*r;
		double r_3 = r_2*r;
		double r_5 = r_3*r_2;
		Mat outer_product = r_n * r_n.transpose();

		Mat a2 = B1(r) * (3.0/r_5 * outer_product - Mat::Identity()/r_3);
		a2 += B2(r) * (4.0 / r_2 * outer_product);

		return a2;
	}
	Mat a2_fourier(const Vec& h, const Vec& r0, double h_norm) const{
		Mat outer_product = h * h.transpose();
		double L = bounds_.getSize();
		return (-4.0*M_PI/(L*L*L) * outer_product) * A(h_norm) * cos(2.0*M_PI/L * h.dot(r0));
	}
	const Configuration<Vec>& conf_;
	const BoundaryConditions<Vec>& bounds_;
protected:
	double alpha_;
	int real_space_iterations_, fourier_space_iterations_;
};

} /* namespace potentials */
} /* namespace treecode */
#endif /* EWALDFORCE_H_ */
