///@file

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <Eigen/Dense>
#include "../Particle.h"
#include "../Node.h"

namespace treecode {

namespace potentials {

/**
 * @brief Precision to work with.
 * This governs the order of the force expansion.
 */
enum Precision {
	monopole,		///< Expand only to monopole order.
	dipole, 		///< Expand to dipole order.
	quadrupole		///< Expand to quadrupole order.
};

/**
 * @class Potential "potentials/Potential.h"
 */

template <class Vec, class Mat>
class Potential {
public:

	/**
	 * @brief Pure virtual method that must return force on a Particle.
	 * @param part			Particle on which the force is calculated.
	 * @param node			Node that is interacting with the particle.
	 * @param precision		Level of precision used (monopole, dipole, quadrupole)
	 * @return	Force on particle.
	 */
	virtual Vec getForce(const Particle<Vec>& part, const Node<Vec, Mat>& node, Precision precision) const = 0;
	/**
	 * @brief Pure virtual method that must return the electric potential at particle's position.
	 * @param part			Particle to calculate the potential at.
	 * @param node			Node generating potential.
	 * @param precision		Level of precision.
	 * @return	Electric potential at part's location.
	 */
	virtual double getPotential(const Particle<Vec>& part, const Node<Vec, Mat>& node, Precision precision) const = 0;
};

} /* namespace potentials */
} /* namespace treecode */
#endif /* POTENTIAL_H_ */
