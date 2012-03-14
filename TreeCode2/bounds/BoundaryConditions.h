///@file

#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_

#include <Eigen/Dense>
#include "../Particle.h"

namespace treecode {

template <class Vec, class Mat>
class BoundaryConditions {
	/**
	 * @class BoundaryConditions
	 * @brief Abstract base class representing boundary conditions of the system.
	 */

public:
	/**
	 * @brief Called whenever a particle is moved.
	 * This can be useful for various things, such as resizing the simulation
	 * region or reflecting particles at boundaries.
	 *
	 * @param p Particle that has moved.
	 */
	virtual void particleMoved(treecode::Particle<Vec,Mat>* p) = 0;

	/**
	 * @brief Called whenever a timestep is over.
	 */
	virtual void timestepOver() = 0;

	/**
	 * @brief Get this boundary conditions' definition of a displacement vector.
	 * @param r1	Position vector 1.
	 * @param r2	Position vector 2.
	 * @return		Displacement vector, adjusted for boundary conditions.
	 */
	virtual Vec getDisplacementVector(const Vec& r1, const Vec& r2) const = 0;

	/**
	 * @brief Get origin of system.
	 * @return	Origin of system.
	 */
	virtual Vec getOrigin() const = 0;

	/**
	 * @brief Get size of each side of the system.
	 * @return Length of each side.
	 */
	virtual double getSize() const = 0;
};

} /* namespace treecode */
#endif /* BOUNDARYCONDITIONS_H_ */
