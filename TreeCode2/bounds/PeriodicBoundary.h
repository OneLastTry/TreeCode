/*
 * PeriodicBoundary.h
 *
 *  Created on: 20 Jan 2012
 *      Author: stefans
 */

#ifndef PERIODICBOUNDARY_H_
#define PERIODICBOUNDARY_H_

#include "BoundaryConditions.h"
#include "../Configuration.h"

namespace treecode {

template <class Vec, class Mat>
class PeriodicBoundary : public BoundaryConditions<Vec,Mat> {
public:
	/**
	 * @class PeriodicBoundary "bounds/PeriodicBoundary.h"
	 * @brief Class representing a periodic boundary condition.
	 * Whenever a particle moves out of the system, it is translated back into the system.
	 *
	 */

	/**
	 * @brief Create a new set of periodic boundary conditions.
	 * @param conf		Configuration.
	 * @param origin	Origin of system.
	 * @param length	Length of each side of the system.
	 */
	PeriodicBoundary(const Configuration<Vec>& conf, const Vec origin, double length):
		conf_(conf), origin_(origin), length_(length){}

	/**
	 * @brief Destructor (does nothing).
	 */
	virtual ~PeriodicBoundary() {
	}

	/**
	 * @brief Shift particle at edges of box.
	 * When a particle is moved, this method should be called. It translates
	 * the particle to the opposite edge of the simulation region when it
	 * moves out.
	 *
	 * @param p	Particle to move.
	 */
	virtual void particleMoved(treecode::Particle<Vec,Mat>* p){
		Vec translation_vector = Vec::Zero();

		for (int i = 0; i < p->getPosition().rows(); i++) {
			double component = p->getPosition()[i];
			if(component < origin_[i])
				translation_vector[i] = length_;
			else if(component > origin_[i] + length_)
				translation_vector[i] = -length_;
		}
		p->updatePosition(translation_vector);
	}

	/**
	 * @brief In the case of a periodic boundary, this does not do anything.
	 */
	virtual void timestepOver(){}

	/**
	 * @brief Find the shortest displacement vector between two points, taking periodic images into account.
	 *
	 * Simply, this checks each component of the displacement vector. If that component is greater
	 * than half the ``length'' of the system, then it is reduced by the length of the system.
	 *
	 * @param r1	Position vector to find the displacement vector from.
	 * @param r2	Position vector at end of displacement vector.
	 * @return		Shortest displacement vector from r1 to r2.
	 */
	virtual Vec getDisplacementVector(const Vec& r1, const Vec& r2) const{
		Vec disp_vec = r1 - r2;
		for (int i = 0; i < disp_vec.rows(); i++) {
			if(disp_vec[i] > length_ / 2)
				disp_vec[i] -= length_;
			else if(disp_vec[i] < -length_ / 2)
				disp_vec[i] += length_;
		}
		return disp_vec;
	}

	/**
	 * @brief Get the origin of the boundary conditions.
	 * @return	Origin of the boundary.
	 */
	virtual Vec getOrigin() const{
		return origin_;
	}

	/**
	 * @brief Get the size of the system.
	 * @return	Length of each side of the system.
	 */
	virtual double getSize() const {
		return length_;
	}

protected:
	const Configuration<Vec>& conf_;
	Vec origin_;
	double length_;
};

} /* namespace treecode */
#endif /* PERIODICBOUNDARY_H_ */
