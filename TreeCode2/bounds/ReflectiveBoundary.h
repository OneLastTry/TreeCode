/*
 * PeriodicBoundary.h
 *
 *  Created on: 20 Jan 2012
 *      Author: stefans
 */

#ifndef PERIODICBOUNDARY_H_
#define PERIODICBOUNDARY_H_

#include "BoundaryConditions.h"

namespace treecode {

template <int D>
class ReflectiveBoundary : public BoundaryConditions<D> {
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:

	/**
	 * @brief Create a new set of periodic boundary conditions.
	 * @param origin	Origin of system.
	 * @param length	Length of each side of the system.
	 */
	ReflectiveBoundary(const Vec origin, double length):
		origin_(origin), length_(length){}

	/**
	 * @brief Destructor (does nothing).
	 */
	virtual ~PeriodicBoundary() {
	}

	virtual void particleMoved(treecode::Particle<D>* p){
		Vec pos_vec = p->getPosition();
		Vec vel_vec = p->getVelocity();
		for(int i = 0; i < pos_vec.rows(); i++){
			if(pos_vec[i] < origin_[i]){
				pos_vec[i] += 2*(pos_vec[i] - origin_[i]);
				vel_vec[i] *= -1;
			}else if(pos_vec[i] > origin_[i] + length_){
				pos_vec[i] -= 2*(pos_vec[i] - (origin_[i] + length));
				vel_vec[i] *= -1;
			}
		}
		p->setPosition(pos_vec);
		p->setVelocity(vel_vec);
	}

	virtual void timestepOver(){}

	virtual Vec getDisplacementVector(const Vec& r1, const Vec& r2) const{
		Vec disp_vec = r1 - r2;
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
	Vec origin_;
	double length_;
};

} /* namespace treecode */
#endif /* PERIODICBOUNDARY_H_ */
