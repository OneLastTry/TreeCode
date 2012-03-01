/*
 * OpenBoundary.h
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */

#ifndef OPENBOUNDARY_H_
#define OPENBOUNDARY_H_

#include "BoundaryConditions.h"

#include <vector>
#include <Eigen/Dense>
#include "../Particle.h"
#include "../Configuration.h"

namespace treecode {

template <class Vec>
class OpenBoundary : public BoundaryConditions<Vec> {
public:
	/**
	 * @class OpenBoundary
	 * @brief Class representing open boundaries.
	 *
	 * The most relevant feature of this class is that the
	 * size and position change dynamically with the particles.
	 * That means it can reliably be used to model an open-boundary
	 * system.
	 */

	/**
	 * @brief Construct a new boundary system.
	 * @param conf	Configuration.
	 */
	OpenBoundary(const Configuration<Vec>& conf) :
			configuration(conf){}

	/**
	 * @brief Initialise with list of particles.
	 * This <em>must</em> be called before starting tree builds.
	 * Otherwise, the tree will not know how big the system is initially.
	 *
	 * @param parts Particles in system.
	 */
	void init(const std::vector<treecode::Particle<Vec>*>& parts){
		reset(parts.front());
		for(Particle<Vec>* p : parts){
			particleMoved(p);
		}
	}

	/**
	 * @brief Gets displacement vector, as one would expect.
	 * @param r1	Position vector 1.
	 * @param r2	Position vector 2.
	 * @return		@f$ \vec{r}_1 - \vec{r}_2 @f$
	 */
	Vec getDisplacementVector(const Vec& r1, const Vec& r2) const{
		return r1 - r2;
	}

	/**
	 * @brief Resets system to sane initial values.
	 * @param p	An abitrary particle in the system.
	 */
	void reset(Particle<Vec>* p){
		//Just reinit so that particleMoved() works properly.
		minimum = p->getPosition().array();
		maximum = p->getPosition().array();
		flag = 0;
	}

	/**
	 * @brief Resize bounds when particle moves.
	 * @param p	Particle moved.
	 */
	void particleMoved(treecode::Particle<Vec>* p){
		//If we are at the start of a new timestep, reset.
		if(flag != 0)
			reset(p);
		minimum = minimum.array().min(p->getPosition().array());
		maximum = maximum.array().max(p->getPosition().array());
	}

	/**
	 * @brief Signal that the timestep is over.
	 */
	void timestepOver(){
		//Just signal that we need to be reset next particle movement.
		flag = 1;
	}

	/**
	 * @brief Get origin.
	 * @return Origin.
	 */
	Vec getOrigin() const{
		return (minimum.array() - 0.01).matrix();
	}

	/**
	 * @brief Get length of a side of the system.
	 * @return System size.
	 */
	double getSize() const{
		return (maximum - minimum).maxCoeff() + 0.02;
	}

	/**
	 * @brief Destructor.
	 */
	~OpenBoundary() {}


private:

	const Configuration<Vec>& configuration;
	Vec minimum, maximum;
	short flag;	//Needs updating if non-zero
};

} /* namespace treecode */
#endif /* OPENBOUNDARY_H_ */
