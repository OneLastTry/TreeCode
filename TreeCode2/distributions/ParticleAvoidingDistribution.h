///@file

#ifndef UNIFORMDISTRIBUTION_H_
#define UNIFORMDISTRIBUTION_H_

#include <boost/random/uniform_01.hpp>

#include "Distribution.h"

#include <vector>
#include <boost/foreach.hpp>
#include "../drude/BigParticle.h"

namespace treecode {
namespace distribution {

template <class RNG, int D, class T>
class ParticleAvoidingDistribution : public VectorDistribution<RNG, D>{
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	ParticleAvoidingDistribution(const Vec& min, const Vec& max, const std::vector<T*>& big_particles) :
		minimum(min), maximum(max), big_particles_(big_particles){}

	virtual Vec getVector(RNG& rng) const  {
		boost::uniform_01<double> dist;
		bool accept_particle = true;
		Vec v;

		do{
			for (unsigned int j = 0; j < v.rows(); j++)
				v[j] = dist(rng) * (maximum[j] - minimum[j]) + minimum[j];

			//Reject vector if it intersects with any particles.
			BOOST_FOREACH(T* bp, big_particles_){
				if( (v - bp->getPosition()).norm() < bp->getRadius())
					continue;
			}
		}while(!accept_particle);

		return v;
	}
	/**
	 * Destructor. Does nothing.
	 */
	virtual ~ParticleAvoidingDistribution(){}
private:
	Vec minimum, maximum;
	const std::vector<T*> big_particles_;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* UNIFORMDISTRIBUTION_H_ */
