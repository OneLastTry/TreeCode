/*
 * TemperatureTracker.h
 *
 *  Created on: 8 Mar 2012
 *      Author: stefans
 */

#ifndef TEMPERATURETRACKER_H_
#define TEMPERATURETRACKER_H_

#include <vector>
#include "Histogram.h"
#include "ParticleTracker.h"

namespace treecode{
namespace output{

template <class Vec>
class TemperatureTracker : public ParticleTracker<Vec>{
public:
	TemperatureTracker(std::string filename, const std::vector<Particle<Vec>* >& parts, double ring_length, Vec origin):
		ParticleTracker<Vec>(filename, parts, ParticleTracker<Vec>::POSITION),
		ring_length_(ring_length),
		origin_(origin),
		temperature_(0), number_(0)
		{}

	virtual ~TemperatureTracker(){}

	virtual void output(){
		typedef std::vector<double> bins;
		typedef Particle<Vec> Part;

		bins b;

		BOOST_FOREACH(Part* p, this->parts_){
			double r = (p->getPosition() - this->origin).norm();
			int bin = r / ring_length_;
			temperature_.set(bin, temperature_.get(bin) + p->getVelocity()*p->getVelocity()*p->getMass()*0.5);
			number_.set(bin, number_.get(bin) + 1);
		}

		(*(this->out_)) << std::endl;
	}
private:
	double ring_length_;
	Vec origin_;
	FilledVector<double> temperature_;
	FilledVector<int> 	 number_;
};

}
}

#endif /* TEMPERATURETRACKER_H_ */
