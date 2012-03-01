/*
 * CountingPeriodicBounds.h
 *
 *  Created on: 20 Feb 2012
 *      Author: stefans
 */

#ifndef COUNTINGPERIODICBOUNDS_H_
#define COUNTINGPERIODICBOUNDS_H_

#include "../Configuration.h"
#include "PeriodicBoundary.h"
#include <fstream>
#include <iostream>

namespace treecode {

template <class Vec>
class CountingPeriodicBounds : public PeriodicBoundary<Vec> {
public:
	CountingPeriodicBounds(Configuration<Vec> c, Vec origin, double length, std::ofstream& output) :
		PeriodicBoundary<Vec>(c, origin, length),
		output_(output),
		left_edge_(Vec::Zero()), right_edge_(Vec::Zero()){}

	/**
	 * @brief Shift particle at edges of box.
	 * When a particle is moved, this method should be called. It translates
	 * the particle to the opposite edge of the simulation region when it
	 * moves out.
	 *
	 * @param p	Particle to move.
	 */
	void particleMoved(treecode::Particle<Vec>* p){
		Vec translation_vector = Vec::Zero();

		for (int i = 0; i < p->getPosition().rows(); i++) {
			double component = p->getPosition()[i];
			if(component < PeriodicBoundary<Vec>::origin_[i]){
				if(p->getCharge() == -1)
					left_edge_[i] += 1;
				translation_vector[i] = PeriodicBoundary<Vec>::length_;
			}
			else if(component > PeriodicBoundary<Vec>::origin_[i] + PeriodicBoundary<Vec>::length_){
				if(p->getCharge() == -1)
					right_edge_[i] += 1;
				translation_vector[i] = -PeriodicBoundary<Vec>::length_;
			}
		}
		p->updatePosition(translation_vector);
	}

	/**
	 * @brief Output the number of particles that have passed each edge.
	 */
	void timestepOver(){
		for(int i = 0; i < left_edge_.rows(); i++){
			output_ << left_edge_[i] << "\t" << right_edge_[i] << "\t";
		}
		output_ << std::endl;
		output_.flush();

		left_edge_ = Vec::Zero();
		right_edge_ = Vec::Zero();
	}

	virtual ~CountingPeriodicBounds(){}
private:
	std::ofstream& output_;
	Vec left_edge_, right_edge_;
};

} /* namespace treecode */
#endif /* COUNTINGPERIODICBOUNDS_H_ */
