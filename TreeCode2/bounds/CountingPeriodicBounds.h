/*
 * CountingPeriodicBounds.h
 *
 *  Created on: 20 Feb 2012
 *      Author: stefans
 */

#ifndef COUNTINGPERIODICBOUNDS_H_
#define COUNTINGPERIODICBOUNDS_H_

#include "PeriodicBoundary.h"
#include <fstream>
#include <iostream>

namespace treecode {

template <int D>
class CountingPeriodicBounds : public PeriodicBoundary<D> {
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	CountingPeriodicBounds(Vec origin, double length, std::ostream& output) :
		PeriodicBoundary<D>(origin, length),
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
	void particleMoved(treecode::Particle<D>* p){
		Vec translation_vector = Vec::Zero();

		for (int i = 0; i < p->getPosition().rows(); i++) {
			double component = p->getPosition()[i];
			if(component < PeriodicBoundary<D>::origin_[i]){
				if(p->getCharge() == -1)
					left_edge_[i] += 1;
				translation_vector[i] = PeriodicBoundary<D>::length_;
			}
			else if(component > PeriodicBoundary<D>::origin_[i] + PeriodicBoundary<D>::length_){
				if(p->getCharge() == -1)
					right_edge_[i] += 1;
				translation_vector[i] = -PeriodicBoundary<D>::length_;
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
	std::ostream& output_;
	Vec left_edge_, right_edge_;
};

} /* namespace treecode */
#endif /* COUNTINGPERIODICBOUNDS_H_ */
