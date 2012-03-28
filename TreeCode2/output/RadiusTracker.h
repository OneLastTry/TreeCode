/*
 * ParticleTracker.h
 *
 *  Created on: 3 Mar 2012
 *      Author: stefans
 */

#ifndef RADIUSTRACKER_H_
#define RADIUSTRACKER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include <map>
#include <cmath>
#include "../Particle.h"

#define CUBE(x) (x)*(x)*(x)

namespace treecode{
namespace output{

template <int D>
class RadiusTracker : public ParticleTracker<D>{
public:

	RadiusTracker(std::string filename, const std::vector<Particle<D>*>& parts, const Vec& origin, double bin_width):
		ParticleTracker<D>(filename, parts), origin_(origin), bin_width_(bin_width){}

	virtual ~RadiusTracker(){	}


	virtual void output(){
		typedef ParticleTracker<D> parent;
		typedef Particle<D> part_t;
		typedef std::map<int, double>::iterator map_it;

		bins_.erase(bins_.begin(), bins_.end());
		BOOST_FOREACH(part_t* p, parent::parts_){
			double r = (p->getPosition() - origin_).norm();
			int bin = r / bin_width_;
			bins_[bin] ++;
		}
		for(map_it it = bins_.begin(); it != bins_.end(); it++){
			double bin_start = (*it).first*bin_width_;
			double shell_volume = (4.0/3)*M_PI * (CUBE(bin_start + bin_width_) - CUBE(bin_start));
			double density = (*it).second / shell_volume;
			(*(parent::out_)) << bin_start << "\t" << density << std::endl;
		}

		(*(parent::out_)) << "============================" << std::endl;
	}
protected:
	const Vec& origin_;
	double bin_width_;
	std::map<int, double> bins_;
};

}//output namespace
}//treecode namespace

#endif /* RADIUSTRACKER_H_ */
