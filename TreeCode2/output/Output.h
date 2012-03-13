/*
 * Output.h
 *
 *  Created on: 3 Feb 2012
 *      Author: stefans
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <string>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include "../Particle.h"

namespace treecode {

namespace output {

template <class Vec>
class Output {
public:
	enum Record {
		COORDINATES,
		VELOCITIES,
		TIME,
		KINETIC_ENERGY,
		POTENTIAL_ENERGY
	};

	Output(std::string filename, bool append = false, std::string separator = "\t"):separator_(separator){
		std::ios_base::openmode open_mode = std::ios_base::out;
		if(append)
			open_mode |= std::ios_base::app;
		out_ = new std::ofstream(filename, open_mode);
	}

	void particleOutput(const Particle<Vec>& p) const{
		BOOST_FOREACH(Output::Record rec, records_){
			switch(rec){
			case Record::COORDINATES:
				for(int i=0;i<p.getPosition().rows();i++)
					(*out_) << p.getPosition()[i] << separator_;
				break;
			case Record::VELOCITIES:
				for(int i=0;i<p.getVelocity().rows();i++)
					*out_ << p.getVelocity()[i] << separator_;
				break;
			default:
				break;
			}
		}
	}

	void timestepOutput(double time, double ke, double pe) const{
		BOOST_FOREACH(Record rec, records_){
			switch(rec){
			case Record::TIME:
				*out_ << time << separator_;
				break;
			case Record::KINETIC_ENERGY:
				*out_ << ke << separator_;
				break;
			case Record::POTENTIAL_ENERGY:
				*out_ << pe << separator_;
				break;
			default:
				break;
			}
		}
	}

	void recordComplete() const{
		*out_ << std::endl;
	}

	void addRecord(Record rec){
		records_.push_back(rec);
	}

	~Output(){
		out_->close();
		delete out_;
	}

private:
	std::vector<Record> records_;
	std::ofstream* out_;
	std::string separator_;
};

} /* namespace output */
} /* namespace treecode */
#endif /* OUTPUT_H_ */
