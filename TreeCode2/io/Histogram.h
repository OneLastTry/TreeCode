/*
 * Histogram.h
 *
 *  Created on: 8 Mar 2012
 *      Author: stefans
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <vector>

namespace treecode{

template <typename T>
class FilledVector {
public:
	FilledVector(T init):init_(init){}

	void set(int position, T data){
		int current_max_pos = bins_.size() - 1;
		for(int i=current_max_pos; i < position;i++)
			bins_.push_back(init_);
		bins_[position] = data;
	}

	T get(int position){
		int current_max_pos = bins_.size() - 1;
		for(int i=current_max_pos; i < position;i++)
			bins_.push_back(init_);
		return bins_[position];
	}

	std::vector<T>& getData(){
		return bins_;
	}

private:
	T init_;
	std::vector<T> bins_;
};

}

#endif /* HISTOGRAM_H_ */
