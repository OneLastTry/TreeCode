/*
 * AcceptanceCriterion.h
 *
 *  Created on: 10 Mar 2012
 *      Author: stefans
 */

#ifndef ACCEPTANCECRITERION_H_
#define ACCEPTANCECRITERION_H_

#include "../Particle.h"

namespace treecode{

template <int D> class Node;

template <int D>
class AcceptanceCriterion{
public:
	enum result {
		ACCEPT,
		CONTINUE,
		REJECT
	};

	virtual result accept(const Particle<D>& p, const Node<D>& n) const = 0;
};
}

#endif /* ACCEPTANCECRITERION_H_ */
