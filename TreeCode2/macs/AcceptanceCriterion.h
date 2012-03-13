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

template <class V, class M> class Node;

template <class Vec, class Mat>
class AcceptanceCriterion{
public:
	virtual bool accept(const Particle<Vec>& p, const Node<Vec,Mat>& n) const = 0;
};
}

#endif /* ACCEPTANCECRITERION_H_ */
