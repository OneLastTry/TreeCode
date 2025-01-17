///@file

#ifndef PUSHER_H_
#define PUSHER_H_

#include <vector>
#include "../potentials/Potential.h"
#include "../Particle.h"
#include "../Node.h"
#include "../Tree.h"
#include "../bounds/BoundaryConditions.h"
#include "../macs/AcceptanceCriterion.h"

namespace treecode {

namespace pusher{

template <int D>
class Pusher{
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

	/**
	 * @class Pusher
	 * @brief Abstract base class for particle pushers.
	 */


public:
	/**
	 * @brief Pure virtual base function for pushing a particle.
	 * @param parts	List of particles to push.
	 * @param tree	Tree to use to push particles.
	 * @param bc	Boundary conditions, to be updated as particles move.
	 * @param prec	Precision to use in force calculation.
	 * @return std::pair<ke,pe> containg kinetic energy and potential energy
	 */
	virtual std::pair<double, double> push_particles(
			std::vector<Particle<D>*> parts,
			Tree<D>& tree,
			BoundaryConditions<D>& bc,
			potentials::Precision prec,
			const AcceptanceCriterion<D>& mac) = 0;

};

}
}

#endif /* PUSHER_H_ */
