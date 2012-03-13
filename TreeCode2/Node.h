///@file

#ifndef NODE_H_
#define NODE_H_

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <boost/foreach.hpp>

#include "Configuration.h"
#include "Particle.h"
#include "bounds/BoundaryConditions.h"
#include "macs/AcceptanceCriterion.h"

namespace treecode {

template <class Vec, class Mat>
class Node {
public:
	/**
	 * @brief Enum representing status of node.
	 */
	enum tree_status {
		ROOT,	///< Node is root node.
		BRANCH, ///< Node is branch node.
		LEAF, 	///< Node is leaf node.
		EMPTY	///< Node contains no particles.
	};

	/**
	 * @class Node "Node.h"
	 * @brief A class representing a tree node.
	 */

	/**
	 * @brief Construct a new Node object.
	 * @param conf
	 * @param pos
	 * @param sz
	 */
	Node(const Configuration<Vec>& conf, const Vec& pos, double sz) :
			configuration(conf), position(pos), size(sz), charge(0), abs_charge(0) {}

	/**
	 * @brief Split node into separate daughters.
	 *
	 * If there are no Particle%s within the current node, then
	 * the Node's status is set to tree_status::EMPTY, and the
	 * method returns with no more action.
	 *
	 * If there is one particle within the node, the status is
	 * set to tree_status::LEAF, and the moments of the node are
	 * trivially calcualted from the one particle in the node (That is,
	 * the centre of charge, absolute charge and charge are just that of
	 *  the particle, and higher moments are zero).
	 *
	 * If there is more than one particle in the node, the status is set
	 * to tree_status::BRANCH. The method then creates @f$ 2^{d} @f$ daughter
	 * nodes, where @f$ d @f$ is the dimensions of system. The daughter nodes
	 * are arranged in quadrants/octants/@f$ 2^d @f$-thants, and the list of
	 * Particle%s in the current node is iterated over. If a Particle falls within
	 * a daughter, it is then assigned to that daughter. Finally, the split()
	 * method is called on each daughter in a recursive manner.
	 */
	void split(){
		if(getNumParticles() == 0){
			setStatus(EMPTY);
			return;
		}else if(getNumParticles() == 1){
			setStatus(LEAF);
			calculateMonopoleMoment();
			calculateDipoleMoment();
			calculateQuadrupoleMoment();
			return;
		}else{
			if(getStatus() != ROOT)
				setStatus(BRANCH);

			//Create daughters
			for (unsigned int i = 0; i < configuration.getMaxChildren(); i++) {
				Vec pos(configuration.getNumDimensions());
				//Set each component of the position
				for (unsigned int j = 0; j < configuration.getNumDimensions(); j++) {
					//Each component is either this->position[j] or
					//this->position[j] + this->size / 2.
					//This can be representing by just shifting the daughter
					//index down by the component index and ANDing with 0x01.
					pos[j] = this->position[j] + this->size / 2 * ((i >> j) & 1);
				}
				if(daughters.size() != configuration.getMaxChildren()){
					daughters.push_back(new Node<Vec,Mat>(configuration, pos, this->size/2));
				}else {
					daughters[i]->clearParticles();
					daughters[i]->setPosition(pos);
					daughters[i]->setSize(this->size/2);
				}
			}

			BOOST_FOREACH(Particle<Vec>* p, particles){
				//This finds which quadrant/octant of the parent the particle is in
				unsigned int index = 0;
				for (unsigned int j = 0; j < configuration.getNumDimensions(); j++) {
					//This will be either 0 or 1, depending on which side
					//of the parent the particle is. To flatten this into
					//A one-dimensional array, shift it up into index.
					//ie, index = 0bxyz in 3 dimensions
					int px = (p->getPosition()[j] - this->position[j]) / (this->size/2);
					index |= (px << j);
				}
				//Get the node corresponding the particle, and add the particle.
				Node* has_p = this->daughters[index];
				has_p->addParticle(p);
			}

			//Recurse into daughters.
			BOOST_FOREACH(Node* d, daughters){
				d->split();
			}
			//Finally, calculate moments from daughters' moments.
			calculateMonopoleMoment();
			calculateDipoleMoment();
			calculateQuadrupoleMoment();
		}
	}

	/**
	 * @brief Add this Node (or daughters) to interaction list for a particle.
	 *
	 * This method adds this node to the interaction list for Particle p if
	 * the multipole acceptance criteria (MAC) are fulfilled. If the MAC is
	 * not fulfilled, then this method is then called for each daughter of
	 * the current node, and so on until the MAC is either fulfilled or a
	 * leaf node is reached.
	 *
	 * @param p		Particle to generate interaction list for.
	 * @param[out] ilist	Interaction list for particle.
	 * @param bounds	Boundary conditions of system.
	 */
	void addToInteractionList(const Particle<Vec>& p, std::vector<Node*>& ilist, const BoundaryConditions<Vec>& bounds,
			const AcceptanceCriterion<Vec,Mat>& mac){

		//If the current node is empty, do nothing.
		//Similarly, if the node is a tree node and the particle is in it, do nothing.
		if(status == EMPTY || (status == LEAF && particles.front() == &p))
			return;

		//Either add to ilist or recurse into daughters.
		if(mac.accept(p, *this)){
			ilist.push_back(this);
		}else{
			BOOST_FOREACH(Node* d, daughters)
				d->addToInteractionList(p, ilist, bounds, mac);
		}
	}

	/**
	 * @brief Calculate total charge, absolute charge and centre of charge.
	 *
	 * @see Node::calculateDipoleMoment()
	 * @see Node::calculateQuadrupoleMoment()
	 */
	void calculateMonopoleMoment(){

		//If we are a leaf, just use the particle params.
		if(status == LEAF){
			charge = particles.front()->getCharge();
			abs_charge = abs(charge);
			centre_of_charge = particles.front()->getPosition();
		}else{
			//If we're a branch, use daughter moments.
			charge = 0;
			abs_charge = 0;
			centre_of_charge = Vec::Zero();
			//Iterate through daughters
			BOOST_FOREACH(Node* daughter, daughters){
				if(daughter->status != EMPTY){
					//Use daughter params
					charge += daughter->getCharge();
					abs_charge += daughter->getAbsCharge();
					centre_of_charge += daughter->getAbsCharge() * daughter->getCentreOfCharge();
				}
			}
			//Can't forget this... That would be silly.
			centre_of_charge /= getAbsCharge();
		}
	}

	/**
	 * @brief Calculate dipole moments by shifting daughter moments.
	 *
	 * See "Many Body Tree Methods" by S. Pfalzner and  P. Gibbon.
	 *
	 * @see Node::calculateMonopoleMoment()
	 * @see Node::calculateQuadrupoleMoment()
	 */
	void calculateDipoleMoment(){

		if(status == LEAF){
			//If we're a leaf, we have no dipole moment.
			dipole_moments = Vec::Zero();
		}else{
			//Otherwise, shift them from the daughters.
			dipole_moments = Vec::Zero();
			BOOST_FOREACH(Node* daughter, daughters){
				if(daughter->status != EMPTY){
					Vec disp_vec = getCentreOfCharge() - daughter->getCentreOfCharge();
					dipole_moments += daughter->getDipoleMoments() - disp_vec * daughter->getCharge();
				}
			}
		}
	}

	/**
	 * @brief Calculate quadrupole moments by shifting daughter's moments.
	 *
	 * In @f$ d @f$ dimensions, the quadrupole moment matrix of this Node,
	 * @f$ \mathbf{Q} @f$ can be
	 * calculated from the quadrupole and dipole moments of the daughters,
	 * @f$ \mathbf{Q_{di}} @f$ and @f$ \mathbf{D_{di}} @f$, by the following:
	 * @f[
	 * \mathbf{Q} = \sum_{i = 0}^{2^d} \mathbf{Q_{di}} -
	 * \vec{r}_i D_{di}^T - D_{di} \vec{r}_i^T + \vec{r}\vec{r}^T q_{di},
	 * @f]
	 * where @f$ \vec{r}_i @f$ and @f$ q_i @f$ are respectively
	 * the displacement vector and charge of daughter @f$ i @f$ (with respect
	 * to the centre of charge of the current node).
	 *
	 * @see Node::calculateMonopoleMoment()
	 * @see Node::calculateDipoleMoment()
	 */
	void calculateQuadrupoleMoment(){
		if(status == LEAF){
			quadrupole_moments = Mat::Zero();
		}else{
			quadrupole_moments = Mat::Zero();
			BOOST_FOREACH(Node* daughter, daughters){
				if(daughter->status != EMPTY){
					Vec disp_vec = getCentreOfCharge() - daughter->getCentreOfCharge();
					quadrupole_moments += (daughter->getQuadrupoleMoments()
							- disp_vec * daughter->getDipoleMoments().transpose()
							- daughter->getDipoleMoments() * disp_vec.transpose()
							+ (disp_vec*disp_vec.transpose()) * daughter->getCharge());
				}
			}
		}
	}

	/**
	 * @brief Dump nodes to file (as gnuplot rectangle objects).
	 *
	 * @param fout	Output stream.
	 * @param index	Index of current rectangle (set to zero when calling).
	 */
	void outputNodes(std::ofstream& fout, unsigned int& index){

	    fout << "set object " << ++index << " rect from "
	    		<< this->position[0] << ","
	    		<< this->position[1] << " to "
	    		<< this->position[0] + this->size << ","
	    		<< this->position[1] + this->size << std::endl;
	    if(status == LEAF)
	    	return;

	    BOOST_FOREACH(Node* d, daughters){
	    	if(d->status != EMPTY)
	    		d->outputNodes(fout, index);
	    }
	}

	/**
	 * @brief Node destructor -- destroys all daughter nodes.
	 */
	~Node() {
		BOOST_FOREACH(Node* n, daughters){
			delete n;
		}
	}

	/**
	 * @brief Add particle to node.
	 * @param p	Particle to add.
	 */
	void addParticle(Particle<Vec>* p){particles.push_back(p);};

	/**
	 * @brief Get number of particles contained within this Node.
	 * @return	Number Particle%s in node.
	 */
	unsigned int getNumParticles() const {return particles.size(); }

	/**
	 * @brief Set particles in bulk form.
	 * @param p	Vector of Particle%s to assign.
	 */
	void setParticles(std::vector<Particle<Vec>*> p){particles = p;}
	/**
	 * @brief Set position (origin) of node.
	 *
	 * This method accepts a const reference, and copies it into
	 * the position vector of the current Node object.
	 * @param pos	New position.
	 */
	void setPosition(const Vec& pos){position = pos;}
	/**
	 * @brief Set size of Node.
	 * @param sz	Size of node.
	 */
	void setSize(double sz){size = sz;}

	/**
	 * @brief Get rank of current Node in the tree.
	 * @return	Status/rank of current Node.
	 */
	Node::tree_status getStatus() const{return status;}
	/**
	 * @brief Set status/rank in tree.
	 * @param s	Status.
	 */
	void setStatus(Node::tree_status s){this->status = s;}

	/**
	 * @brief Get total charge of all Particle%s in current Node.
	 * @return	Total charge contained in current Node.
	 */
	double getCharge() const {return charge;}
	/**
	 * @brief Get sum of modulus of all Particle%s in current Node.
	 * @return	Absolute charge of node.
	 */
	double getAbsCharge() const {return abs_charge;}

	/**
	 * @brief Get centre of charge of all particles in node.
	 * @return	Centre of charge.
	 */
	const Vec& getCentreOfCharge() const {return centre_of_charge;}
	/**
	 * @brief Get dipole moments abount centre of charge, from all Particle%s in node.
	 * @return	Dipole moments.
	 */
	const Vec& getDipoleMoments() const {return dipole_moments; }
	/**
	 * @brief Get quadrupole moments about centre of charge, from all Particle%s in node.
	 * @return	Quadrupole moments.
	 */
	const Mat& getQuadrupoleMoments() const {return quadrupole_moments; }

	/**
	 * @brief Delete all particles from node.
	 */
	void clearParticles() { particles.clear(); }

	/**
	 * @brief Get list of particles.
	 * @return Particles.
	 */
	const std::vector<Particle<Vec>* >& getParticles() const {return particles;}

	/**
	 * @brief Get size of node.
	 * @returns Size of node.
	 */
	double getSize() const {return size; }

private:
	const Configuration<Vec>& configuration;
	Vec position;
	double size;
	std::vector<Particle<Vec>*> particles;
	std::vector<Node<Vec,Mat>*> daughters;
	Node::tree_status status;

	double charge, abs_charge;

	Vec centre_of_charge;
	Vec dipole_moments;
	Mat quadrupole_moments;
};

} /* namespace treecode */
#endif /* NODE_H_ */
