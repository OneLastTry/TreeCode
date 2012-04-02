/*
 * Tree.h
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */

#ifndef TREE_H_
#define TREE_H_

#include <vector>
#include <fstream>
#include <Eigen/Dense>

#include "bounds/BoundaryConditions.h"
#include "Particle.h"
#include "Node.h"
#include "macs/AcceptanceCriterion.h"

namespace treecode {

template <int D>
class Tree {
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	/**
	 * @class Tree "Tree.h"
	 * @brief A representation of a tree.
	 *
	 * This class represents a tree. That is, it has both a root node and
	 * some boundary conditions. These are represented by a Node and
	 * BoundaryConditions class, respectively.
	 */

	/**
	 * @brief Construct a Tree.
	 * @param bc		Boundary conditions the tree will obey.
	 * @param parts		Vector of particles that the tree is constructed around.
	 */
	Tree(const BoundaryConditions<D>& bc, const std::vector<Particle<D>*>& parts):
		boundary(bc), particles(parts){
		//Create a root node using config and boundary conditions
		root = new Node<D>(bc.getOrigin(), bc.getSize());
		//Populate with particles, and gently inform the node it is root
		root->setParticles(parts);
		root->setStatus(Node<D>::ROOT);
	}

	/**
	 * @brief Causes a new quad/oct-tree to be generated.
	 */
	void rebuild(){
		delete root;
		//Test memory usage if deleting root at each timestep
		root = new Node<D>(boundary.getOrigin(), boundary.getSize());
		//Populate with particles, and gently inform the node it is root
		root->setParticles(particles);
		root->setStatus(Node<D>::ROOT);
		root->split();
	}

	/**
	 * @brief Dump nodes to a file.
	 * @param fout	Output stream.
	 */
	void dumpNodes(std::ofstream& fout){
		unsigned int index = 0;
		root->outputNodes(fout, index);
	}

	/**
	 * @brief Generate interaction list for a particle.
	 * @param p		Particle to generate the list for.231
	 * @param[out] ilist	Nodes to interact with will be placed here.
	 */
	void getInteractionList(const Particle<D>& p, std::vector<Node<D>*>& ilist, const AcceptanceCriterion<D>& mac) const{
		//Make sure the list is clear, and then delegate to the root node.
		ilist.clear();
		root->addToInteractionList(p, ilist, boundary, mac);
	}

	/**
	 * @brief Get the centre of charge of the entire system.
	 * @return	Centre of charge.
	 */
	Vec getTotalCentreOfCharge() const{
		return root->getCentreOfCharge();
	}

	/**
	 * @brief Get the dipole moments of the entire system.
	 * @return	The dipole moment.
	 */
	Vec getTotalDipoleMoments() const{
		return root->getDipoleMoments();
	}

	/**
	 * @brief Get the quadrupole moments of the entire system.
	 *
	 * @note Note that all quadrupole moments used in this program
	 * are <em>intrinsic</em> quadrupole moments, <em>not</em>
	 * traceless quadrupole moments.
	 * @return	The quadrupole moments.
	 */
	Mat getTotalQuadrupoleMoments() const{
		return root->getQuadrupoleMoments();
	}

	/**
	 * @brief Get total charge of the system.
	 * @return	Total charge of the entire system (usually approximately zero).
	 */
	int getTotalCharge() const {
		return root->getCharge();
	}

	/**
	 * @brief Get the sum of the modulus of the charge on every particle.
	 * @return	Total absolute charge of the system.
	 */
	int getTotalAbsCharge() const {
		return root->getAbsCharge();
	}

	/**
	 * @brief Destructor for Tree class.
	 *
	 * This deletes the root node, which cascades through every
	 * child node.
	 */
	~Tree() {
		delete root;
	}

	/**
	 * Get root node.
	 * @return	Root node.
	 */
	const Node<D>& getRoot() const{return *root; }

private:
	const BoundaryConditions<D>& boundary;
	const std::vector<Particle<D>*>& particles;
	Node<D>* root;
};

} /* namespace treecode */
#endif /* TREE_H_ */
