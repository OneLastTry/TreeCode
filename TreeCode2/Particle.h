#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <Eigen/Dense>
#include <vector>

#include <iostream>

#include "distributions/ChargeDistribution.h"
#include "distributions/Distribution.h"

namespace treecode {

template <int D> class Node;

template <int D>
class Particle{
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	/**
	 * @class Particle "Particle.h"
	 * @brief Class representing a particle.
	 *
	 * This class contains only the necessary information:
	 * Charge, mass, position and velocity. There is also
	 * a configuration object that must be set to the
	 * global configuration, but this is stored as a const
	 * reference, and so should take negligible resources.
	 */

	/**
	 * @brief Construct a new particle.
	 * @param q		Charge.
	 * @param m		Mass.
	 * @param pos	Position.
	 * @param vel	Velocity.
	 * @param id	A unique global ID (for data storage).
	 */
	Particle(double q, double m, const Vec& pos, const Vec& vel) :
			charge(q), mass(m), position(pos), velocity(vel){}

	/**
	 * @brief Destructor. Nothing needs to be done.
	 */
	virtual ~Particle() {}


	/**
	 * @brief Get charge on particle.
	 * @return	Particle charge.
	 */
    double getCharge() const{return charge;}
    /**
     * @brief Set charge on particle.
     * @param charge Charge to set to.
     */
    void setCharge(double charge){this->charge = charge;}

    /**
     * @brief Get mass of particle.
     * @return	Mass of particle.
     */
    double getMass() const{return mass;}
    /**
     * @brief Set mass of particle.
     * @param mass	Mass to set to.
     */
    void setMass(double mass){this->mass = mass;}

    /**
     * @brief Get position of particle.
     * @return	Position of particle.
     */
    Vec getPosition() const{return position;}	//TODO: Should this return a reference?
    /**
     * @brief Get velocity of particle.
     * @return	Velocity of particle.
     */
    Vec getVelocity() const{return velocity;}

    /**
     * @brief Set position.
     * @param pos	Position to move to.
     */
    void setPosition(const Vec& pos){ position = pos; }

    /**
     * @brief Set velocity.
     * @param vel New velocity.
     */
    void setVelocity(const Vec& vel){ velocity = vel; }

    /**
     * @brief Increment velocity of particle.
     * @param dv Amount to increase by.
     */
    void updateVelocity(Vec dv){ velocity += dv; }
    /**
     * @brief Increment position of particle.
     * @param dx Amount to increase by.
     */
    void updatePosition(Vec dx){ position += dx; }

    /**
     * @brief Generate particles corresponding to supplied distributions.
     *
     * @tparam RNG				A boost random number generator, such as mt19937.
     * @param c					Global configuration.
     * @param num_particles		Number of particles to generate.
     * @param mass				Mass of particles to generate.
     * @param rng				Boost random number generator.
     * @param position_dist		Position distribution.
     * @param velocity_dist		Velocity distribution.
     * @param charge_dist		Charge distribution.
     * @param[in,out] id		Sequential ID number. Incremented with each particle.
     * @return					Vector of Particle%s obeying supplied parameters.
     */
    template <class RNG>
    static std::vector<Particle*>  generateParticles(unsigned int num_particles,
    		double mass, RNG& rng,
    		const distribution::VectorDistribution& position_dist,
    		const distribution::VectorDistribution& velocity_dist,
    		const distribution::ChargeDistribution& charge_dist){
    	//Create the vector and reserve the number of particles.
    	std::vector<Particle<D>*> parts;
    	parts.reserve(num_particles);
    	for (unsigned int i = 0; i < num_particles; i++) {
    		//Get charge, position and velocity from supplied distributions.
    		int charge = charge_dist.getCharge();
    		Vec pos = position_dist.getVector();
    		Vec vel = velocity_dist.getVector();
    		//Generate a new particle, and add to vector.
    		Particle* p = new Particle<D>(charge, mass, pos, vel);
    		parts.push_back(p);
    	}
    	return parts;
    }

    /**
     * @brief Delete (freeing memory) a Vector of Particle pointers.
     * @param parts	Vector of Particle pointers to free.
     */
    static void deleteParticles(std::vector<Particle*>& parts){
    	//Just go through all particles and baleet them all.
    	while(parts.size() > 0){
    		Particle* p = parts.back();
    		delete p;
    		parts.pop_back();
    	}
    }

private:
	double charge, mass;
	Vec position, velocity;
};

} /* namespace treecode */
#endif /* PARTICLE_H_ */
