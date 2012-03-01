/*
 * Configuration.h
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

namespace treecode {

template<class Vec>
class Configuration {
public:

	/**
	 * @class Configuration "Configuration.h"
	 * @brief Global configuration class.
	 */

	/**
	 * @brief Construct new Configuration object.
	 * @param num_dims	Number of dimensions we are working in.
	 * @param mac		Multipole acceptance criterion (@f$ \theta @f$).
	 * @param dt		Timestep.
	 * @param mt		Maximum time to simulate until.
	 * @param dens		Number of electrons in a Debye sphere/circle.
	 * @param fs		Force softening param.
	 */
	Configuration(unsigned int num_dims, double mac, double dt, double mt, double dens, double fs) :
			num_dimensions(num_dims), theta(mac), delta_t(dt), max_time(mt), density(dens), force_softening(fs){}

	///Destructor. Does nothing.
	~Configuration() {}


	/**
	 * @brief Get number of dimensions we are working in.
	 * @return	Number of dimensions.
	 */
    unsigned int getNumDimensions() const{return num_dimensions;}
    /**
     * @brief Set number of dimensions.
     * @param num_dimensions	Number of dimensions.
     */
    void setNumDimensions(unsigned int num_dimensions){this->num_dimensions = num_dimensions;}

    /**
     * @brief Get @f$ \theta @f$.
     * @return	@f$ \theta @f$
     */
    double getTheta() const {return theta;}
    /**
     * @brief Set @f$ \theta @f$.
     * @param t	@f$ \theta @f$
     */
    void setTheta(double t){this->theta = t;}

    /**
     * @brief Get maximum number of children in these dimensions (@f$ 2^d @f$).
     * @return	Maximum number of children.
     */
    unsigned int getMaxChildren() const {return (1u << num_dimensions); }

    /**
     * @brief Get timestep.
     * @return	Timestep.
     */
    double getTimestep() const {return delta_t; }
    /**
     * @brief Set timestep.
     * @param dt	Timestep.
     */
    void setTimestep(double dt) { delta_t = dt; }

    /**
     * @brief Get time to simulate until.
     * @return Maximum time.
     */
    double getMaxTime() const {return max_time;}

    /**
     * @brief Set maximum time to simulate until.
     * @param mt Maximum time.
     */
    void setMaxTime(double mt) { max_time = mt; }

    /**
     * @brief Get the number of <em>electrons</em> in a Debye sphere/circle.
     * @return Number of electrons in a Debye sphere.
     */
    double getDensity() const {return density;}

    /**
     * @brief Set the number of <em>electrons</em> in a Debye sphere/circle.
     * @param d Number of electrons in a Debye sphere.
     */
    void setDensity(double d) { density = d; }

    /**
     * @brief Get force softening parameter.
     * @return Force softening parameter.
     */
    double getForceSoftening() const {return force_softening;}

    /**
     * @brief Set force softening parameter.
     * @param fs Force softening parameter.
     */
    void setForceSoftening(double fs){force_softening = fs;}

private:
	unsigned int num_dimensions; 	///< Number of dimensions
	double theta;					///< MAC
	double delta_t;					///< Timestep
	double max_time;				///< Maximum time to simulate until.
	double density;					///< Number of electrons in a Debye sphere/circle.
	double force_softening;			///< Force softening parameter;
};

} /* namespace treecode */
#endif /* CONFIGURATION_H_ */
