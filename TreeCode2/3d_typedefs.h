/*
 * 3d_typedefs.h
 *
 *  Created on: 10 Feb 2012
 *      Author: stefans
 */

#ifndef THREE_D_TYPEDEFS_H_
#define THREE_D_TYPEDEFS_H_
#include <Eigen/Dense>

#include "bounds/OpenBoundary.h"
#include "bounds/PeriodicBoundary.h"
#include "bounds/BoundaryConditions.h"

#include "distributions/UniformDistribution.h"
#include "distributions/SphericalDistribution.h"
#include "distributions/ConstantChargeDistribution.h"
#include "distributions/ConstDistribution.h"
#include "distributions/MaxwellDistribution.h"
#include "distributions/SinusoidalDistribution.h"

#include "potentials/CoulombForce.h"
#include "potentials/Potential.h"
#include "potentials/InterpolatedEwaldSum.h"

#include "pushers/LeapfrogPusher.h"

#include "TimeIntegrator.h"

typedef Eigen::Vector3d Vec;
typedef Eigen::Matrix3d Mat;

namespace treecode {
typedef Particle<Vec> 				Particle3d;
typedef Node<Vec,Mat>				Node3d;
typedef Tree<Vec,Mat>				Tree3d;
typedef TimeIntegrator<Vec, Mat> 	TimeIntegrator3d;
typedef BoundaryConditions<Vec>		BoundaryConditions3d;
typedef OpenBoundary<Vec>			OpenBoundary3d;
typedef PeriodicBoundary<Vec>		PeriodicBoundary3d;
typedef Configuration<Vec> 			Configuration3d;

namespace pusher{
typedef LeapfrogPusher<Vec, Mat> 							LeapfrogPusher3d;
typedef Pusher<Vec, Mat>									Pusher3d;
}

namespace potentials{
typedef CoulombForceThreeD<Vec, Mat> 						CoulombForce3d;
typedef EwaldForce<Vec, Mat>								EwaldForce3d;
typedef InterpolatedEwaldSum<Vec,Mat>						InterpolatedEwaldSum3d;
}

namespace distribution{
typedef UniformDistribution<boost::mt19937, Vec> 	UniformDistribution3d;
typedef SphericalDistribution<boost::mt19937, Vec>	SphericalDistribution3d;
typedef ConstDistribution<boost::mt19937, Vec>		ConstDistribution3d;
typedef MaxwellDistribution<boost::mt19937, Vec>	MaxwellDistribution3d;
typedef SinusoidalDistribution<boost::mt19937, Vec>	SinusoidalDistribution3d;
typedef ConstantChargeDistribution<boost::mt19937>	ConstantChargeDistribution3d;
}

}

#endif /* THREE_D_TYPEDEFS_H_ */
