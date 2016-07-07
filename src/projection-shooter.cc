//crabe8pince
// Copyright (c) 2015 CNRS
// Authors: Mylene Campana
//
// This file is part of hpp-core
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

# include <hpp/rbprm/projection-shooter.hh>
# include <sstream>
# include <hpp/util/debug.hh>
# include <hpp/fcl/distance.h>
# include <hpp/model/collision-object.hh>
# include <hpp/model/configuration.hh>
# include <hpp/model/device.hh>
# include <hpp/model/joint.hh>
# include <hpp/model/joint-configuration.hh>
# include <hpp/core/configuration-shooter.hh>
# include <hpp/core/config-validations.hh>
# include <hpp/core/distance-between-objects.hh>
# include <hpp/rbprm/fullbodyBallistic/parabola-library.hh>
# include <hpp/rbprm/rbprm-validation.hh>
# include <hpp/rbprm/rbprm-validation-report.hh>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::vector_t;
    using model::size_type;
    
    core::ConfigurationPtr_t ProjectionShooter::shoot () const
    {
      core::ConfigValidationsPtr_t configValidations
	(core::ConfigValidations::create ()); // full robot validation
      core::ValidationReportPtr_t validationReport;
      core::RbprmValidationReportPtr_t rbValidationReport = boost::dynamic_pointer_cast<core::RbprmValidationReport> (validationReport);
      RbPrmValidationPtr_t rbValidation
	(RbPrmValidation::create(rbRobot_, filter_, normalFilter_));
      core::ConfigurationPtr_t config
	(new core::Configuration_t (robot_->configSize ()));
      //bool contactValid = rbprmPathValidation->getValidator()->validateRoms(*q_new);
      
      do{
	// Sample a collision-free configuration
	do {
	  *config = uniformlySample ();
	}
	while (!configValidations->validate (*config, validationReport));

	// Project on nearest obstacle and shift away
	// TODO: shift to reach rom
	*config = project (*config);
      }while (!rbValidation->validate (*config));
      return config;
    }

    core::Configuration_t ProjectionShooter::uniformlySample () const {
      core::Configuration_t q (robot_->configSize ());
      core::JointVector_t jv = robot_->getJointVector ();
      for (core::JointVector_t::const_iterator itJoint = jv.begin ();
	   itJoint != jv.end (); itJoint++) {
	std::size_t rank = (*itJoint)->rankInConfiguration ();
	(*itJoint)->configuration ()->uniformlySample (rank, q);
      }
      // Shoot extra configuration variables
      size_type extraDim = robot_->extraConfigSpace ().dimension ();
      size_type offset = robot_->configSize () - extraDim;
      for (size_type i=0; i<extraDim; ++i) {
	value_type lower = robot_->extraConfigSpace ().lower (i);
	value_type upper = robot_->extraConfigSpace ().upper (i);
	value_type range = upper - lower;
	if ((range < 0) ||
	    (range == std::numeric_limits<double>::infinity())) {
	  std::ostringstream oss
	    ("Cannot uniformy sample extra config variable ");
	  oss << i << ". min = " <<lower<< ", max = " << upper << std::endl;
	  throw std::runtime_error (oss.str ());
	}
	q [offset + i] = lower + (upper - lower) * rand ()/RAND_MAX;
      }
      return q;
    }

    core::Configuration_t ProjectionShooter::project
    (const core::Configuration_t q) const {
      core::Configuration_t qout = q;
      const size_type ecsDim = robot_->extraConfigSpace ().dimension ();
      const size_type index = robot_->configSize() - ecsDim; // ecs index
      fcl::Vec3f pi, pj, dir, n; // fcl nearest points of collision pairs
      value_type minDistance = std::numeric_limits <value_type>::infinity();
      value_type distance = minDistance;
      core::CollisionObjectPtr_t nearestRob, nearestObst;
      core::DistanceBetweenObjectsPtr_t distanceBetweenObjects
	(problem_.distanceBetweenObjects ());
      core::ConfigValidationsPtr_t configValidations (problem_.configValidations());
      core::ValidationReportPtr_t validationReport;

      // Step 1: get nearest obstacle and surface-normale at config q
      robot_->currentConfiguration (q);
      robot_->computeForwardKinematics ();
      distanceBetweenObjects->computeDistances (); // only outers !
      const model::DistanceResults_t& dr =
	distanceBetweenObjects->distanceResults ();

      for (model::DistanceResults_t::const_iterator itDistance = 
	     dr.begin (); itDistance != dr.end (); itDistance++) {
	distance = itDistance->distance ();
	if (distance < minDistance){
	  minDistance = distance;
	  nearestRob = itDistance->innerObject;
	  nearestObst = itDistance->outerObject;
	  pi = itDistance->closestPointInner (); // point Body
	  pj = itDistance->closestPointOuter (); // point Obst
	  dir = pi - pj; // obstacle normale direction
	}
      }
      // WARN: nearestObst is not necessarily an effector of the ROM (TODO)
      const value_type dir_norm = sqrt (dir [0]*dir [0] + dir [1]*dir [1]
					+ dir [2]*dir [2]);
      n = dir/dir_norm;

      // Step 2: set orientation with n
      qout (index) = n [0];
      qout (index + 1) = n [1];
      qout (index + 2) = n [2];
      qout = setOrientation (robot_, qout); //PB linking FCL quat.fromRotation (Mat)
      if (!configValidations->validate (qout, validationReport))
	return qout; // thrown by planner

      // Step 3: re-compute new distance to nearestObst
      fcl::DistanceRequest distanceRequest (true, 0, 0, fcl::GST_INDEP);
      model::DistanceResult dr1;
      robot_->currentConfiguration (qout);
      robot_->computeForwardKinematics (); // may not be necessary
      fcl::distance (nearestRob->fcl ().get (), nearestObst->fcl ().get (),
		     distanceRequest, dr1.fcl);
      const value_type dist = dr1.distance ();

      // Step 4: project and shift
      // make shiftDist negative to have ROM COLLIDING
      const value_type shiftDist = -shiftDistance_;
      qout (0) += - dist * n [0] + shiftDistance_ * n [0]; // x
      qout (1) += - dist * n [1] + shiftDistance_ * n [1]; // y
      qout (2) += - dist * n [2] + shiftDistance_ * n [2]; // z
      return qout;
    }
    /// \}
  } //   namespace rbprm
} // namespace hpp

