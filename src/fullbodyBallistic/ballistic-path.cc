//
// Copyright (c) 2014 CNRS
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

#include <hpp/util/debug.hh>
#include <hpp/model/configuration.hh>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/joint-configuration.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-path.hh>
#include <hpp/core/straight-path.hh>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::value_type;
    using core::vector_t;
    using core::interval_t;
    using core::size_type;

    BallisticPath::BallisticPath (const core::DevicePtr_t& device,
				  core::ConfigurationIn_t init,
				  core::ConfigurationIn_t end,
				  value_type length,
				  vector_t coefs) :
      parent_t (interval_t (0, length), device->configSize (),
                device->numberDof ()), device_ (device), initial_ (init),
      end_ (end), coefficients_ (vector_t(coefs.size ())), length_ (length)
    {
      assert (device);
      coefficients (coefs);
    }

    BallisticPath::BallisticPath (const BallisticPath& path) :
      parent_t (path), device_ (path.device_), initial_ (path.initial_),
      end_ (path.end_), coefficients_ (path.coefficients_),
      length_ (path.length_)
    {
    }

    bool BallisticPath::impl_compute (core::ConfigurationOut_t result,
				      value_type param) const
    {
      if (param == 0 || initial_(0) == end_(0)) {
        result = initial_;
        return true;
      }
      if (param >= length_) {
        result = end_;
        return true;
      }

      const size_type nbConfig = device_->configSize();
      const size_type ecsDim = device_->extraConfigSpace ().dimension ();
      const value_type u = param/length_;
      const value_type theta = coefficients_(3);
      const value_type x_theta_max = - 0.5 *
	coefficients_ (1) / coefficients_ (0);
      const value_type x_theta_initial = cos(theta)*initial_ (0) +
	sin(theta)*initial_ (1);
      const value_type x_theta_end = cos(theta)*end_ (0) +
	sin(theta)*end_ (1);
      const value_type u_max = (x_theta_max - x_theta_initial)
	/ (x_theta_end - x_theta_initial);
      const bool tanThetaNotDefined = (theta < M_PI/2 + 1e-2 && theta > M_PI/2 - 1e-2) || (theta > -M_PI/2 - 1e-2 && theta < -M_PI/2 + 1e-2);

      if (!tanThetaNotDefined) { //theta != +- pi/2
	const value_type tanTheta = tan(theta);
	result (0) = (1 - u)*initial_(0) + u*end_(0);
	result (1) = tanTheta*result (0) -tanTheta*initial_(0) + initial_(1);
	const value_type x_theta = cos(theta)*result (0) +
	  sin(theta)*result (1);
	result (2) = coefficients_(0)*x_theta*x_theta
	  + coefficients_(1)*x_theta + coefficients_(2);
      }
      else { //theta = +- pi/2
	result (0) = initial_ (0);
	result (1) = (1 - u)*initial_(1) + u*end_(1);
	const value_type x_theta = cos(theta)*result (0) +
	  sin(theta)*result (1);
	result (2) = coefficients_(0)*x_theta*x_theta
	  + coefficients_(1)*x_theta + coefficients_(2);
      }

      /* Quaternions interpolation */
      const core::JointPtr_t SO3joint = device_->getJointByName ("base_joint_SO3");
      const std::size_t rank = SO3joint->rankInConfiguration ();
      const core::size_type dimSO3 = SO3joint->configSize ();
      SO3joint->configuration ()->interpolate
	(initial_, end_, u, rank, result);

      /* if robot has internal DoF (except freeflyer ones) */
      // translation dimension of freeflyer hardcoded...
      // min value (to reach for u = u_max) hardcoded...
      // manual interpolation since joint not available with index...
      const std::size_t freeflyerDim = 3 + dimSO3;
      const bool hasInternalDof = nbConfig > ecsDim + freeflyerDim;
      const value_type maxVal = 0; // because here initial_ = end_ ...
      if (hasInternalDof) {
	for (core::size_type i = freeflyerDim; i<nbConfig-ecsDim; i++) {
	  /* monopod leg interpolation
	    if (u <= u_max) {
	    const value_type u_prime = u / u_max;
	    result (i) = (1 - u_prime) * initial_ (i) + u_prime * maxVal;
	  }
	  else {
	    const value_type u_prime = (u - u_max) / (1 - u_max);
	    result (i) = (1 - u_prime) * maxVal + u_prime * end_ (i);
	    }*/
	  /* classical interpolation for robot trunk and limbs */
	  result (i) = (1 - u) * initial_ (i) + u * end_ (i);
	}
      }

      /* Normal vector interpolation
	 result (nbConfig-ecsDim) = (1 - u) *
	 initial_(nbConfig-ecsDim) + u*end_(nbConfig-ecsDim);
	 result (nbConfig-ecsDim+1) = (1 - u) *
	 initial_(nbConfig-ecsDim+1) + u*end_(nbConfig-ecsDim+1);
	 result (nbConfig-ecsDim+2) = (1 - u) *
	 initial_(nbConfig-ecsDim+2) + u*end_(nbConfig-ecsDim+2);*/
      return true;
    }


    core::PathPtr_t BallisticPath::extract (const interval_t& subInterval) const throw (hpp::core::projection_error)
    {
      bool success;
      core::Configuration_t q1 ((*this) (subInterval.first, success)); // straight
      core::Configuration_t q2 ((*this) (subInterval.second, success)); // straight
      core::PathPtr_t result = rbprm::BallisticPath::create(device_,q1,q2,computeLength(q1,q2),coefficients_);
      return result;
    }

    core::PathPtr_t BallisticPath::reverse () const{
      hppDout(notice, "reverse path parabola !!!!!!!!!!!!!!!!!!!!!!!!");
      bool success;
      core::Configuration_t q1 ((*this) (length_, success));
      core::Configuration_t q2 ((*this) (0, success));
      core::PathPtr_t result = BallisticPath::create (device_, q1, q2, length_,
                                                     coefficients_);
      return result;
    }

    core::DevicePtr_t BallisticPath::device () const
    {
      return device_;
    }

    value_type BallisticPath::computeLength
    (const core::ConfigurationIn_t q1, const core::ConfigurationIn_t q2) const {
      const int N = 6; // number -1 of interval sub-divisions
      // for N = 4, computation error ~= 1e-5.
      // for N = 20, computation error ~= 1e-11.
      value_type length = 0;
      const value_type theta = coefficients_ (3);
      value_type x1 = cos(theta) * q1 (0)  + sin(theta) * q1 (1); // x_theta_0
      value_type x2 = cos(theta) * q2 (0) + sin(theta) * q2 (1); // x_theta_imp
      hppDout(notice,"xTheta0 = "<<x1<<"   xThetaImp = "<<x2);
      // Define integration bounds
      if (x1 > x2) { // re-order integration bounds
        const value_type xtmp = x1;
        x1 = x2;
        x2 = xtmp;
      }

      const value_type dx = (x2 - x1) / N; // integration step size
      for (int i=0; i<N; i++) {
        length += dx*( 0.166666667*lengthFunction (x1 + i*dx)
                       + 0.666666667*lengthFunction (x1 + (i+0.5)*dx )
                       + 0.166666667*lengthFunction (x1 + (i+1)*dx ));
        // apparently, 1/6 and 2/3 are not recognized as floats ...
      }
      hppDout (notice, "length = " << length);
      return length;
    }

    // Function equivalent to sqrt( 1 + f'(x)^2 )
    value_type BallisticPath::lengthFunction (const value_type x)
    const {
      const value_type y = sqrt (1+(2*coefficients_ (0)*x+coefficients_(1))
                                 * (2*coefficients_ (0)*x+coefficients_(1)));
      return y;
    }

    vector_t BallisticPath::evaluateVelocity (const value_type t) const {
      vector_t vel (3);
      bool success;
      const value_type theta = coefficients_(3);
      const value_type alpha = coefficients_(4);
      const value_type x_theta_0_dot = coefficients_(5);
      const value_type inv_x_theta_0_dot_sq = 1/(x_theta_0_dot*x_theta_0_dot);
      const value_type x_theta_0 = coefficients_(6);
      const core::Configuration_t q = (*this) (t, success);
      const value_type x_theta = q [0]*cos(theta) + q [1]*sin(theta);
      vel [0] = x_theta_0_dot * cos(theta);
      vel [1] = x_theta_0_dot * sin(theta);
      vel [2] = x_theta_0_dot * (-9.81 * (x_theta - x_theta_0)*inv_x_theta_0_dot_sq + tan(alpha));
      return vel;
    }

  } //   namespace rbprm
} // namespace hpp

