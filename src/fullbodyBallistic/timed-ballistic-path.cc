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
#include <hpp/rbprm/fullbodyBallistic/timed-ballistic-path.hh>
#include <hpp/core/straight-path.hh>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::value_type;
    using core::vector_t;
    using core::interval_t;
    using core::size_type;
    
    TimedBallisticPath::TimedBallisticPath (const core::DevicePtr_t& device,
                                            core::ConfigurationIn_t init,
                                            core::ConfigurationIn_t end,
                                            value_type length,
                                            value_type theta,
					    value_type xTheta0dot) :
      parent_t (interval_t (0, length), device->configSize (),
                device->numberDof ()), device_ (device), initial_ (init),
      end_ (end), length_ (length), theta_(theta), xTheta0dot_ (xTheta0dot),
      g_(9.81)
    {
      assert (device);
    }
    
    TimedBallisticPath::TimedBallisticPath
    (const rbprm::BallisticPathPtr_t ballisticPath) :
      parent_t (interval_t(0,ballisticPath->length()),
		ballisticPath->device()->configSize(),
		ballisticPath->device()->numberDof ()),
      device_(ballisticPath->device()), initial_(ballisticPath->initial()),
      end_(ballisticPath->end()), length_(ballisticPath->length()), g_(9.81)
    {
      theta_ = ballisticPath->coefficients()[3];
      xTheta0dot_ = ballisticPath->coefficients()[5];
      length_ = computeLength(initial_,end_); // replace by duration
      times_.resize (1);
      times_ [0] = length_;
      ballisticPath_ = ballisticPath;
    }
    
    TimedBallisticPath::TimedBallisticPath (PathVectorBP pathPair) :
      parent_t (interval_t (0,pathPair.first->length()),
		pathPair.second->device()->configSize(),
		pathPair.second->device()->numberDof ()),
      device_(pathPair.second->device()),
      initial_(pathPair.first->initial()),
      end_(pathPair.first->end()), length_(),
      g_(9.81), pathVector_ (pathPair.first),
      ballisticPath_ (pathPair.second)
    {
      theta_ = ballisticPath_->coefficients()[3];
      xTheta0dot_ = ballisticPath_->coefficients()[5];
      length_ = computeLength(initial_,end_);
      hppDout (info, "length_= " << length_);
      times_.resize (pathVector_->numberPaths ());
      hppDout (info, "number of subpath in pathVector= " << pathVector_->numberPaths ());
      
      // compute and store duration of each subpath:
      value_type XTheta_0 = cos(theta_)*initial_[0] + sin(theta_)*initial_[1];
      hppDout (info, "XTheta_0= " << XTheta_0);
      hppDout (info, "previous config= " << displayConfig(initial_));
      for (std::size_t i = 0; i < pathVector_->numberPaths () ; i++) {
	const core::PathPtr_t path_i = pathVector_->pathAtRank (i);
	const value_type XTheta_i = cos(theta_)*path_i->end()[0] + sin(theta_)*path_i->end()[1] - XTheta_0;
	hppDout (info, "config_i= " << displayConfig(path_i->end()));
	times_ [i] = XTheta_i / xTheta0dot_;
	hppDout (info, "i= " << i);
	hppDout (info, "XTheta_i= " << XTheta_i);
	hppDout (info, "times_ [i]= " << times_ [i]);
      }
    }
    
    TimedBallisticPath::TimedBallisticPath (const TimedBallisticPath& path) :
      parent_t (path), device_ (path.device_), initial_ (path.initial_),
      end_ (path.end_),length_ (path.length_),theta_(path.theta_),
      xTheta0dot_ (path.xTheta0dot_), g_ (path.g_),
      pathVector_ (path.pathVector_),
      ballisticPath_ (path.ballisticPath_), times_ (path.times_)
    {
    }
    
    bool TimedBallisticPath::impl_compute (core::ConfigurationOut_t result,
                                           value_type t) const
    {
      if (t == 0 || initial_(0) == end_(0)) {
        result = initial_;
        return true;
      }
      if (t >= length_) {
        result = end_;
        return true;
      }

      // convertion of BP param u into TBP param t, then
      // interpolation for the joints are done by ballisticPath
      std::size_t i = 0;
      while (t > times_ [i] && i < times_.size () - 1) i++;
      //hppDout (info, "t= " << t << ", times_["<<i<<"]= " << times_ [i]);
      // t < times_[i]
      value_type u;
      const core::PathPtr_t pv_i = pathVector_->pathAtRank (i);
      if (i == 0)
	u = t/times_ [i] * pv_i->length();
      else
	u = (t - times_[i-1])/(times_[i] - times_[i-1]) * pv_i->length();
      (*pv_i)(result,u);
      // hppDout (info, "u= " << u);
      //hppDout (info, "tmp config= " << displayConfig(result));
      
      // replace with the correct position / orientation for the center
      //const value_type xTheta0 =cos(theta_)*initial[0]+sin(theta_)*initial[1];
      const value_type alpha = ballisticPath_->coefficients()[4];
      const value_type xTheta = xTheta0dot_ * t;
      result[0] = xTheta*cos(theta_) + initial_[0];
      result[1] = xTheta*sin(theta_) + initial_[1];     
      result[2] = -0.5*g_*t*t + xTheta0dot_*tan(alpha)*t + initial_[2];

      /* Quaternions interpolation */
      u = t/length_;
      const core::JointPtr_t SO3joint = device_->getJointByName ("base_joint_SO3");
      const std::size_t rank = SO3joint->rankInConfiguration ();
      const core::size_type dimSO3 = SO3joint->configSize ();
      SO3joint->configuration ()->interpolate (initial_, end_, u, rank, result);
      return true;
    }
    
    core::PathPtr_t TimedBallisticPath::extract (const interval_t& subInterval) const throw (hpp::core::projection_error)
    {
      bool success;
      core::Configuration_t q1 ((*this) (subInterval.first, success)); // straight
      core::Configuration_t q2 ((*this) (subInterval.second, success)); // straight
      core::PathPtr_t result = rbprm::TimedBallisticPath::create(device_, q1, q2, computeLength(q1,q2), theta_, xTheta0dot_);
      return result;
    }
    
    core::PathPtr_t TimedBallisticPath::reverse () const{
      hppDout(notice, "reverse timed-ballistic path");
      bool success;
      core::Configuration_t q1 ((*this) (length_, success));
      core::Configuration_t q2 ((*this) (0, success));
      core::PathPtr_t result = TimedBallisticPath::create (device_, q1, q2, length_, theta_, -xTheta0dot_);
      return result;
    }
    
    // duration
    value_type TimedBallisticPath::computeLength
    (const core::ConfigurationIn_t q1, const core::ConfigurationIn_t q2) const {
      const value_type XTheta = cos(theta_)*(q2[0] - q1[0]) + sin(theta_)*(q2[1] - q1[1]);
      const value_type lenght = XTheta / xTheta0dot_;
      return lenght; // path duration
    }
    
    
    
  } //   namespace rbprm
} // namespace hpp
