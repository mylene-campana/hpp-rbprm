//
// Copyright (c) 2015 CNRS
// Authors: Florent Lamiraux
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

#ifndef HPP_RBPRM_LIMB_RRT_STEERING_HH
# define HPP_RBPRM_LIMB_RRT_STEERING_HH

# include <hpp/core/discretized-path-validation.hh>
# include <hpp/core/steering-method-straight.hh>
# include <hpp/core/straight-path.hh>
# include <hpp/rbprm/interpolation/limb-rrt-path.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-path.hh>

namespace hpp {
  namespace rbprm {
  namespace interpolation {
    HPP_PREDEF_CLASS(LimbRRTSteering);

    class LimbRRTSteering;
    typedef boost::shared_ptr <LimbRRTSteering> LimbRRTSteeringPtr_t;
    /// \addtogroup validation
    /// \{


    /// Discretized validation of a path for the LimbRRT algorithm
    ///
    /// Apply some configuration validation algorithms at discretized values
    /// of the path parameter.
    class HPP_CORE_DLLAPI LimbRRTSteering : public hpp::core::SteeringMethod
    {
    public:
      /// Create instance and return shared pointer
      static LimbRRTSteeringPtr_t create (const core::ProblemPtr_t& problem,
                                          const std::size_t pathDofRank,const BallisticPathPtr_t bp)
      {
    LimbRRTSteering* ptr = new LimbRRTSteering (problem, pathDofRank,bp);
    LimbRRTSteeringPtr_t shPtr (ptr);
    ptr->init (shPtr);
    return shPtr;
      }
      
      static LimbRRTSteeringPtr_t create (const core::ProblemPtr_t& problem,
                                          const std::size_t pathDofRank)
      {
    LimbRRTSteering* ptr = new LimbRRTSteering (problem, pathDofRank);
    LimbRRTSteeringPtr_t shPtr (ptr);
    ptr->init (shPtr);
    return shPtr;
      }
      
      
      /// Create instance and return shared pointer
      static LimbRRTSteeringPtr_t create
    (const core::DevicePtr_t& device, const core::WeighedDistancePtr_t& distance,
      const std::size_t pathDofRank)
        HPP_CORE_DEPRECATED
      {
    LimbRRTSteering* ptr = new LimbRRTSteering (device,
                                  distance, pathDofRank);
    LimbRRTSteeringPtr_t shPtr (ptr);
    ptr->init (shPtr);
    return shPtr;
      }
      /// Copy instance and return shared pointer
      static LimbRRTSteeringPtr_t createCopy
    (const LimbRRTSteeringPtr_t& other)
      {
    LimbRRTSteering* ptr = new LimbRRTSteering (*other);
    LimbRRTSteeringPtr_t shPtr (ptr);
    ptr->init (shPtr);
    return shPtr;
      }
      /// Copy instance and return shared pointer
      virtual core::SteeringMethodPtr_t copy () const
      {
    return createCopy (weak_.lock ());
      }

      /// create a path between two configurations
      virtual core::PathPtr_t impl_compute (core::ConfigurationIn_t q1,
                      core::ConfigurationIn_t q2) const
      {

        core::value_type length;
        if(bp_){
            length = bp_->computeLength(q1,q2);
        }
        else{
            length = (*problem_->distance()) (q1, q2);
        }

        core::ConstraintSetPtr_t c;
        if (constraints() && constraints()->configProjector ()) {
          c = HPP_STATIC_PTR_CAST (core::ConstraintSet, constraints()->copy ());
          c->configProjector()->rightHandSideFromConfig (q1);
        } else {
          c = constraints ();
        }
        core::PathPtr_t path;
        if(bp_){
          size_t rankParamRoot = q1.size() - 1;
          core::Configuration_t q11(q1);
          q11[rankParamRoot] = 0.;
          core::Configuration_t q22(q2);
          q22[rankParamRoot] = length;
          hppDout(notice,"create path with ballistic root");
          BallisticPathPtr_t bpExtract =  BallisticPath::create(bp_->device(),q11,q22,length,bp_->coefficients());
          bpExtract->lastRootIndex(bp_->lastRootIndex());
          path = LimbRRTPath::create
            (problem_->robot(), q11, q22, length, c, pathDofRank_,bpExtract);
        }
        else{
          hppDout(notice,"create path without root path");
          path = LimbRRTPath::create
            (problem_->robot(), q1, q2, length, c, pathDofRank_);
        }
            return path;
      }
    protected:
      /// Constructor with robot
      /// Weighed distance is created from robot
      LimbRRTSteering (const core::ProblemPtr_t& problem,
                       const std::size_t pathDofRank,BallisticPathPtr_t bp) :
    SteeringMethod (problem), pathDofRank_(pathDofRank), weak_ (),bp_(bp)
      {
      }
      
      /// Constructor with robot
      /// Weighed distance is created from robot
      LimbRRTSteering (const core::ProblemPtr_t& problem,
                       const std::size_t pathDofRank) :
    SteeringMethod (problem), pathDofRank_(pathDofRank), weak_ (),bp_()
      {
      }
      
      
      /// Constructor with weighed distance
      LimbRRTSteering (const core::DevicePtr_t& device,
                  const core::WeighedDistancePtr_t& distance,
                       const std::size_t pathDofRank) :
    SteeringMethod (new core::Problem (device)), pathDofRank_(pathDofRank), weak_ ()
      {
        problem_->distance (distance);
      }
      /// Copy constructor
      LimbRRTSteering (const LimbRRTSteering& other) :
    SteeringMethod (other), pathDofRank_(other.pathDofRank_), weak_ ()
      {
      }

      /// Store weak pointer to itself
      void init (LimbRRTSteeringWkPtr_t weak)
      {
    SteeringMethod::init (weak);
    weak_ = weak;
      }

    private:      
      const std::size_t pathDofRank_;
      LimbRRTSteeringWkPtr_t weak_;
      const BallisticPathPtr_t bp_;
    }; // SteeringMethodStraight
    /// \}
  } // namespace interpolation
  } // namespace rbprm
} // namespace hpp

#endif // HPP_RBPRM_LIMB_RRT_STEERING_HH
