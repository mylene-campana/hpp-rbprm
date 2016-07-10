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
#include <hpp/rbprm/fullbodyBallistic/timed-ballistic-path.hh>


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
                                          const std::size_t pathDofRank,const core::PathPtr_t rootPath)
      {
    LimbRRTSteering* ptr = new LimbRRTSteering (problem, pathDofRank,rootPath);
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
        else if(tbp_){
            length = tbp_->computeLength(q1,q2);
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
          BallisticPathPtr_t bpExtract =  BallisticPath::create(bp_->device(),q1,q2,length,bp_->coefficients(),bp_->lastRootIndex());
          path = LimbRRTPath::create
            (problem_->robot(), q1, q2, length, c, pathDofRank_,bpExtract);
        }
        else if(tbp_){
          hppDout(notice,"create path with ballistic root");
	      length = tbp_->ballisticPath()->computeLength(q1,q2);
          BallisticPathPtr_t bpExtract =  BallisticPath::create(tbp_->device(),q1,q2,length,tbp_->coefficients(),tbp_->lastRootIndex());
          TimedBallisticPathPtr_t tbpExtract =  TimedBallisticPath::create(bpExtract);
 		  length = tbpExtract->length();
          size_t rankParamRoot = q1.size() - 1;
          core::Configuration_t q11(q1);
          q11[rankParamRoot] = 0.;
          core::Configuration_t q22(q2);
          q22[rankParamRoot] = length;          
	      path = LimbRRTPath::create
            (problem_->robot(), q1, q2, length, c, pathDofRank_,tbpExtract);
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
                       const std::size_t pathDofRank,core::PathPtr_t root) :
    SteeringMethod (problem), pathDofRank_(pathDofRank), weak_ (),bp_(),tbp_()
      {
          BallisticPathPtr_t bp = boost::dynamic_pointer_cast<BallisticPath>(root);
          if(bp)
              bp_=bp;

          TimedBallisticPathPtr_t tbp = boost::dynamic_pointer_cast<TimedBallisticPath>(root);
          if(tbp)
              tbp_=tbp;

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
      BallisticPathPtr_t bp_;
      TimedBallisticPathPtr_t tbp_;
    }; // SteeringMethodStraight
    /// \}
  } // namespace interpolation
  } // namespace rbprm
} // namespace hpp

#endif // HPP_RBPRM_LIMB_RRT_STEERING_HH
