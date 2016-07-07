//
// Copyright (c) 2016 CNRS
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

#ifndef HPP_RBPRM_TIMED_BALLISTIC_PATH_HH
# define HPP_RBPRM_TIMED_BALLISTIC_PATH_HH

# include <hpp/core/fwd.hh>
# include <hpp/core/config.hh>
# include <hpp/core/path.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-path.hh>

namespace hpp {
  namespace rbprm {
    
    using core::value_type;
    // forward declaration of class
    HPP_PREDEF_CLASS (TimedBallisticPath);
    // Planner objects are manipulated only via shared pointers
    typedef boost::shared_ptr <TimedBallisticPath> TimedBallisticPathPtr_t;
    
    
    /// Linear interpolation between two configurations
    ///
    /// Degrees of freedom are interpolated depending on the type of
    /// \link hpp::model::Joint joint \endlink
    /// they parameterize:
    ///   \li linear interpolation for translation joints, bounded rotation
    ///       joints, and translation part of freeflyer joints,
    ///   \li angular interpolation for unbounded rotation joints,
    ///   \li constant angular velocity for SO(3) part of freeflyer joints.
    class TimedBallisticPath : public core::Path
    {
    public:
      typedef Path parent_t;
      /// Destructor
      virtual ~TimedBallisticPath () throw () {}
      
      /// Create instance and return shared pointer
      /// \param device Robot corresponding to configurations
      /// \param init, end Start and end configurations of the path
      /// \param length Distance between the configurations.
      static TimedBallisticPathPtr_t create (const core::DevicePtr_t& device,
                                             core::ConfigurationIn_t init,
                                             core::ConfigurationIn_t end,
                                             core::value_type length,
                                             value_type alpha,
                                             value_type theta,
                                             value_type v0)
      {
        TimedBallisticPath* ptr = new TimedBallisticPath (device, init, end, length,
                                                           alpha, theta, v0);
        TimedBallisticPathPtr_t shPtr (ptr);
        ptr->init (shPtr);
        return shPtr;
      }
      
      static TimedBallisticPathPtr_t create (const rbprm::BallisticPathPtr_t ballisticPath)
      {
        TimedBallisticPath* ptr = new TimedBallisticPath (ballisticPath);
        TimedBallisticPathPtr_t shPtr (ptr);
        ptr->init (shPtr);
        return shPtr;
      }
      
      static TimedBallisticPathPtr_t create (const BallisticPathPtr_t bp1,const BallisticPathPtr_t bp1Max,const BallisticPathPtr_t bp2Max,const BallisticPathPtr_t bp2)
      {
        TimedBallisticPath* ptr = new TimedBallisticPath (bp1,bp1Max,bp2Max,bp2);
        TimedBallisticPathPtr_t shPtr (ptr);
        ptr->init (shPtr);
        return shPtr;
      }
      
      /// Create copy and return shared pointer
      /// \param path path to copy
      static TimedBallisticPathPtr_t createCopy (const TimedBallisticPathPtr_t& path)
      {
        TimedBallisticPath* ptr = new TimedBallisticPath (*path);
        TimedBallisticPathPtr_t shPtr (ptr);
        ptr->initCopy (shPtr);
        return shPtr;
      }
      
      /// Create copy and return shared pointer
      /// \param path path to copy
      /// \param constraints the path is subject to
      /// <!> constraints part NOT IMPLEMENTED YET
      static TimedBallisticPathPtr_t createCopy
      (const TimedBallisticPathPtr_t& path, const core::ConstraintSetPtr_t& /*constraints*/)
      {
        //TimedBallisticPath* ptr = new TimedBallisticPath (*path, constraints);
        TimedBallisticPath* ptr = new TimedBallisticPath (*path);
        TimedBallisticPathPtr_t shPtr (ptr);
        ptr->initCopy (shPtr);
        return shPtr;
      }
      
      /// Return a shared pointer to this
      ///
      /// As TimedBallisticPath are immutable, and refered to by shared pointers,
      /// they do not need to be copied.
      virtual core::PathPtr_t copy () const
      {
        return createCopy (weak_.lock ());
      }
      
      /// Return a shared pointer to a copy of this and set constraints
      ///
      /// \param constraints constraints to apply to the copy
      /// \precond *this should not have constraints.
      virtual core::PathPtr_t copy (const core::ConstraintSetPtr_t& constraints) const
      {
        return createCopy (weak_.lock (), constraints);
      }
      
      /// Extraction/Reversion of a sub-path
      /// \param subInterval interval of definition of the extract path
      /// If upper bound of subInterval is smaller than lower bound,
      /// result is reversed.
      virtual core::PathPtr_t extract (const core::interval_t& subInterval) const throw (core::projection_error);
      
      /// Reversion of a path
      /// \return a new path that is this one reversed.
      virtual core::PathPtr_t reverse () const;
      
      
      /// Modify initial configuration
      /// \param initial new initial configuration
      /// \pre input configuration should be of the same size as current initial
      /// configuration
      void initialConfig (core::ConfigurationIn_t initial)
      {
        assert (initial.size () == initial_.size ());
        initial_ = initial;
      }
      
      /// Modify end configuration
      /// \param end new end configuration
      /// \pre input configuration should be of the same size as current end
      /// configuration
      void endConfig (core::ConfigurationIn_t end)
      {
        assert (end.size () == end_.size ());
        end_ = end;
      }
      
      /// Return the internal robot.
      core::DevicePtr_t device () const;
      
      /// Get the initial configuration
      core::Configuration_t initial () const
      {
        return initial_;
      }
      
      /// Get the final configuration
      core::Configuration_t end () const
      {
        return end_;
      }
      
      /// Get previously computed length
      virtual core::value_type length () const {
        return length_;
      }
      

      

      
      core::value_type computeLength (const core::ConfigurationIn_t q1,
                                      const core::ConfigurationIn_t q2) const;
      
      /// Evaluate velocity vector at path abcissa t
      
    protected:
      /// Print path in a stream
      virtual std::ostream& print (std::ostream &os) const
      {
        os << "TimedBallisticPath:" << std::endl;
        os << "interval: [ " << timeRange ().first << ", "
           << timeRange ().second << " ]" << std::endl;
        os << "initial configuration: " << initial_.transpose () << std::endl;
        os << "final configuration:   " << end_.transpose () << std::endl;
        return os;
      }
      /// Constructor
      TimedBallisticPath (const core::DevicePtr_t& robot,
                          core::ConfigurationIn_t init,
                          core::ConfigurationIn_t end, core::value_type length,
                          value_type alpha,
                          value_type theta,
                          value_type v0);
      
      
      TimedBallisticPath (const BallisticPathPtr_t ballisticPath);
      
      TimedBallisticPath (const BallisticPathPtr_t bp1,const BallisticPathPtr_t bp1Max,const BallisticPathPtr_t bp2Max,const BallisticPathPtr_t bp2);
      
      
      /// Copy constructor
      TimedBallisticPath (const TimedBallisticPath& TimedBallisticPath);
      
      
      void init (TimedBallisticPathPtr_t self)
      {
        parent_t::init (self);
        weak_ = self;
      }
      
      void initCopy (TimedBallisticPathPtr_t self)
      {
        parent_t::initCopy (self);
        weak_ = self;
      }
      
      /// Same as parabola-path for trunk (freeflyer + internal trunk DOF)
      /// Param is the curvilinear abcissa \in [0 : pathLength]
      /// The pathLength can be computed as long as the coefficients_ are known
      /// Finally:
      /// config(0) = x(param) = (1 - param/length)*x1 + param/length*x2
      /// config(1) = coefs(0)*x(param)^2 + coefs(1)*x(param) + coefs(2)
      virtual bool impl_compute (core::ConfigurationOut_t result,
                                 core::value_type t) const;
      
    private:
      core::DevicePtr_t device_;
      core::Configuration_t initial_;
      core::Configuration_t end_;
      TimedBallisticPathWkPtr_t weak_;
      mutable core::value_type length_;      
      value_type alpha_;
      value_type theta_;
      value_type v0_;
      value_type g_;
      value_type xTheta0_;
      // we use the 4 subpath for joint interpolation
      BallisticPathPtr_t bp1_;
      BallisticPathPtr_t bp1Max_;
      BallisticPathPtr_t bp2Max_;
      BallisticPathPtr_t bp2_;
      value_type t1_,t1Max_,t2Max_,t2_;
      
    }; // class TimedBallisticPath
  } //   namespace rbprm
} // namespace hpp
#endif // HPP_CORE_BALLISTIC_PATH_HH
