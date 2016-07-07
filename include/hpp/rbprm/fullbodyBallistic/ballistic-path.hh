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

#ifndef HPP_RBPRM_BALLISTIC_PATH_HH
# define HPP_RBPRM_BALLISTIC_PATH_HH

# include <hpp/core/fwd.hh>
# include <hpp/core/config.hh>
# include <hpp/core/path.hh>

namespace hpp {
  namespace rbprm {


    // forward declaration of class
    HPP_PREDEF_CLASS (BallisticPath);
    // Planner objects are manipulated only via shared pointers
    typedef boost::shared_ptr <BallisticPath> BallisticPathPtr_t;


    /// Linear interpolation between two configurations
    ///
    /// Degrees of freedom are interpolated depending on the type of
    /// \link hpp::model::Joint joint \endlink
    /// they parameterize:
    ///   \li linear interpolation for translation joints, bounded rotation
    ///       joints, and translation part of freeflyer joints,
    ///   \li angular interpolation for unbounded rotation joints,
    ///   \li constant angular velocity for SO(3) part of freeflyer joints.
    class BallisticPath : public core::Path
    {
    public:
      typedef Path parent_t;
      /// Destructor
      virtual ~BallisticPath () throw () {}

      /// Create instance and return shared pointer
      /// \param device Robot corresponding to configurations
      /// \param init, end Start and end configurations of the path
      /// \param length Distance between the configurations.
      static BallisticPathPtr_t create (const core::DevicePtr_t& device,
                                       core::ConfigurationIn_t init,
                                       core::ConfigurationIn_t end,
                                       core::value_type length,
                                       core::vector_t coefficients)
      {
	BallisticPath* ptr = new BallisticPath (device, init, end, length,
					      coefficients);
	BallisticPathPtr_t shPtr (ptr);
	ptr->init (shPtr);
	return shPtr;
      }

      /// Create copy and return shared pointer
      /// \param path path to copy
      static BallisticPathPtr_t createCopy (const BallisticPathPtr_t& path)
      {
	BallisticPath* ptr = new BallisticPath (*path);
	BallisticPathPtr_t shPtr (ptr);
	ptr->initCopy (shPtr);
	return shPtr;
      }

      /// Create copy and return shared pointer
      /// \param path path to copy
      /// \param constraints the path is subject to
      /// <!> constraints part NOT IMPLEMENTED YET
      static BallisticPathPtr_t createCopy
        (const BallisticPathPtr_t& path, const core::ConstraintSetPtr_t& /*constraints*/)
      {
	//BallisticPath* ptr = new BallisticPath (*path, constraints);
        BallisticPath* ptr = new BallisticPath (*path);
	BallisticPathPtr_t shPtr (ptr);
	ptr->initCopy (shPtr);
	return shPtr;
      }

      /// Return a shared pointer to this
      ///
      /// As BallisticPath are immutable, and refered to by shared pointers,
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

      /// Set the three parabola coefficients
      void coefficients (core::vector_t coefs) const {
	for (size_t i = 0; i < coefs.size () ; i++)
	  coefficients_ (i) = coefs (i);
      }

      /// Get parabola coefficients
      core::vector_t coefficients () const {
	return coefficients_;
      }

      core::value_type computeLength (const core::ConfigurationIn_t q1,
				      const core::ConfigurationIn_t q2) const;

      /// Evaluate velocity vector at path abcissa t
      core::vector_t evaluateVelocity (const core::value_type t) const;

    protected:
      /// Print path in a stream
      virtual std::ostream& print (std::ostream &os) const
      {
	os << "BallisticPath:" << std::endl;
	os << "interval: [ " << timeRange ().first << ", "
	   << timeRange ().second << " ]" << std::endl;
	os << "initial configuration: " << initial_.transpose () << std::endl;
	os << "final configuration:   " << end_.transpose () << std::endl;
	return os;
      }
      /// Constructor
      BallisticPath (const core::DevicePtr_t& robot,
		     core::ConfigurationIn_t init,
		     core::ConfigurationIn_t end, core::value_type length,
		     core::vector_t coefficients);

      /// Copy constructor
      BallisticPath (const BallisticPath& Ballisticpath);

      core::value_type lengthFunction (const core::value_type x)const;

      void init (BallisticPathPtr_t self)
      {
	parent_t::init (self);
	weak_ = self;
      }

      void initCopy (BallisticPathPtr_t self)
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
                                 core::value_type param) const;

    private:
      core::DevicePtr_t device_;
      core::Configuration_t initial_;
      core::Configuration_t end_;
      BallisticPathWkPtr_t weak_;
      mutable core::vector_t coefficients_; // 4 parabola coefficients
      mutable core::value_type length_;
    }; // class BallisticPath
  } //   namespace rbprm
} // namespace hpp
#endif // HPP_CORE_BALLISTIC_PATH_HH
