//crabe8pince
// Copyright (c) 2015-2016 CNRS
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

#ifndef HPP_RBPRM_PROJECTION_SHOOTER_HH
# define HPP_RBPRM_PROJECTION_SHOOTER_HH

# include <sstream>
# include <hpp/model/device.hh>
# include <hpp/core/problem.hh>
# include <hpp/core/configuration-shooter.hh>
# include <hpp/rbprm/rbprm-rom-validation.hh>
# include <hpp/rbprm/rbprm-device.hh>
# include <hpp/rbprm/config.hh>

namespace hpp {
  namespace rbprm {
    using core::value_type;
    /// \addtogroup configuration_sampling
    /// \{

    // forward declaration of class
    HPP_PREDEF_CLASS (ProjectionShooter);
    // Planner objects are manipulated only via shared pointers
    typedef boost::shared_ptr <ProjectionShooter> ProjectionShooterPtr_t;

    /// Uniformly sample with bounds of degrees of freedom.
    class HPP_CORE_DLLAPI ProjectionShooter : public core::ConfigurationShooter
    {
    public:
      static ProjectionShooterPtr_t 
	create (const core::DevicePtr_t& robot, const core::Problem &problem,
		const value_type shiftDist,
		const model::RbPrmDevicePtr_t& rbRobot,
		const std::vector<std::string>& filter,
		const std::map<std::string, rbprm::NormalFilter>& normalFilter)
      {
	ProjectionShooter* ptr = new ProjectionShooter (robot, problem, 
							shiftDist, rbRobot,
							filter, normalFilter);
	ProjectionShooterPtr_t shPtr (ptr);
	ptr->init (shPtr);
	return shPtr;
      }

      /// Process successively an uniform sample and a projection
      virtual core::ConfigurationPtr_t shoot () const;

      /// Uniformly sample a configuration
      core::Configuration_t uniformlySample () const;

      /// Project the given configuration on the nearest obstacle,
      /// then shift it to avoid contact
      core::Configuration_t project (const core::Configuration_t q) const;
      
    protected:
      /// Uniformly sample configuration space, then project robot on 
      /// nearest obstacle, at given distance.
      ///
      /// Note that translation joints have to be bounded.
      ProjectionShooter (const core::DevicePtr_t& robot,
			 const core::Problem &problem,
			 const value_type& shiftDist,
			 const model::RbPrmDevicePtr_t& rbRobot,
			 const std::vector<std::string>& filter,
			 const std::map<std::string, rbprm::NormalFilter>& normalFilter) :
	problem_ (problem), robot_ (robot), shiftDistance_ (shiftDist),
	rbRobot_ (rbRobot), filter_ (filter), normalFilter_ (normalFilter)
	
      {
      }
      void init (const ProjectionShooterPtr_t& self)
      {
	core::ConfigurationShooter::init (self);
	weak_ = self;
      }

    private:
      const core::Problem& problem_;
      const core::DevicePtr_t& robot_;
      ProjectionShooterWkPtr_t weak_;
      const value_type& shiftDistance_;
      const model::RbPrmDevicePtr_t& rbRobot_;
      const std::vector<std::string>& filter_;
      const std::map<std::string, rbprm::NormalFilter>& normalFilter_;
    }; // class ProjectionShooter
    /// \}
  } //   namespace rbprm
} // namespace hpp

#endif // HPP_RBPRM_CONFIGURATION_PROJECTION_SHOOTER_HH
