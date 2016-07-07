//
// Copyright (c) 2015 CNRS
// Authors: Mylene Campana
//
// This file is part of hpp-rbprm
// hpp-rbprm is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-rbprm is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-rbprm  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_RBPRM_PRM_PLANNER_HH
# define HPP_RBPRM_PRM_PLANNER_HH

# include <boost/tuple/tuple.hpp>
# include <hpp/core/steering-method.hh>
# include <hpp/core/path-planner.hh>
# include <hpp/rbprm/planner/steering-method-parabola.hh>
# include <hpp/rbprm/planner/rbprm-roadmap.hh>

namespace hpp {
  namespace rbprm {
    /// \addtogroup path_planning
    /// \{

    // forward declaration of class Planner
    HPP_PREDEF_CLASS (PrmPlanner);
    // Planner objects are manipulated only via shared pointers
    typedef boost::shared_ptr <PrmPlanner> PrmPlannerPtr_t;

    /// Generic implementation of PRM algorithm
    class HPP_CORE_DLLAPI PrmPlanner : public core::PathPlanner
    {
      typedef boost::tuple <core::NodePtr_t, core::NodePtr_t,
			    core::PathPtr_t> DelayedEdge_t;
      typedef std::vector <DelayedEdge_t> DelayedEdges_t;
    public:
      /// Return shared pointer to new object.
      static PrmPlannerPtr_t createWithRoadmap
	(const core::Problem& problem, const core::RoadmapPtr_t& roadmap)
      {
	PrmPlanner* ptr = new PrmPlanner (problem, roadmap);
	return PrmPlannerPtr_t (ptr);
      }

      /// Return shared pointer to new object.
      static PrmPlannerPtr_t create (const core::Problem& problem)
      {
	PrmPlanner* ptr = new PrmPlanner (problem);
	return PrmPlannerPtr_t (ptr);
      }

      /// One step of extension.
      virtual void oneStep ();
      
      // disabled during testing
      virtual void tryDirectPath();

      /// Set configuration shooter.
      void configurationShooter (const core::ConfigurationShooterPtr_t& shooter)
      {
	configurationShooter_ = shooter;
      }

      virtual const core::RoadmapPtr_t& roadmap () const{
        return roadmap_;
      }

      const core::RbprmRoadmapPtr_t& rbprmRoadmap() const {
        return rbRoadmap_;
      }

    protected:
      /// Constructor
      PrmPlanner (const core::Problem& problem,
		  const core::RoadmapPtr_t& roadmap);

      /// Constructor with roadmap
      PrmPlanner (const core::Problem& problem);

      /// Store weak pointer to itself
      void init (const PrmPlannerWkPtr_t& weak)
      {
	core::PathPlanner::init (weak);
	weakPtr_ = weak;
      }

    private:
      core::ConfigurationShooterPtr_t configurationShooter_;
      mutable core::Configuration_t qProj_;
      PrmPlannerWkPtr_t weakPtr_;
      SteeringMethodParabolaPtr_t smParabola_;
      const core::RbprmRoadmapPtr_t rbRoadmap_;
      const core::RoadmapPtr_t roadmap_;
    };
    /// \}
  } // namespace rbprm
} // namespace hpp
#endif // HPP_RBPRM_PRM_PLANNER_HH
