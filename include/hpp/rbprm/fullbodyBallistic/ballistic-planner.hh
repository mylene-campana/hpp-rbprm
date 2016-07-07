//
// Copyright (c) 2016 CNRS
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

#ifndef HPP_RBPRM_BALLISTIC_PLANNER_HH
# define HPP_RBPRM_BALLISTIC_PLANNER_HH

# include <boost/tuple/tuple.hpp>
# include <hpp/core/steering-method.hh>
# include <hpp/core/path-planner.hh>
# include <hpp/rbprm/rbprm-fullbody.hh>
# include <hpp/rbprm/planner/steering-method-parabola.hh>
# include <hpp/rbprm/planner/rbprm-roadmap.hh>

namespace hpp {
  namespace rbprm {
    /// \addtogroup path_planning
    /// \{

    // forward declaration of class Planner
    HPP_PREDEF_CLASS (BallisticPlanner);
    // Planner objects are manipulated only via shared pointers
    typedef boost::shared_ptr <BallisticPlanner> BallisticPlannerPtr_t;

    /// Generic implementation of PRM algorithm
    class HPP_CORE_DLLAPI BallisticPlanner : public core::PathPlanner
    {
      typedef boost::tuple <core::NodePtr_t, core::NodePtr_t,
			    core::PathPtr_t> DelayedEdge_t;
      typedef std::vector <DelayedEdge_t> DelayedEdges_t;
    public:
      /// Return shared pointer to new object.
      static BallisticPlannerPtr_t createWithRoadmap
	(const core::Problem& problem, const core::RoadmapPtr_t& roadmap)
      {
	BallisticPlanner* ptr = new BallisticPlanner (problem, roadmap);
	return BallisticPlannerPtr_t (ptr);
      }

      /// Return shared pointer to new object.
      static BallisticPlannerPtr_t create (const core::Problem& problem)
      {
	BallisticPlanner* ptr = new BallisticPlanner (problem);
	return BallisticPlannerPtr_t (ptr);
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

      virtual const core::RbprmRoadmapPtr_t& rbprmRoadmap() const {
        return rbRoadmap_;
      }

    protected:
      /// Constructor with roadmap
      BallisticPlanner (const core::Problem& problem,
			const core::RoadmapPtr_t& roadmap);

      /// Constructor
      BallisticPlanner (const core::Problem& problem);

      /// Store weak pointer to itself
      void init (const BallisticPlannerWkPtr_t& weak)
      {
	core::PathPlanner::init (weak);
	weakPtr_ = weak;
      }

    private:

      /**
       * @brief computeGIWC compute the GIWC for the configuration and fill the node attribut, get validation report and call the second method
       * @param x the node
       */
      void computeGIWC (const core::Configuration_t q);

      /// return "average" direction of CoM cone, as average of contact cone 
      /// directions, if these cones are in the middle of the contact areas
      fcl::Vec3f computeMiddleContacts (const core::Configuration_t q) const;

      core::ProblemPtr_t problem_;
      core::ConfigurationShooterPtr_t configurationShooter_;
      BallisticPlannerWkPtr_t weakPtr_;
      SteeringMethodParabolaPtr_t smParabola_;
      core::RbprmRoadmapPtr_t rbRoadmap_;
      const core::RoadmapPtr_t roadmap_;
      const RbPrmFullBodyPtr_t fullRobot_; // for contact generation
      core::vector_t contactSize_; // should depend on the ROM
    };
    /// \}
  } // namespace rbprm
} // namespace hpp
#endif // HPP_RBPRM_BALLISTIC_PLANNER_HH
