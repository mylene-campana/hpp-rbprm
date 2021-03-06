//
// Copyright (c) 2016 CNRS
// Authors: Mylene Campana, Steve Tonneau (mcampana@laas.fr)
//
// This file is part of hpp-rbprm.
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
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_BALLISTIC_INTERPOLATION_HH
# define HPP_BALLISTIC_INTERPOLATION_HH

# include <hpp/rbprm/config.hh>
# include <hpp/rbprm/rbprm-fullbody.hh>
# include <hpp/core/path-vector.hh>
# include <hpp/core/problem.hh>
# include <hpp/model/device.hh>
#include <hpp/rbprm/planner/parabola-path.hh>

# include <vector>

namespace hpp {
  namespace rbprm {
    typedef std::pair<core::PathVectorPtr_t, BallisticPathPtr_t> PathVectorBP;
    typedef std::vector<PathVectorBP> T_PathVectorBP;

    core::Configuration_t computeContactPose(const State &state, core::Configuration_t contactPose,rbprm::RbPrmFullBodyPtr_t robot);

    HPP_PREDEF_CLASS(BallisticInterpolation);

    /// Interpolation class for transforming a path computed by RB-PRM into
    /// a discrete sequence of balanced contact configurations.
    ///
    class BallisticInterpolation;
    typedef boost::shared_ptr <BallisticInterpolation>
    BallisticInterpolationPtr_t;

    class HPP_RBPRM_DLLAPI BallisticInterpolation
    {
    public:
      /// Creates a smart pointer to the Interpolation class
      ///
      /// \param path the path returned by RB-PRM computation
      /// \param robot the FullBody instance considered for extending the part
      /// \param start the start full body configuration of the problem
      /// \param end the end full body configuration of the problem
      /// \return a pointer to the created BallisticInterpolation instance
      static BallisticInterpolationPtr_t create
	(const core::ProblemPtr_t& problem, const RbPrmFullBodyPtr_t robot,
	 const State& start, const State& end,
	 const core::PathVectorConstPtr_t path = core::PathVectorConstPtr_t()) {
	BallisticInterpolation* rbprmDevice =
	  new BallisticInterpolation (problem, robot, start, end, path);
	BallisticInterpolationPtr_t res (rbprmDevice);
	res->init (res);
	return res;
      }

      ~BallisticInterpolation(){}

      /// Transforms the path computed by RB-PRM into
      /// a discrete sequence of configurations.
      /// Do not return states anymore since no contact or stability
      /// TODO: collision avoidance ??
      ///
      /// \param timeStep the discretization step of the path.
      /// \return a pointer to the created BallisticInterpolation instance
      /*std::vector<model::Configuration_t> InterpolateConfigs
	(const double timeStep = 0.01);*/

      /// On each subpath, tranform the trunk path and the start-goal 
      /// states so that:
      /// the trunk is still following the parabola path
      /// the limbs are linearly interpolated
      /// \param u_offset: normalized curvilinear abcissa of contact maintain 
      /// for an cushion takeoff/landing. If u_offset = 0, this effect is
      /// disable.
      T_PathVectorBP InterpolateFullPath (const core::value_type u_offset = 0,
					  T_StateFrame* stateFramesRef = NULL);

      /// Between Start and End states, tranform the trunk path and the 
      /// start-goal states so that:
      /// the trunk is still following the parabola path
      /// the limbs are linearly interpolated
      /// \param u_offset: normalized curvilinear abcissa of contact maintain
      /// for an cushion takeoff/landing. If u_offset = 0, this effect is
      /// disable.
      T_PathVectorBP InterpolateDirectPath
	(const core::value_type u_offset = 0,
	 T_StateFrame* stateFramesRef = NULL);

      void extendingPose (const core::Configuration_t extendingPose) {
	extendingPose_ = extendingPose;
      }

      core::Configuration_t extendingPose () {
	return extendingPose_;
      }

      void flexionPose (const core::Configuration_t flexionPose) {
	flexionPose_ = flexionPose;
      }

      core::Configuration_t flexionPose () {
	return flexionPose_;
      }
      
      void takeoffContactPose (const core::Configuration_t takeoffContactPose) {
	takeoffContactPose_ = takeoffContactPose;
      }

      void landingContactPose (const core::Configuration_t landingContactPose) {
	landingContactPose_ = landingContactPose;
      }

      core::Configuration_t takeoffContactPose () {
	return takeoffContactPose_;
      }

      core::Configuration_t landingContactPose () {
	return landingContactPose_;
      }

      void flexionFinalPose (const core::Configuration_t flexionFinalPose) {
	flexionFinalPose_ = flexionFinalPose;
      }

      core::Configuration_t flexionFinalPose () {
	return flexionFinalPose_;
      }

      std::map<std::string, std::vector<std::string> > affordanceFilters () {
	return affFilters_;
      }
      
      affMap_t affordanceMap () { // often called "affordances"
	return affMap_;
      }

      void affordanceFilters (const std::map<std::string, std::vector<std::string> > affFilters) {
	affFilters_ = affFilters;
      }
      
      void affordanceMap (const affMap_t affMap) { // often called "affordances"
	affMap_ = affMap;
      }

      void setMaxIterMaintainContacts (const size_t maxIter) {
	maxIterMaintainContacts_ = maxIter;
      }

    protected:
      BallisticInterpolation (const core::ProblemPtr_t &problem,
			      const RbPrmFullBodyPtr_t robot,
			      const State& start, const State& end,
			      const core::PathVectorConstPtr_t path);

      ///
      /// \brief Initialization.
      ///
      void init (const BallisticInterpolationWkPtr_t& weakPtr) {
      weakPtr_ = weakPtr;
    }

      /// Compute direction vector for EFORT relatively to a config
      /// It is basically the difference between the landing velocity to the 
      /// config and the takeoff from the config.
      fcl::Vec3f computeDir (const core::vector_t V0,
			     const core::vector_t Vimp);

      /// Fill current trunk-config with 0 for limbs, and copy ECS part
      core::Configuration_t fillConfiguration
	(const core::Configuration_t config, const std::size_t configSize);

      /// Copy trunk and extra-config parts of trunkConfig, and
      /// limbs part (the rest) of refConfig
      core::Configuration_t fillConfiguration
	(const core::Configuration_t trunkConfig,
	 const core::Configuration_t refConfig);

      /// Compute the number of joints in the given limb (starting from limb
      /// and finishing to effector
      std::size_t computeLimbLength (const model::JointPtr_t limb,
				     const model::JointPtr_t effector);

      /// Replace limb-in-contact parts of configuration from q_contact
      /// to q_full
      core::Configuration_t replaceLimbConfigsInFullConfig
	(const core::Configuration_t q_full,
	 const core::Configuration_t q_contact,
	 const std::vector <RbPrmLimbPtr_t> limbs);

      /// Get a list of limbs that are active, i.e. which ROMs are in collision
      /// Comparison is made between the ROM name and the limb effector-name
      std::vector <RbPrmLimbPtr_t> activeLimbsFromROMs
	(const ParabolaPathPtr_t pp, const bool startConfig);

      std::vector<bool> setMaintainRotationConstraints(const fcl::Vec3f&);
      std::vector<bool> setTranslationConstraints(const fcl::Vec3f&);

      /// Add to the projector a LockJoint component
      void LockJointRec (const std::string& limb, const model::JointPtr_t joint,
			 core::ConfigProjectorPtr_t& projector);

      /*bool ComputeCollisionFreeConfiguration
	(State& current, core::CollisionValidationPtr_t validation,
	 const hpp::rbprm::RbPrmLimbPtr_t& limb,
	 model::ConfigurationOut_t configuration,
     const double robustnessTreshold, bool stability = true);*/

      /// Compute a state that *if possible* keeps the same contacts as
      /// in the previous state given.
  /*    State MaintainPreviousContacts
	(const State& previous,
	 std::map<std::string,core::CollisionValidationPtr_t>& limbValidations,
	 model::ConfigurationIn_t configuration, bool& contactMaintained,
     bool& multipleBreaks, std::vector <RbPrmLimbPtr_t>& successLimbs, const double robustnessTreshold = 0);*/

      /// Return the configuration based at u_offset on subpath, trying to 
      /// apply same contacts as in previousState. u is decreased 
      /// (or increased) of u_offset to get closer to previousState 
      /// until maxIterMaintainContacts is reached or 
      /// when contact cannot be maintained anymore.
      State computeOffsetContactConfig
	(const BallisticPathPtr_t bp, const State &previousState,
	 T_StateFrame& stateFrames, const core::value_type u_offset,
	 const bool increase_u_offset, core::value_type &length);

      /// Return the configuration at the top of the parabola (path),
      /// using extendingPose_ for limbs part if defined,
      /// otherwise, just unsing interpolation (bp)
      State computeTopExtendingPose(const core::PathPtr_t path,
				    const BallisticPathPtr_t bp,
				    core::value_type &length);

      
      /**
       * @brief computeTakeoffContactPose compute the configuration of a waypoint. In contact with the environnement.
       * Must be called after ComputeContact for the given State.
       * It return the configuration of the state exept for the trunk DOF and the limb which aren't in contact
       * @param state 
       * @return the configuration
       */
      core::Configuration_t computeTakeoffContactPose(const State& state) {
	if(takeoffContactPose_.size() == 0)
	  return state.configuration_;
	else
	  return rbprm::computeContactPose(state,takeoffContactPose_,robot_);
      }

/**
       * @brief computeLandingContactPose compute the configuration of a waypoint. In contact with the environnement.
       * Must be called after ComputeContact for the given State.
       * It return the configuration of the state exept for the trunk DOF and the limb which aren't in contact
       * @param state 
       * @return the configuration
       */
      core::Configuration_t computeLandingContactPose (const State& state) {
	if(landingContactPose_.size() == 0)
	  return state.configuration_;
	else
	  return rbprm::computeContactPose (state,landingContactPose_,robot_);
      }
      /**
       * @brief computeFlexionContactPose compute the configuration of a waypoint. In contact with the environnement.
       * Must be called after ComputeContact for the given State.
       * It return the configuration of the state exept for the trunk DOF and the limb which aren't in contact
       * @param state 
       * @return the configuration
       */
      core::Configuration_t computeFlexionContactPose (const State& state);
      
      /// Blend the two configurations with a as ratio:
      /// result = r*q1 + (1-r)*q2
      core::Configuration_t blendPoses (const core::Configuration_t q1,
					const core::Configuration_t q2,
					const core::value_type r);

      /// Find which limb is detected in collision and return it
      rbprm::RbPrmLimbPtr_t findLimbInCollision
	(const core::ValidationReportPtr_t validationReport);

      /// NO ! Replace the limbs not used for contact with their configuration in flexionPose_
      /// Replace the trunkDOF with contactPose value : (we suppose we always work with freeflyer as root)
      /// set EC of lastState from trunk config
      State replaceAccurateTrunkAndFakeLimbs
	(const State& refState,  const core::Configuration_t refTrunk,
	 const bool increase_u_offset,
	 const std::vector<std::string> contactingLimbs,
	 const core::Configuration_t lastfixedLimb);

      /// Create a state which contacts corresponds to the given contact-cones
      State createStateFromCone (const State previousState,
				 library::ContactCones cone);

      /// Replace all non-contacting limb configs with given reference config
      core::Configuration_t replaceNonContactingLimbConfig(State state, const core::Configuration_t refLimbConfig);
      
    private:
      const core::PathVectorConstPtr_t path_;
      State start_;
      State end_;
      core::ProblemPtr_t problem_;
      RbPrmFullBodyPtr_t robot_; // device of fullbody
      BallisticInterpolationWkPtr_t weakPtr_;
      // order of poses in parabola: flexion -> takeoffContact -> extending -> landingContact -> flexion
      core::Configuration_t extendingPose_; // parabola apex
      core::Configuration_t flexionPose_; // parabola extremity
      core::Configuration_t takeoffContactPose_; // when releasing contact
      core::Configuration_t landingContactPose_; // when starting contact
      core::Configuration_t flexionFinalPose_; // last flexion pose
      size_t lastRootIndex_;
      affMap_t affMap_;
      std::map<std::string, std::vector<std::string> > affFilters_;
      size_t trunkConfigSize_;
      size_t maxIterMaintainContacts_;
      int numberOfAdds_; // number of additionnal states between top and contacts, for manual interpolation of the COM
    }; // class BallisticInterpolation
    
  } // namespace rbprm
} // namespace hpp

#endif // HPP_BALLISTIC_INTERPOLATION_HH
