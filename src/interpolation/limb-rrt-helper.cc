// Copyright (c) 2014, LAAS-CNRS
// Authors: Steve Tonneau (steve.tonneau@laas.fr)
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
// hpp-rbprm. If not, see <http://www.gnu.org/licenses/>.

#include <hpp/util/debug.hh>
#include <hpp/model/configuration.hh>
#include <hpp/rbprm/interpolation/limb-rrt-helper.hh>
#include <hpp/rbprm/interpolation/limb-rrt-shooter.hh>
#include <hpp/rbprm/interpolation/limb-rrt-path-validation.hh>
#include <hpp/rbprm/interpolation/limb-rrt-steering.hh>
#include <hpp/core/steering-method-straight.hh>
#include <hpp/core/problem-target/goal-configurations.hh>
#include <hpp/core/bi-rrt-planner.hh>
#include <hpp/core/config-validations.hh>
#include <hpp/core/random-shortcut.hh>
#include <hpp/core/constraint-set.hh>>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/constraints/position.hh>
#include <hpp/constraints/orientation.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/core/locked-joint.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/core/subchain-path.hh>
#include <hpp/model/joint.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-path.hh>
#include <hpp/rbprm/tools.hh>

namespace hpp {
  using namespace core;
  using namespace model;
  namespace rbprm {
    namespace interpolation {
      using model::displayConfig;

      namespace{
	core::DevicePtr_t DeviceFromLimb(const std::string& name, RbPrmLimbPtr_t limb)
	{
	  DevicePtr_t limbDevice = Device::create(name);
	  limbDevice->rootJoint(limb->limb_->clone());
	  JointPtr_t current = limb->limb_, clone = limbDevice->rootJoint();
	  while(current->name() != limb->effector_->name())
	    {
	      current = current->childJoint(0);
	      clone->addChildJoint(current->clone());
	      clone = clone->childJoint(0);
	    }
	  limbDevice->setDimensionExtraConfigSpace(1);
	  return limbDevice;
	}
      }

      //find contact creation

      LimbRRTHelper::LimbRRTHelper(RbPrmFullBodyPtr_t fullbody, hpp::core::ProblemPtr_t referenceProblem, hpp::core::PathPtr_t rootPath)
        : fullbody_(fullbody)
        , fullBodyDevice_(fullbody->device_->clone())
        , rootProblem_(fullBodyDevice_)
        , rootPath_(rootPath)
      {
      
        // adding extra DOF for including time in sampling

        fullBodyDevice_->setDimensionExtraConfigSpace(fullBodyDevice_->extraConfigSpace().dimension()+1);
        /*for(ObjectVector_t::const_iterator cit = referenceProblem->collisionObstacles().begin() ; cit != referenceProblem->collisionObstacles().end() ; ++cit){
          rootProblem_.addObstacle(*cit);
	  }*/
        rootProblem_.collisionObstacles(referenceProblem->collisionObstacles());
        hppDout(notice,"REFERENCE PROBLEM OBSTACLE :"<<rootProblem_.collisionObstacles().size()<<" ; "<<referenceProblem->collisionObstacles().size());
	const core::value_type error_treshold = 0.001;
	proj_ = core::ConfigProjector::create(rootProblem_.robot(),"proj", error_treshold, 1000);
        BallisticPathPtr_t bp = boost::dynamic_pointer_cast<BallisticPath>(rootPath);
        hppDout(notice,"test last root index helper = "<<bp->lastRootIndex());
        if(bp){
          rootProblem_.steeringMethod(LimbRRTSteering::create(&rootProblem_,fullBodyDevice_->configSize()-1,bp));
          hppDout(notice,"create steering method with ballistic path.");
        }
        else{
          rootProblem_.steeringMethod(LimbRRTSteering::create(&rootProblem_,fullBodyDevice_->configSize()-1));
          hppDout(notice,"create steering method without ballistic path.");
        }
	rootProblem_.plannerIterLimit (referenceProblem->plannerIterLimit ());
	hppDout (info, "initial plannerIterLimit_ (from referenceProblem)= " << rootProblem_.plannerIterLimit ());
          
      }

      namespace
      {
	core::PathPtr_t generateRootPath(const Problem& problem, const State& from, const State& to)
	{
	  Configuration_t startRootConf(from.configuration_);
	  Configuration_t endRootConf(to.configuration_);
	  return (*(problem.steeringMethod()))(startRootConf, endRootConf);
	}

	void DisableUnNecessaryCollisions(core::Problem& problem, rbprm::RbPrmLimbPtr_t limb)
	{
	  // TODO should we really disable collisions for other bodies ?
	  hppDout(notice,"REFERENCE PROBLEM OBSTACLE :"<<problem.collisionObstacles().size());
	  tools::RemoveNonLimbCollisionRec<core::Problem>(problem.robot()->rootJoint(),                                                        limb->limb_->name(), problem.collisionObstacles(),problem);

	  if(limb->disableEndEffectorCollision_)
	    {
	      hpp::tools::RemoveEffectorCollision<core::Problem>(problem,
								 problem.robot()->getJointByName(limb->effector_->name()),
								 problem.collisionObstacles());
	    }
	}

	ConfigurationPtr_t limbRRTConfigFromDevice(const LimbRRTHelper& helper, const State& state, const double time)
	{
	  Configuration_t config(helper.fullBodyDevice_->currentConfiguration());
	  config.head(state.configuration_.rows()) = state.configuration_;
	  config[config.rows()-1] = time;
	  return ConfigurationPtr_t(new Configuration_t(config));
	}

	void SetConfigShooter(LimbRRTHelper& helper, RbPrmLimbPtr_t limb, core::PathPtr_t& rootPath)
	{
	  ConfigurationShooterPtr_t limbRRTShooter = LimbRRTShooter::create(limb, rootPath,
									    helper.fullBodyDevice_->configSize()-1);
	  helper.rootProblem_.configurationShooter(limbRRTShooter);
	}

	std::vector<bool> setMaintainRotationConstraints() // direction)
	{
	  std::vector<bool> res;
	  for(std::size_t i =0; i <3; ++i)
	    {
	      res.push_back(true);
	    }
	  return res;
	}

	void LockRootAndNonContributingJoints
	(model::DevicePtr_t device, core::ConfigProjectorPtr_t& projector,
	 const std::vector<std::string>& fixedContacts,
	 const State& from, const State& to)
	{
	  std::vector<std::string> spared = fixedContacts;
	  to.contactCreations(from, spared);
	  to.contactBreaks(from, spared);
	  for(std::vector<std::string>::const_iterator cit = spared.begin(); cit != spared.end(); ++cit)
	    {
	      //std::cout << "spared " << *cit << std::endl;
	    }
	  tools::LockJointRec(spared, device->rootJoint(), projector);
	}

	std::vector<std::string> extractEffectorsName(const rbprm::T_Limb& limbs)
	{
	  std::vector<std::string> res;
	  for(rbprm::T_Limb::const_iterator cit = limbs.begin(); cit != limbs.end(); ++cit)
	    {
	      hppDout(info,"all effector names + "<<cit->first);
	      res.push_back(cit->first);
	    }
	  return res;
	}

    
	void AddContactConstraints(LimbRRTHelper& helper, const State& from,
				   const State& to)
	{
	  std::vector<bool> cosntraintsR = setMaintainRotationConstraints();
	  //const rbprm::T_Limb& limbs = helper.fullbody_->GetLimbs();
	  hppDout (info, "fix from contacts");
	  hppDout (info, "from= " << displayConfig(from.configuration_));
	  hppDout (info, "to= " << displayConfig(to.configuration_));
	  //std::vector<std::string> fixed = to.allFixedContacts(from,extractEffectorsName(limbs)); // problem: seems to return all limb in contact
	  std::vector<std::string> fixed = to.fixedContacts(from);
	  core::Problem& problem = helper.rootProblem_;
	  model::DevicePtr_t device = problem.robot();
	  core::ConstraintSetPtr_t cSet = core::ConstraintSet::create(device,"");
	  core::ConfigProjectorPtr_t proj = core::ConfigProjector::create(device,"proj", 1e-2, 30);
	  bool disableConstr = true;
	  bool atLeastOneConstr = false;

	  // For verification:
	  Configuration_t initialLong (from.configuration_.size () + 1);
	  for (std::size_t i = 0; i < from.configuration_.size (); i++)
	    initialLong [i] = from.configuration_ [i];
	  initialLong [from.configuration_.size ()] = 0;
	  Configuration_t endLong (to.configuration_.size () + 1);
	  for (std::size_t i = 0; i < to.configuration_.size (); i++)
	    endLong [i] = to.configuration_ [i];
	  endLong [to.configuration_.size ()] = 0;
	  //hppDout (info, "endLong config= " << displayConfig(endLong));

	  for(std::vector<std::string>::const_iterator cit = fixed.begin();
	      cit != fixed.end(); ++cit)
	    {
	      hppDout (info, "create projector");
	      proj = core::ConfigProjector::create(device,"proj", 1e-2, 30);
	      hppDout (info, "fixed contact: " << *cit);
	      RbPrmLimbPtr_t limb = helper.fullbody_->GetLimbs().at(*cit);
	      const bool isInContact = from.contacts_.at(*cit);
	      if (isInContact) {
		const fcl::Vec3f& ppos  = from.contactPositions_.at(*cit);
		const fcl::Vec3f& ppos_to  = to.contactPositions_.at(*cit);
	      
		const JointPtr_t effectorJoint = device->getJointByName(limb->effector_->name()); // limb->effector_ directly NOT WORKING for ConfigProjector
	      
		const fcl::Transform3f& transform =  effectorJoint->currentTransformation ();

		hppDout (info, "effector joint= " << limb->effector_->name());
		hppDout (info, "pos from= " << ppos);
		hppDout (info, "pos to= " << ppos_to);

		/// DEBUG test configs
		device->currentConfiguration(initialLong);
		device->computeForwardKinematics();
		fcl::Vec3f configPos = effectorJoint->currentTransformation ().getTranslation();
		const value_type normFrom = (ppos - configPos).norm ();
		hppDout (info, "real effector pos config from= " << configPos);
		hppDout (info, "normFrom= " << normFrom);

		device->currentConfiguration(endLong);
		device->computeForwardKinematics();
		configPos = effectorJoint->currentTransformation ().getTranslation();
		const value_type normTo = (ppos_to - configPos).norm ();
		hppDout (info, "real effector pos config to= " << configPos);
		hppDout (info, "normTo= " << normTo);

		if (normFrom < 1e-4 && normTo < 1e-4) {
		  // create constraint(s)
		  hppDout (info, "create contact position constraint");
		  proj->add(core::NumericalConstraint::create (constraints::Position::create("", device, effectorJoint, fcl::Vec3f(0,0,0), ppos)));
		  disableConstr = false;
		  if (!proj->isSatisfied (endLong)) {
		    disableConstr = true;
		    hppDout (error, "End configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION ONLY");
		    //throw projection_error ("End configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION ONLY");
		  }
		  if (!proj->isSatisfied (initialLong)) {
		    disableConstr = true;
		    hppDout (error, "Initial configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION ONLY");
		    //throw projection_error ("Initial configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION ONLY");
		  }
		  const bool ignoreRotation = to.ignoreRotation_.at(*cit);
		  hppDout (info, "is rot constraint bypassed ?= " <<ignoreRotation);
		  if(limb->contactType_ == hpp::rbprm::_6_DOF && !to.ignore6DOF 
		     && !ignoreRotation) {
		    hppDout (info, "create contact rotation constraint");
		    const fcl::Matrix3f& rotation = from.contactRotation_.at(*cit);
		    hppDout (info, "rot from= " << rotation);
		    hppDout (info, "rot to= " << from.contactRotation_.at(*cit));
		    proj->add(core::NumericalConstraint::create (constraints::Orientation::create("", device, effectorJoint, rotation, cosntraintsR)));
		    if (!proj->isSatisfied (endLong)) {
		      disableConstr = true;
		      hppDout (error, "End configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION + ROTATION");
		      //throw projection_error ("End configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION + ROTATION");
		    }
		    if (!proj->isSatisfied (initialLong)) {
		      disableConstr = true;
		      hppDout (error, "Initial configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION + ROTATION");
		      //throw projection_error ("Initial configuration of limbRRT does not satisfy the contact-constraints of the initial config-- POSITION + ROTATION");
		    }
		  }// if rotation
		}// if norm OK
		else
		  hppDout (info, "constraint not added to limbRRT because of norm");
	      }// if contact
	      else
		hppDout (info, "is actually NOT in contact in from !!");

	      if (!disableConstr) {
		hppDout (info, "constr not disabled");
		//cSet->addConstraint(proj);     // DEBUG !!!
		atLeastOneConstr = true;
		hppDout (info, "constr not added because DEBUG");
	      } else
		hppDout (info, "constr disabled");

	    }//for fixed limbs

	  if (atLeastOneConstr) {
	    hppDout (info, "add all new constraints to problem");
	    //problem.constraints(cSet);         // DEBUG !!!
	  } else {
	    hppDout (info, "NO new constraint to add to problem");
	  }

	  hppDout (info, "test end config with configProjector build with initial config");
	  if (!proj->isSatisfied (endLong)) {
	    //if (!problem.constraints()->isSatisfied (endLong)) {
	    hppDout (error, "End configuration of limbRRT does not satisfy the contact-constraints of the initial config");
	    throw projection_error ("End configuration of limbRRT does not satisfy the contact-constraints of the initial config");
	  }

	}

	void SetPathValidation(LimbRRTHelper& helper)
	{
	  LimbRRTPathValidationPtr_t pathVal =
	    LimbRRTPathValidation::create(helper.fullBodyDevice_, 0.05,helper.fullBodyDevice_->configSize()-1);
	  helper.rootProblem_.pathValidation(pathVal);
	}
      }

      void InitConstraints (LimbRRTHelper& helper)
      {
        core::ConstraintSetPtr_t cSet = core::ConstraintSet::create(helper.rootProblem_.robot(),"");
        cSet->addConstraint(helper.proj_);
        helper.rootProblem_.constraints(cSet);
      }

      PathVectorPtr_t interpolateStates(LimbRRTHelper& helper, const State& from, const State& to)
      {
        PathVectorPtr_t res;
        core::PathPtr_t rootPath = helper.rootPath_;
        const rbprm::T_Limb& limbs = helper.fullbody_->GetLimbs();
        // get limbs that moved
        std::vector<std::string> variations = to.allVariations(from,extractEffectorsName(limbs));
	hppDout (info, "number of variations: " << variations.size ());
        core::ValidationReportPtr_t validationReport;

	if (variations.size () != 0) {
	  //AddContactConstraints(helper, from, to);  // DEBUG !!
	
	  //std::vector<std::string> variations = extractEffectorsName(limbs);
	  for(std::vector<std::string>::const_iterator cit = variations.begin();
	      cit != variations.end(); ++cit)
	    {
	      hppDout (info, "variation: " << *cit);
	      SetPathValidation(helper);
	      //DisableUnNecessaryCollisions(helper.rootProblem_, limbs.at(*cit));
	      SetConfigShooter(helper,limbs.at(*cit),rootPath);

	      ConfigurationPtr_t start = limbRRTConfigFromDevice(helper, from, 0.);
	      ConfigurationPtr_t end   = limbRRTConfigFromDevice(helper, to  ,rootPath->length());
	      hppDout (info, "BiRRT planner start: " << displayConfig (from.configuration_));
	      hppDout (info, "BiRRT planner end: " << displayConfig (to.configuration_));
	      const model::DevicePtr_t robot = helper.fullbody_->device_;

	      if (!helper.rootProblem_.configValidations ()->validate (*start, validationReport))
		hppDout (info, "start config NOT valid");
	      if (!helper.rootProblem_.configValidations ()->validate (*end, validationReport))
		hppDout (info, "end config NOT valid");
	      helper.rootProblem_.initConfig(start);
	      hppDout (info, "create BiRRT planner");
	      BiRRTPlannerPtr_t planner = BiRRTPlanner::create(helper.rootProblem_);
	      ProblemTargetPtr_t target = problemTarget::GoalConfigurations::create (planner);
	      helper.rootProblem_.target (target);
	      helper.rootProblem_.addGoalConfig(end);

	      hppDout (info, "before BiRRT planner solve");
	      res = planner->solve();
	      hppDout (info, "after BiRRT planner solve");
	      helper.rootProblem_.resetGoalConfigs();
	      hppDout (info, "helper.rootProblem_.nbPathPlannerFails_= " << helper.rootProblem_.nbPathPlannerFails_);
	    }// for variations

	  // if no variation, still solve the problem to keep the contact along the path
	} else {
	  hppDout (info, "variation null, try straight path!");

	  res = directInterpolation (helper, from, to);
	  PathValidationPtr_t pathValidation = helper.rootProblem_.pathValidation ();
	  core::PathPtr_t validPart;
	  core::PathValidationReportPtr_t pathReport;
	    
	  if (!pathValidation->validate (res, false, validPart, pathReport)) {
	    hppDout (info, "straight path invalid, limbRRT for each limb without constraints");
	    for(T_Limb::const_iterator lit = limbs.begin(); lit != limbs.end(); ++lit)
	      {
		hppDout (info, "limb= " << lit->first);
		SetPathValidation(helper);
		//DisableUnNecessaryCollisions(helper.rootProblem_, limbs.at(*cit));
		SetConfigShooter(helper,lit->second,rootPath);
		ConfigurationPtr_t start = limbRRTConfigFromDevice(helper, from,0.);
		ConfigurationPtr_t end   = limbRRTConfigFromDevice(helper, to ,rootPath->length());
		hppDout (info, "BiRRT planner start: " << displayConfig (from.configuration_));
		hppDout (info, "BiRRT planner end: " << displayConfig (to.configuration_));
		const model::DevicePtr_t robot = helper.fullbody_->device_;
		if (!helper.rootProblem_.configValidations ()->validate (*start, validationReport))
		  hppDout (info, "start config NOT valid");
		if (!helper.rootProblem_.configValidations ()->validate (*end, validationReport))
		  hppDout (info, "end config NOT valid");
		helper.rootProblem_.initConfig(start);
		hppDout (info, "create BiRRT planner");
		BiRRTPlannerPtr_t planner = BiRRTPlanner::create(helper.rootProblem_);
		ProblemTargetPtr_t target = problemTarget::GoalConfigurations::create (planner);
		helper.rootProblem_.target (target);
		helper.rootProblem_.addGoalConfig(end);

		hppDout (info, "before BiRRT planner solve");
		res = planner->solve();
		hppDout (info, "after BiRRT planner solve");
		helper.rootProblem_.resetGoalConfigs();
		hppDout (info, "helper.rootProblem_.nbPathPlannerFails_= " << helper.rootProblem_.nbPathPlannerFails_);
	      }// for limbs
	  }// if straight-path invalid
	}// if no variation
	hppDout (info, "number of subpaths in res= " << res->numberPaths ());
        return res;
      }

      namespace
      {
        PathVectorPtr_t optimize(LimbRRTHelper& helper,
				 PathVectorPtr_t partialPath,
				 const std::size_t numOptimizations)
        {
	  core::RandomShortcutPtr_t rs = core::RandomShortcut::create(helper.rootProblem_);
	  for(std::size_t j=0; j<numOptimizations;++j)
            {
	      partialPath = rs->optimize(partialPath);
            }
	  return partialPath;
        }

        std::size_t checkPath(const std::size_t& distance, std::vector<bool> valid)
        {
	  std::size_t numValid(distance);
	  for(std::size_t i = 0; i < distance; ++i)
            {
	      if (!valid[i])
		{
		  hppDout (info, "not valid= " << i);
		  numValid= i;
		  break;
		}
            }
	  if (numValid==0)
	    throw std::runtime_error("No path found at state 0");
	  else if(numValid != distance)
            {
	      std::cout << "No path found at state " << numValid << std::endl;
            }
	  return numValid;
        }

	PathPtr_t ConcatenateAndResizePath (std::vector<PathVectorPtr_t> res,
					    std::size_t numValid)
        {
	  PathVectorPtr_t completePath = res[0];
	  for(std::size_t i = 1; i < numValid; ++i)
	    completePath->concatenate(*res[i]);
	  // reducing path
	  core::SizeInterval_t interval(0, completePath->initial().rows()-1);
	  core::SizeIntervals_t intervals;
	  intervals.push_back(interval);
	  PathPtr_t reducedPath = core::SubchainPath::create(completePath,
							     intervals);
	  return reducedPath;
	}

	PathVectorPtr_t ConcatenatePathInPathVector
	(std::vector<PathVectorPtr_t> res, std::size_t numValid)
	{
	  //Note: res pathVectors has configs of size +1
	  const std::size_t ReducedConfigSize = res[0]->outputSize () - 1;
	  hppDout (info, "ReducedConfigSize= " << ReducedConfigSize);
	  PathVectorPtr_t result = core::PathVector::create (ReducedConfigSize, res[0]->outputDerivativeSize ()-1);
	  core::SizeInterval_t interval (0, ReducedConfigSize);
	  core::SizeIntervals_t intervals;
	  intervals.push_back(interval);
	  
	  for(std::size_t i = 0; i < numValid; i++) {
	    PathPtr_t reducedPath_i = core::SubchainPath::create(res[i],
								 intervals);
	    result->appendPath (reducedPath_i);
	  }
	  hppDout (info, "result number of paths= " << result->numberPaths ());

	  /*PathVectorPtr_t completePath = res[0];
	  hppDout (info, "size of res[0]= " << res[0]->numberPaths ());
	  hppDout (info, "init = " <<displayConfig (res[0]->initial ()));
	  hppDout (info, "end = " <<displayConfig (res[0]->end ()));
	  for(std::size_t i = 1; i < numValid; ++i)
	    {
	      hppDout (info, "size of res ["<<i<<"]= " << res[i]->numberPaths ());
	      hppDout (info, "init = " <<displayConfig (res[i]->initial ()));
	      hppDout (info, "end = " <<displayConfig (res[i]->end ()));
	      completePath->concatenate(*res[i]); // bout a bout
	    }
	  hppDout (info, "size of completePath= " << completePath->numberPaths ());
	  return completePath;*/
	  return result;
	}
      }

      PathVectorPtr_t interpolateStatesinPathVector
      (RbPrmFullBodyPtr_t fullbody, core::ProblemPtr_t referenceProblem,
       const PathPtr_t rootPath, const CIT_StateFrame &startState,
       const CIT_StateFrame &endState, const std::size_t numOptimizations)
      {
	std::size_t distance = std::distance(startState,endState);
	std::vector<PathVectorPtr_t> res (distance + 1);
	std::vector<bool> valid (distance + 1);
	hppDout(notice,"InterpolateState rootpath distance = "<<distance);
	// treat each interpolation between two states separatly
	// in a different thread5
	//#pragma omp parallel for
	for(std::size_t i = 0; i < distance; ++i)
	  {
	    hppDout (info, "limbRRT interpolation from state= " << i);
	    CIT_StateFrame a, b;
	    a = (startState+i);
	    b = (startState+i+1);
	    hppDout (info, "extract path");
	    PathPtr_t extractedPath = rootPath->extract(core::interval_t(a->first, b->first));
	    hppDout(notice," extracted path length = "<<extractedPath->length());
	    hppDout (info, "create local helper");
	    LimbRRTHelper helper(fullbody, referenceProblem,
				 rootPath->extract(core::interval_t(a->first, b->first)));
	    hppDout (info, "create partial path from interpolateStates");
	    PathVectorPtr_t partialPath = interpolateStates(helper, a->second, b->second);
	    if(partialPath)
	      {
		if (helper.rootProblem_.nbPathPlannerFails_ > 0)
		  referenceProblem->nbPathPlannerFails_++;
		hppDout (info, "referenceProblem->nbPathPlannerFails_= " << referenceProblem->nbPathPlannerFails_);
		//hppDout (info, "optimize partial path");
		//res[i] = optimize(helper,partialPath, numOptimizations); // this will create subpaths in res[i]
		hppDout (info, "partial path from interpolateStates exists");
		hppDout (info, "numberPaths= " << partialPath->numberPaths ());
		if (!(partialPath->numberPaths ())) {
		  hppDout (info, "problem with partialPath from interpolateStates");
		  hppDout (info, "numberPaths= " << partialPath->numberPaths ());
		  // use direct interpolation instead ?
		  partialPath = directInterpolation (helper, a->second, b->second);
		}
		res[i] = partialPath;
		valid[i]=true;
	      }
	    else
	      {
		valid[i] = false;
		hppDout (info, "partial path from interpolateStates is empty");
	      }
	  }
	hppDout (info, "check path");
	std::size_t numValid = checkPath(distance, valid);
	hppDout (info, "ConcatenatePath");
	hppDout (info, "distance= " << distance << " numValid= " << numValid);
	hppDout (info, "referenceProblem->nbPathPlannerFails_= " << referenceProblem->nbPathPlannerFails_);
	// numValid is the index of the last existing path in pathVector
	return ConcatenatePathInPathVector (res, numValid);
      }

      PathPtr_t interpolateStates
      (RbPrmFullBodyPtr_t fullbody, core::ProblemPtr_t referenceProblem,
       const PathPtr_t rootPath, const CIT_StateFrame &startState,
       const CIT_StateFrame &endState, const std::size_t numOptimizations)
      {
	std::size_t distance = std::distance(startState,endState);
	std::vector<PathVectorPtr_t> res (distance + 1);
	std::vector<bool> valid (distance + 1);
	// treat each interpolation between two states separatly
	// in a different thread
	//#pragma omp parallel for
	for(std::size_t i = 0; i < distance; ++i)
	  {
	    CIT_StateFrame a, b;
	    a = (startState+i);
	    b = (startState+i+1);
	    PathPtr_t extractedPath = rootPath->extract(core::interval_t(a->first, b->first));
	    LimbRRTHelper helper(fullbody, referenceProblem,
				 rootPath->extract(core::interval_t(a->first, b->first)));
	    PathVectorPtr_t partialPath = interpolateStates(helper, a->second, b->second);
	    if(partialPath)
	      {
		res[i] = optimize(helper,partialPath, numOptimizations);
		res[i] =partialPath;
		valid[i]=true;
	      }
	    else
	      {
		valid[i] = false;
		hppDout (info, "partial path from interpolateStates is empty");
	      }
	  }
	std::size_t numValid = checkPath(distance, valid);
	return ConcatenateAndResizePath(res, numValid);
      }

      PathPtr_t interpolateStates(RbPrmFullBodyPtr_t fullbody,
				  core::ProblemPtr_t referenceProblem,
				  const CIT_State &startState,
				  const CIT_State &endState,
				  const std::size_t numOptimizations)
      {
	std::size_t distance = std::distance(startState,endState);
	std::vector<PathVectorPtr_t> res (distance + 1);
	std::vector<bool> valid (distance + 1);
	hppDout(notice,"InterpolateState distance = "<<distance);
	// treat each interpolation between two states separatly
	// in a different thread
#pragma omp parallel for
	for(std::size_t i = 0; i < distance; ++i)
	  {
	    CIT_State a, b;
	    a = (startState+i);
	    b = (startState+i+1);
	    LimbRRTHelper helper(fullbody, referenceProblem, generateRootPath(*referenceProblem, *a, *b));
	    PathVectorPtr_t partialPath = interpolateStates(helper, *a, *b);
	    if(partialPath)
	      {
		res[i] = optimize(helper,partialPath, numOptimizations);
		valid[i]=true;
	      }
	    else
	      {
		valid[i] = false;
	      }
	  }
	std::size_t numValid = checkPath(distance, valid);
	return ConcatenateAndResizePath(res, numValid);
      }

      PathVectorPtr_t directInterpolation (LimbRRTHelper& helper,
					   const State& from, const State& to)
      {
	ConfigurationPtr_t start = limbRRTConfigFromDevice(helper, from, 0.);
	ConfigurationPtr_t end = limbRRTConfigFromDevice(helper, to, helper.rootPath_->length());
	PathVectorPtr_t res = core::PathVector::create ((*start).size (), (*start).size () - 1);
	Problem& problem = helper.rootProblem_;
	const DistancePtr_t& distance = problem.distance();
	const value_type length = (*distance) (from.configuration_, to.configuration_);
	// build direct interpolation
	hppDout (info, "before straightPath creation");
	StraightPathPtr_t path = StraightPath::create (problem.robot (), *start, *end, length);
	res->appendPath (path);
	//assert(0);
	hppDout (info, "direct interpolation for partialPath");
	return res;
	}

    }// namespace interpolation
  }// namespace rbprm
}// namespace hpp
