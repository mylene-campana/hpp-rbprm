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

#include <hpp/rbprm/interpolation/limb-rrt-helper.hh>
#include <hpp/rbprm/interpolation/limb-rrt-shooter.hh>
#include <hpp/rbprm/interpolation/limb-rrt-path-validation.hh>
#include <hpp/rbprm/interpolation/limb-rrt-steering.hh>
#include <hpp/core/steering-method-straight.hh>
#include <hpp/core/problem-target/goal-configurations.hh>
#include <hpp/core/bi-rrt-planner.hh>
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

    void LockRootAndNonContributingJoints(model::DevicePtr_t device, core::ConfigProjectorPtr_t& projector,
                                          const std::vector<std::string>& fixedContacts,
                                          const State& from, const State& to)
    {
        std::vector<std::string> spared = fixedContacts;
        to.contactCreations(from, spared);
        to.contactBreaks(from, spared);
        for(std::vector<std::string>::const_iterator cit = spared.begin(); cit != spared.end(); ++cit)
        {
            std::cout << "spared " << *cit << std::endl;
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

    
    void AddContactConstraints(LimbRRTHelper& helper, const State& from, const State& to)
    {
        std::vector<bool> cosntraintsR = setMaintainRotationConstraints();
        const rbprm::T_Limb& limbs = helper.fullbody_->GetLimbs();        
        std::vector<std::string> fixed = to.allFixedContacts(from,extractEffectorsName(limbs));
        core::Problem& problem = helper.rootProblem_;
        model::DevicePtr_t device = problem.robot();
        core::ConstraintSetPtr_t cSet = core::ConstraintSet::create(device,"");
        core::ConfigProjectorPtr_t proj = core::ConfigProjector::create(device,"proj", 1e-2, 30);
        for(std::vector<std::string>::const_iterator cit = fixed.begin();
            cit != fixed.end(); ++cit)
        {
            std::cout << "constraint " << *cit << std::endl;
            RbPrmLimbPtr_t limb = helper.fullbody_->GetLimbs().at(*cit);
            const fcl::Vec3f& ppos  = from.contactPositions_.at(*cit);

            JointPtr_t effectorJoint = device->getJointByName(limb->effector_->name());
            proj->add(core::NumericalConstraint::create (
                                    constraints::deprecated::Position::create("",device,
                                                                  effectorJoint,fcl::Vec3f(0,0,0), ppos)));
            if(limb->contactType_ == hpp::rbprm::_6_DOF && !to.ignore6DOF)
            {
                const fcl::Matrix3f& rotation = from.contactRotation_.at(*cit);

                proj->add(core::NumericalConstraint::create (constraints::deprecated::Orientation::create("", device,
                                                                                  effectorJoint,
                                                                                  rotation,
                                                                                  cosntraintsR)));

            }
        }

        cSet->addConstraint(proj);
        problem.constraints(cSet);

    }

    void SetPathValidation(LimbRRTHelper& helper)
    {
        LimbRRTPathValidationPtr_t pathVal = LimbRRTPathValidation::create(
                    helper.fullBodyDevice_, 0.05,helper.fullBodyDevice_->configSize()-1);
        helper.rootProblem_.pathValidation(pathVal);
    }


    }

    PathVectorPtr_t interpolateStates(LimbRRTHelper& helper, const State& from, const State& to)
    {
        PathVectorPtr_t res;
        core::PathPtr_t rootPath = helper.rootPath_;
        const rbprm::T_Limb& limbs = helper.fullbody_->GetLimbs();
        // get limbs that moved
        std::vector<std::string> variations = to.allVariations(from,extractEffectorsName(limbs));
        
        //std::vector<std::string> variations = extractEffectorsName(limbs);
        for(std::vector<std::string>::const_iterator cit = variations.begin();
            cit != variations.end(); ++cit)
        {
            SetPathValidation(helper);
            //DisableUnNecessaryCollisions(helper.rootProblem_, limbs.at(*cit));
            SetConfigShooter(helper,limbs.at(*cit),rootPath);

            ConfigurationPtr_t start = limbRRTConfigFromDevice(helper, from, 0.);
            ConfigurationPtr_t end   = limbRRTConfigFromDevice(helper, to  ,rootPath->length());
            helper.rootProblem_.initConfig(start);
            BiRRTPlannerPtr_t planner = BiRRTPlanner::create(helper.rootProblem_);
            ProblemTargetPtr_t target = problemTarget::GoalConfigurations::create (planner);
            helper.rootProblem_.target (target);
            helper.rootProblem_.addGoalConfig(end);
            AddContactConstraints(helper, from, to);

            res = planner->solve();
            helper.rootProblem_.resetGoalConfigs();
        }
        return res;
    }

    namespace
    {
        PathVectorPtr_t optimize(LimbRRTHelper& helper, PathVectorPtr_t partialPath, const std::size_t numOptimizations)
        {
            core::RandomShortcutPtr_t rs = core::RandomShortcut::create(helper.rootProblem_);
            for(std::size_t j=0; j<numOptimizations;++j)
            {
                partialPath = rs->optimize(partialPath);
            }
            return partialPath;
        }

        std::size_t checkPath(const std::size_t& distance, bool valid[])
        {
            std::size_t numValid(distance);
            for(std::size_t i = 0; i < distance; ++i)
            {
               if (!valid[i])
               {
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

        PathPtr_t ConcatenateAndResizePath(PathVectorPtr_t res[], std::size_t numValid)
        {
            PathVectorPtr_t completePath = res[0];
            for(std::size_t i = 1; i < numValid; ++i)
            {
                completePath->concatenate(*res[i]);
            }
            // reducing path
            core::SizeInterval_t interval(0, completePath->initial().rows()-1);
            core::SizeIntervals_t intervals;
            intervals.push_back(interval);
            PathPtr_t reducedPath = core::SubchainPath::create(completePath,intervals);
            return reducedPath;
        }
    }

    PathPtr_t interpolateStates(RbPrmFullBodyPtr_t fullbody, core::ProblemPtr_t referenceProblem, const PathPtr_t rootPath,
                                      const CIT_StateFrame &startState, const CIT_StateFrame &endState, const  std::size_t numOptimizations)
    {
        PathVectorPtr_t res[100];
        bool valid[100];
        std::size_t distance = std::distance(startState,endState);
        hppDout(notice,"InterpolateState distance = "<<distance);
        assert(distance < 100);
        // treat each interpolation between two states separatly
        // in a different thread
        #pragma omp parallel for
        for(std::size_t i = 0; i < distance; ++i)
        {
            CIT_StateFrame a, b;
            a = (startState+i);
            b = (startState+i+1);
            PathPtr_t extractedPath = rootPath->extract(core::interval_t(a->first, b->first));
            hppDout(notice," extracted path length = "<<extractedPath->length());
            LimbRRTHelper helper(fullbody, referenceProblem,
                                 rootPath->extract(core::interval_t(a->first, b->first)));
            PathVectorPtr_t partialPath = interpolateStates(helper, a->second, b->second);
            if(partialPath)
            {
                res[i] = optimize(helper,partialPath, numOptimizations);
                //res[i] =partialPath;
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

    PathPtr_t interpolateStates(RbPrmFullBodyPtr_t fullbody, core::ProblemPtr_t referenceProblem,
                                      const CIT_State &startState, const CIT_State &endState, const std::size_t numOptimizations)
    {
        PathVectorPtr_t res[100];
        bool valid[100];
        std::size_t distance = std::distance(startState,endState);
        assert(distance < 100);
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
  }// namespace interpolation
  }// namespace rbprm
}// namespace hpp
