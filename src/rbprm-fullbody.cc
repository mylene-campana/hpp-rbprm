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

#include <hpp/model/configuration.hh>
#include <hpp/rbprm/rbprm-fullbody.hh>
#include <hpp/model/joint.hh>
#include <hpp/rbprm/tools.hh>
#include <hpp/rbprm/stability/stability.hh>
#include <hpp/rbprm/ik-solver.hh>

#include <hpp/core/constraint-set.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/core/locked-joint.hh>
#include <hpp/model/device.hh>
#include <hpp/constraints/generic-transformation.hh>


#include <hpp/fcl/BVH/BVH_model.h>

#include <stack>

#ifdef PROFILE
    #include "hpp/rbprm/rbprm-profiler.hh"
#endif

namespace hpp {
  namespace rbprm {

    using model::displayConfig;
    const double epsilon = 1e-3;

    bool isCentroidalConeValid (const core::vector_t& n,
				const core::vector_t& V0Dir,
				const core::vector_t& VfDir,
				const core::value_type& mu) {
      const core::vector_t V0dir = V0Dir/V0Dir.norm ();
      const core::vector_t Vfdir = -VfDir/VfDir.norm ();
      const core::value_type a1 = fabs(acos(V0dir.dot(n)));
      const core::value_type a2 = fabs(acos(Vfdir.dot(n)));
      const core::value_type phi = atan(mu);
      hppDout (info, "a1= " << a1);
      hppDout (info, "a2= " << a2);
      hppDout (info, "phi= " << phi);
      if (a1 < phi && a2 < phi)
	return true;
      else
	return false;
    }

    fcl::Vec3f updateConeDirection (const State state) {
      const std::size_t contactNumber = state.contacts_.size ();
      fcl::Vec3f normalAv = (0,0,0);
      for(std::map<std::string, bool>::const_iterator cit = state.contacts_.begin() ; cit != state.contacts_.end() ; ++cit){
	if(cit->second){ // in contact
	  const fcl::Vec3f& normal = state.contactNormals_.at(cit->first);
	  for (std::size_t j = 0; j < 3; j++)
	    normalAv [j] += normal [j]/contactNumber;
	}
      }
      normalAv.normalize ();
      return normalAv;
    }

    // Initialize ExtraConfigs of current state config with:
    // - zeros if there is no contact or if it is not the cone direction (theta)
    // - normalAv (average of contact directions) if there are contacts
    State initializeExtraConfigs (const hpp::rbprm::RbPrmFullBodyPtr_t& body,
				  const State state) {
      State result = state;
      const fcl::Vec3f normalAv = updateConeDirection (state);
      hppDout (info, "normalAv= " << normalAv);
      const std::size_t ecSize = body->device_->extraConfigSpace ().dimension ();
      for (std::size_t k = 0; k < ecSize; k++) {
	if (k < 3)
	  result.configuration_ [body->device_->configSize () - ecSize + k] = normalAv [k];
	else
	  result.configuration_ [body->device_->configSize () - ecSize + k] = 0;
      }
      hppDout (info, "config after  updateConeDirection + initializeExtraConfigs= " << displayConfig (result.configuration_));
      return result;
    }

    RbPrmFullBodyPtr_t RbPrmFullBody::create (const model::DevicePtr_t &device)
    {
        RbPrmFullBody* fullBody = new RbPrmFullBody(device);
        RbPrmFullBodyPtr_t res (fullBody);
        res->init (res);
        return res;
    }

    RbPrmFullBody::~RbPrmFullBody()
    {
        // NOTHING
    }

    bool RbPrmFullBody::AddHeuristic(const std::string& name, const sampling::heuristic func)
    {
        return factory_.AddHeuristic(name, func);
    }


    void RbPrmFullBody::AddLimbPrivate(rbprm::RbPrmLimbPtr_t limb, const std::string& id, const std::string& name,
                        const model::ObjectVector_t &collisionObjects, const bool disableEffectorCollision)
    {
        core::CollisionValidationPtr_t limbcollisionValidation_ = core::CollisionValidation::create(this->device_);
        // adding collision validation
        for(model::ObjectVector_t::const_iterator cit = collisionObjects.begin();
            cit != collisionObjects.end(); ++cit)
        {
            if(limbs_.empty())
            {
                collisionValidation_->addObstacle(*cit);
            }
            limbcollisionValidation_->addObstacle(*cit);
	    hppDout (info, "limb collisionValidation:" << limb->limb_->name());
	    hppDout (info, "obst collisionValidation:" << (*cit)->name());
            //remove effector collision
            if(disableEffectorCollision)
            {
                hpp::tools::RemoveEffectorCollision<core::CollisionValidation>((*collisionValidation_.get()), limb->effector_, *cit);
                hpp::tools::RemoveEffectorCollision<core::CollisionValidation>((*limbcollisionValidation_.get()), limb->effector_, *cit);
            }
        }
        limbs_.insert(std::make_pair(id, limb));
        tools::RemoveNonLimbCollisionRec<core::CollisionValidation>(device_->rootJoint(),name,collisionObjects,*limbcollisionValidation_.get());
        hpp::core::RelativeMotion::matrix_type m = hpp::core::RelativeMotion::matrix(device_);
        limbcollisionValidation_->filterCollisionPairs(m);
        collisionValidation_->filterCollisionPairs(m);
        limbcollisionValidations_.insert(std::make_pair(id, limbcollisionValidation_));
        // insert limb to root group
        T_LimbGroup::iterator cit = limbGroups_.find(name);
        if(cit != limbGroups_.end())
        {
            cit->second.push_back(id);
        }
        else
        {
            std::vector<std::string> group;
            group.push_back(id);
            limbGroups_.insert(std::make_pair(name, group));
        }
    }

    std::map<std::string, const sampling::heuristic>::const_iterator checkLimbData(const std::string& id, const rbprm::T_Limb& limbs, const rbprm::sampling::HeuristicFactory& factory, const std::string& heuristicName)
    {
        rbprm::T_Limb::const_iterator cit = limbs.find(id);
        std::map<std::string, const sampling::heuristic>::const_iterator hit = factory.heuristics_.find(heuristicName);
        if(cit != limbs.end())
            throw std::runtime_error ("Impossible to add limb for joint "
                                      + id + " to robot; limb already exists");
        else if(hit == factory.heuristics_.end())
            throw std::runtime_error ("Impossible to add limb for joint "
                                      + id + " to robot; heuristic not found " + heuristicName +".");
        return hit;
    }

    void RbPrmFullBody::AddLimb(const std::string& id, const std::string& name, const std::string &effectorName,
                                const fcl::Vec3f &offset,const fcl::Vec3f &normal, const double x,
                                const double y,
                                const model::ObjectVector_t &collisionObjects, const std::size_t nbSamples, const std::string &heuristicName, const double resolution,
                                ContactType contactType, const bool disableEffectorCollision)
    {
        std::map<std::string, const sampling::heuristic>::const_iterator hit = checkLimbData(id, limbs_,factory_,heuristicName);
        model::JointPtr_t joint = device_->getJointByName(name);
        rbprm::RbPrmLimbPtr_t limb = rbprm::RbPrmLimb::create(joint, effectorName, offset,normal,x,y, nbSamples, hit->second, resolution,contactType, disableEffectorCollision);
        AddLimbPrivate(limb, id, name,collisionObjects, disableEffectorCollision);
    }

    void RbPrmFullBody::AddLimb(const std::string& database, const std::string& id,
                                const model::ObjectVector_t &collisionObjects,
                                const std::string& heuristicName,
                                const bool loadValues, const bool disableEffectorCollision)
    {
        std::map<std::string, const sampling::heuristic>::const_iterator hit = checkLimbData(id, limbs_,factory_,heuristicName);;
        std::ifstream myfile (database.c_str());
        if (!myfile.good())
            throw std::runtime_error ("Impossible to open database");
        rbprm::RbPrmLimbPtr_t limb = rbprm::RbPrmLimb::create(device_, myfile, loadValues, hit->second, disableEffectorCollision);
        myfile.close();
        AddLimbPrivate(limb, id, limb->limb_->name(),collisionObjects, disableEffectorCollision);
    }

    void RbPrmFullBody::init(const RbPrmFullBodyWkPtr_t& weakPtr)
    {
        weakPtr_ = weakPtr;
    }

    RbPrmFullBody::RbPrmFullBody (const model::DevicePtr_t& device)
        : device_(device)
        , collisionValidation_(core::CollisionValidation::create(device))
        , weakPtr_(), noStability_ (true)
    {
        // NOTHING
    }

    // assumes unit direction
    std::vector<bool> setMaintainRotationConstraints(const fcl::Vec3f&) // direction)
    {
        std::vector<bool> res;
        for(std::size_t i =0; i <3; ++i)
        {
            res.push_back(true);
        }
        return res;
    }

    std::vector<bool> setRotationConstraints(const fcl::Vec3f&)// direction)
    {
        std::vector<bool> res;
        for(std::size_t i =0; i <2; ++i)
        {
            res.push_back(true);
        }
        res.push_back(false);
        return res;
    }


    std::vector<bool> setTranslationConstraints(const fcl::Vec3f&)// normal)
    {
        std::vector<bool> res;
        for(std::size_t i =0; i <3; ++i)
        {
            res.push_back(true);
        }
        return res;
    }

    bool ComputeCollisionFreeConfiguration
    (const hpp::rbprm::RbPrmFullBodyPtr_t& body,
     State& current,
     core::CollisionValidationPtr_t validation,
     const hpp::rbprm::RbPrmLimbPtr_t& limb,
     model::ConfigurationOut_t configuration,
     const double robustnessTreshold, bool stability )
    {
      sampling::T_OctreeReport finalSet;
      fcl::Contact contact;
      fcl::Vec3f normal;
      for(sampling::SampleVector_t::const_iterator cit = limb->sampleContainer_.samples_.begin();
          cit != limb->sampleContainer_.samples_.end(); ++cit){
        sampling::OctreeReport report(&(*cit),contact,(*cit).staticValue_,normal );
        finalSet.insert(report);
      }
      for(sampling::T_OctreeReport::const_iterator it = finalSet.begin() ; it != finalSet.end(); ++it){ // ordered by best static value ??
        sampling::Load(*it->sample_, configuration);
        hpp::core::ValidationReportPtr_t valRep (new hpp::core::CollisionValidationReport);
	hppDout (info, "classic validation test= " << validation->validate(configuration, valRep));
	core::CollisionValidationReport* report = static_cast <core::CollisionValidationReport*> (valRep.get());
	if (!report) {
	  hppDout (info, "collision with:"<<report->object1->name ());
	  hppDout (info, "collision with:"<<report->object2->name ());
	}
        if(validation->validate(configuration, valRep) && (!stability || stability::IsStable(body,current) >=robustnessTreshold))
        {
            current.configuration_ = configuration;
	    current = initializeExtraConfigs (body, current);
	    hppDout (info, "coll-free config= " << displayConfig(current.configuration_));
            return true;
        }
      }
      current = initializeExtraConfigs (body, current); // in case
      hppDout (info, "no coll-free config computed");
      return false;
      
      
      //////////////////////////////////////////////
      /*  for(std::vector<sampling::Sample>::const_iterator cit = limb->sampleContainer_.samples_.begin();
            cit != limb->sampleContainer_.samples_.end(); ++cit)
        {
            sampling::Load(*cit, configuration);
            hpp::core::ValidationReportPtr_t valRep (new hpp::core::CollisionValidationReport);
            if(validation->validate(configuration, valRep) && (!stability || stability::IsStable(body,current) >=robustnessTreshold))
            {
                current.configuration_ = configuration;
                return true;
            }
        }*/
        return false;
    }

    // first step
    State MaintainPreviousContacts(const State& previous, const hpp::rbprm::RbPrmFullBodyPtr_t& body,
                                   std::map<std::string,core::CollisionValidationPtr_t>& limbValidations,
                                   model::ConfigurationIn_t configuration, bool& contactMaintained, bool& multipleBreaks, const double robustnessTreshold, bool ignore6DOF = false)
    {
        contactMaintained = true;
        std::vector<std::string> brokenContacts;
        // iterate over every existing contact and try to maintain them
        State current;
        current.configuration_ = configuration;
        model::Configuration_t config = configuration;
        core::ConfigurationIn_t save = body->device_->currentConfiguration();
        // iterate over contact filo list
        std::queue<std::string> previousStack = previous.contactOrder_;
        while(!previousStack.empty())
        {
            const std::string name = previousStack.front();
            previousStack.pop();
            const RbPrmLimbPtr_t limb = body->GetLimbs().at(name);
            // try to maintain contact
            const fcl::Vec3f& ppos  =previous.contactPositions_.at(name);
            core::ConfigProjectorPtr_t proj = core::ConfigProjector::create(body->device_,"proj", 1e-3, 30);
            hpp::tools::LockJointRec(limb->limb_->name(), body->device_->rootJoint(), proj);
            const fcl::Vec3f z = limb->effector_->currentTransformation().getRotation() * limb->normal_;
            const fcl::Matrix3f& rotation = previous.contactRotation_.at(name);
            proj->add(core::NumericalConstraint::create (constraints::Position::create("",body->device_, limb->effector_,fcl::Vec3f(0,0,0), ppos)));
            if(limb->contactType_ == hpp::rbprm::_6_DOF && !ignore6DOF)
            {
                proj->add(core::NumericalConstraint::create (constraints::Orientation::create("rot_maintain_contact",body->device_,
                                                                                  limb->effector_,
                                                                                  rotation,
                                                                                  setMaintainRotationConstraints(z))));
            }
            if(proj->apply(config))
            {
                hpp::core::ValidationReportPtr_t valRep (new hpp::core::CollisionValidationReport);
                if(limbValidations.at(name)->validate(config, valRep) /*true*/)
                {
                    // stable?
                    current.contacts_[name] = true;
                    current.contactPositions_[name] = previous.contactPositions_.at(name);
                    current.contactNormals_[name] = previous.contactNormals_.at(name);
                    current.contactRotation_[name] = previous.contactRotation_.at(name);
                    current.contactOrder_.push(name);
                    current.configuration_ = config;
                }
                else
                {
                    contactMaintained = false;
                    ComputeCollisionFreeConfiguration(body,current,limbValidations.at(name),limb,current.configuration_,robustnessTreshold,false);
                    brokenContacts.push_back(name);
                }
            }
            else
            {
                contactMaintained = false;
                ComputeCollisionFreeConfiguration(body,current,limbValidations.at(name),limb,current.configuration_,robustnessTreshold,false);
                brokenContacts.push_back(name);
            }
        }
        // reload previous configuration
        body->device_->currentConfiguration(save);
        if(brokenContacts.size() > 1)
        {
            contactMaintained = false;
            multipleBreaks = true;
        }
        return current;
    }

    enum ContactComputationStatus
    {
      NO_CONTACT = 0,
      STABLE_CONTACT = 1,
      UNSTABLE_CONTACT = 2
      };

    bool ContactExistsWithinGroup(const hpp::rbprm::RbPrmLimbPtr_t& limb,
                                  const hpp::rbprm::RbPrmFullBody::T_LimbGroup& limbGroups,
                                  const State& current)
    {
        const std::vector<std::string>& group = limbGroups.at(limb->limb_->name());
        for(std::vector<std::string>::const_iterator cit = group.begin();
            cit != group.end(); ++cit)
        {
            if(current.contactPositions_.find(*cit) != current.contactPositions_.end())
            {
                return true;
            }
        }
        return false;
    }


    ContactComputationStatus ComputeStableContact(const hpp::rbprm::RbPrmFullBodyPtr_t& body,
                              State& current,
                              core::CollisionValidationPtr_t validation,
                              const std::string& limbId,
                              const hpp::rbprm::RbPrmLimbPtr_t& limb,
                              model::ConfigurationIn_t rbconfiguration,
                              model::ConfigurationOut_t configuration, const model::ObjectVector_t affordances,
                              const fcl::Vec3f& direction,
                              fcl::Vec3f& position, fcl::Vec3f& normal, const double robustnessTreshold,
                              bool contactIfFails = true, bool stableForOneContact = true,
                              const sampling::heuristic evaluate = 0)
    {
      hppDout (info, "compute stable contacts");
      // state already stable just find collision free configuration
      if(current.stable)
      {
         //if (ComputeCollisionFreeConfiguration(body, current, validation, limb, configuration)) return true;
      }
      fcl::Matrix3f rotation;
      sampling::T_OctreeReport finalSet;

      limb->limb_->robot()->currentConfiguration(rbconfiguration);
      limb->limb_->robot()->computeForwardKinematics ();
      fcl::Transform3f transform = limb->octreeRoot(); // get root transform from configuration


      //#pragma omp parallel for
      // request samples which collide with each of the collision objects
      
		
			sampling::heuristic eval = evaluate; if(!eval) eval =  limb->evaluate_;
      std::size_t i (0);
		  if (affordances.empty ()) {
		  	throw std::runtime_error ("No aff objects found!!!");
		  }
		  std::vector<sampling::T_OctreeReport> reports(affordances.size());
      for(model::ObjectVector_t::const_iterator oit = affordances.begin();
          oit != affordances.end(); ++oit, ++i)
      {
	hppDout (info, "collision obj= " << (*oit)->name ());
          if(eval)
            sampling::GetCandidates(limb->sampleContainer_, transform, *oit, direction, reports[i], eval);
          else
            sampling::GetCandidates(limb->sampleContainer_, transform, *oit, direction, reports[i]);
      }
      // order samples according to EFORT
      for(std::vector<sampling::T_OctreeReport>::const_iterator cit = reports.begin();
          cit != reports.end(); ++cit)
      {
	hppDout (info, "insert octree collision");
          finalSet.insert(cit->begin(), cit->end());
      }
      // pick first sample which is collision free
      bool found_sample(false);
      bool unstableContact(false); //set to true in case no stable contact is found
      core::Configuration_t moreRobust;
      double maxRob = -std::numeric_limits<double>::max();
      sampling::T_OctreeReport::const_iterator it = finalSet.begin();
      hppDout (info, "finalSet size: " << finalSet.size ());
      for(;!found_sample && it!=finalSet.end(); ++it)
      {
	//hppDout (info, "iteration on samples");
          const sampling::OctreeReport& bestReport = *it;
          sampling::Load(*bestReport.sample_, configuration);
          body->device_->currentConfiguration(configuration);
          {
              normal = bestReport.normal_;
              position = bestReport.contact_.pos;
              // the normal is given by the normal of the contacted object
              const fcl::Vec3f z = limb->effector_->currentTransformation().getRotation() * limb->normal_;
              const fcl::Matrix3f alignRotation = tools::GetRotationMatrix(z,normal);
              //rotation = alignRotation * limb->effector_->initialPosition ().getRotation();
	      rotation = alignRotation * limb->effector_->currentTransformation().getRotation();
              // Add constraints to resolve Ik
              core::ConfigProjectorPtr_t proj = core::ConfigProjector::create(body->device_,"proj", 1e-2, 30); // DEBUG, initially: 1e-4, 20
              // get current normal orientation
              hpp::tools::LockJointRec(limb->limb_->name(), body->device_->rootJoint(), proj);
              fcl::Vec3f posOffset = position - rotation * limb->offset_;
              posOffset = posOffset + normal * epsilon;
              fcl::Transform3f localFrame, globalFrame;
              globalFrame.setTranslation(posOffset);
	      //std::cout << "target " << globalFrame << std::endl;
	      proj->add(core::NumericalConstraint::create (constraints::Position::create("pos_stable_contact",body->device_, limb->effector_, localFrame, globalFrame, setTranslationConstraints(normal))));//

	      if(limb->contactType_ == hpp::rbprm::_6_DOF)
              {
		hppDout (info, "6 DOF contact type");
                  proj->add(core::NumericalConstraint::create (constraints::Orientation::create("rot_stable_contact",body->device_, limb->effector_, fcl::Transform3f(rotation),
                //localFrame.getRotation(),
                setRotationConstraints(z))));
              } else hppDout (info, "NOT 6 DOF contact type");
    
              if(proj->apply(configuration))
              {
		hppDout (info, "in projection test");
                hpp::core::ValidationReportPtr_t valRep (new hpp::core::CollisionValidationReport);
		hppDout (info, "tested config for contact= " << displayConfig(configuration));
                if(validation->validate(configuration, valRep)) // true
                {
		  hppDout (info, "config is valid");
		  /*#ifdef PROFILE
    watch.stop("collision");
    #endif*/
                    // test stability of new configuration
                    bool noStability = body->noStability_;
                    body->device_->currentConfiguration(configuration);
                    body->device_->computeForwardKinematics();
                    State tmp (current);
                    tmp.contacts_[limbId] = true;
                    tmp.contactPositions_[limbId] = limb->effector_->currentTransformation().getTranslation();
                    tmp.contactRotation_[limbId] = limb->effector_->currentTransformation().getRotation();
                    tmp.contactNormals_[limbId] = normal;
                    tmp.configuration_ = configuration;
                    ++tmp.nbContacts;
                    double robustness = stability::IsStable(body,tmp);
		    // tests cones Mylene
		    core::vector_t V0Dir = body->V0dir_;
		    core::vector_t VfDir = body->Vfdir_;
		    const std::size_t indexEC = body->device_->configSize () - body->device_->extraConfigSpace ().dimension ();
		    core::value_type mu = body->mu_;
		    bool inCone = true;
		    if (V0Dir.size () > 1 && VfDir.size () > 1 && indexEC >= 3) {
		      fcl::Vec3f n_fcl = updateConeDirection (tmp);
		      hppDout (info, "n_fcl= " << n_fcl);
		      core::vector_t n = n_fcl;
		      /*if (n.norm () < 1.1)
			inCone = isCentroidalConeValid (n, V0Dir, VfDir, mu);*/
		      hppDout (info, "inCone?= " << inCone);
		    } else
		      hppDout (info, "problem initializing isCentroidalConeValid");
		    
		    if (noStability) {
		      hppDout (info, "stability bypassed");
		      robustness = 4000000; // hardcoded to bypass stability
		    } else
		      hppDout (info, "stability computed");
                    if(((tmp.nbContacts == 1 && !stableForOneContact) || robustness>=robustnessTreshold) && inCone)
                    {
                        maxRob = std::max(robustness, maxRob);
                        position = limb->effector_->currentTransformation().getTranslation();
                        rotation = limb->effector_->currentTransformation().getRotation();
                        found_sample = true;
                    }
                    // if no stable candidate is found, select best contact
                    // anyway
                    else if((robustness > maxRob) && contactIfFails)
                    {
                        moreRobust = configuration;
                        maxRob = robustness;
                        position = limb->effector_->currentTransformation().getTranslation();
                        rotation = limb->effector_->currentTransformation().getRotation();
                        unstableContact = true;
                    }
                }else{
		  // Debug: check if collision and get colliding objects
		  core::ValidationReportPtr_t valReport;
		  validation->validate(configuration, valReport);
		  core::CollisionValidationReport* report =
                    static_cast <core::CollisionValidationReport*> (valReport.get());
		  hppDout (info, "config is NOT valid");
		  hppDout (info, "collision with:"<<report->object1->name ());
		  hppDout (info, "collision with:"<<report->object2->name ());
		}                
		/*#ifdef PROFILE
else
       watch.stop("collision");
       #endif*/
              }
	      /*#ifdef PROFILE
else
       watch.stop("ik");
       #endif*/
          }
      }

      ContactComputationStatus status(NO_CONTACT);
      if(found_sample)
      {
          status = STABLE_CONTACT;
          current.stable = true;
	  /*#ifdef PROFILE
    RbPrmProfiler& watch = getRbPrmProfiler();
    watch.add_to_count("contact", 1);
    #endif*/
      }
      else if(unstableContact)
      {          
          status = UNSTABLE_CONTACT;
          configuration = moreRobust;
	  /*#ifdef PROFILE
    RbPrmProfiler& watch = getRbPrmProfiler();
    watch.add_to_count("unstable contact", 1);
    #endif*/
      }
      else
      {
	/*#ifdef PROFILE
    RbPrmProfiler& watch = getRbPrmProfiler();
    watch.add_to_count("no contact", 1);
    #endif*/
          if(!found_sample)
          {
	    hppDout (info, "current.configuration_ (before ComputeCollisionFree)= " << displayConfig(current.configuration_));
              ComputeCollisionFreeConfiguration(body, current, validation, limb, configuration,robustnessTreshold,false);
	      hppDout (info, "sample not found, coll-free config");
	      hppDout (info, "current.configuration_= " << displayConfig(current.configuration_));
	      core::ValidationReportPtr_t valReport;
	      hppDout (info, "collision test= " << validation->validate(current.configuration_,valReport));
          }
      }
      if(found_sample || unstableContact)
      {
          current.contacts_[limbId] = true;
          current.contactNormals_[limbId] = normal;
          current.contactPositions_[limbId] = position;
          current.contactRotation_[limbId] = rotation;
          current.configuration_ = configuration;
          current.contactOrder_.push(limbId);
	  hppDout (info, "found sample or unstableContact -> update state");
	  hppDout (info, "current.configuration_= " << displayConfig(current.configuration_));
      }
      hppDout (info, "ComputeStableContact status= " << status);
      return status;
    }

		hpp::model::ObjectVector_t getAffObjectsForLimb(const std::string& limb,
			const affMap_t& affordances, const std::map<std::string, std::vector<std::string> >& affFilters)
		{
			model::ObjectVector_t affs;
			std::vector<std::string> affTypes;
			bool settingFound = false;
			for (std::map<std::string, std::vector<std::string> >::const_iterator fIt =
				affFilters.begin (); fIt != affFilters.end (); ++fIt) {
				std::size_t found = fIt->first.find(limb);
 			  if (found != std::string::npos) {
				affTypes = fIt->second;
				settingFound = true;
				break;
				}
			}
			if (!settingFound) {
        // TODO: Keep warning or delete it?
				std::cout << "No affordance filter setting found for limb " << limb
					<< ". Has such filter been set?" << std::endl;
				// Use all AFF OBJECTS as default if no filter setting exists
				for (affMap_t::const_iterator affordanceIt = affordances.begin ();
					affordanceIt != affordances.end (); ++affordanceIt) {
					std::copy (affordanceIt->second.begin (), affordanceIt->second.end (),
						std::back_inserter (affs));
				}
            } else {
				for (std::vector<std::string>::const_iterator affTypeIt = affTypes.begin ();
					affTypeIt != affTypes.end (); ++affTypeIt) {
					affMap_t::const_iterator affIt = affordances.find(*affTypeIt);
					std::copy (affIt->second.begin (), affIt->second.end (), std::back_inserter (affs));
				}
			}
			if (affs.empty()) {
			throw std::runtime_error ("No aff objects found for limb " + limb);
			}
			return affs;
		}

    bool RepositionContacts(State& result, const hpp::rbprm::RbPrmFullBodyPtr_t& body,
			core::CollisionValidationPtr_t validation,
      model::ConfigurationOut_t config, const affMap_t& affordances,
      const std::map<std::string, std::vector<std::string> >& affFilters, const fcl::Vec3f& direction,
			const double robustnessTreshold)
    {
        // replace existing contacts
        // start with older contact created
        std::stack<std::string> poppedContacts;
        std::queue<std::string> oldOrder = result.contactOrder_;
        std::queue<std::string> newOrder;
        core::Configuration_t savedConfig = config;
        std::string nContactName ="";
        State previous = result;
        while(!result.stable &&  !oldOrder.empty())
        {
            std::string previousContactName = oldOrder.front();
            std::string groupName = body->GetLimbs().at(previousContactName)->limb_->name();
            const std::vector<std::string>& group = body->GetGroups().at(groupName);
            oldOrder.pop();
            fcl::Vec3f normal, position;
            core::ConfigurationIn_t save = body->device_->currentConfiguration();
            bool notFound(true);
            for(std::vector<std::string>::const_iterator cit = group.begin();
                notFound && cit != group.end(); ++cit)
            {
								hpp::model::ObjectVector_t affs = getAffObjectsForLimb (*cit,
									affordances, affFilters);

                if(ComputeStableContact(body, result, validation, *cit, body->GetLimbs().at(*cit),save, config,
                	affs, direction, position, normal, robustnessTreshold, false)
                  == STABLE_CONTACT)
                {
                    nContactName = *cit;
                    notFound = false;
                }
                else
                {
                    result = previous;
                    config = savedConfig;
                }
            }
            if(notFound)
            {
                config = savedConfig;
                result.configuration_ = savedConfig;
                poppedContacts.push(previousContactName);
                body->device_->currentConfiguration(save);
            }
        }
        while(!poppedContacts.empty())
        {
            newOrder.push(poppedContacts.top());
            poppedContacts.pop();
        }
        while(!oldOrder.empty())
        {
            newOrder.push(oldOrder.front());
            oldOrder.pop();
        }
        if(result.stable)
        {
            newOrder.push(nContactName);
        }
        result.contactOrder_ = newOrder;
        return result.stable;
    }

    hpp::rbprm::State ComputeContacts
    (const hpp::rbprm::RbPrmFullBodyPtr_t& body,
     model::ConfigurationIn_t configuration, const affMap_t& affordances,
     const std::map<std::string, std::vector<std::string> >& affFilters, const fcl::Vec3f& direction,const double robustnessTreshold)
    {
      hppDout (info, "compute contacts");
        const T_Limb& limbs = body->GetLimbs();
        State result;
        // save old configuration
        core::ConfigurationIn_t save = body->device_->currentConfiguration();
        result.configuration_ = configuration;
        body->device_->currentConfiguration(configuration);
        body->device_->computeForwardKinematics();
        for(T_Limb::const_iterator lit = limbs.begin(); lit != limbs.end(); ++lit)
	  {
	    hppDout (info, "limb: " << lit->first);
            if(!ContactExistsWithinGroup(lit->second, body->limbGroups_ ,result))
            {
	      hpp::model::ObjectVector_t affs = getAffObjectsForLimb (lit->first, affordances, affFilters);
	      fcl::Vec3f normal, position;
	      ComputeStableContact(body,result, body->limbcollisionValidations_.at(lit->first), lit->first, lit->second, configuration, result.configuration_, affs, direction, position, normal, robustnessTreshold, true, false);
            }
            result.nbContacts = result.contactNormals_.size();
        }
        // reload previous configuration
        body->device_->currentConfiguration(save);
	result = initializeExtraConfigs (body, result);
        return result;
    }

    hpp::rbprm::State ComputeContacts
    (const hpp::rbprm::State& previous,
     const hpp::rbprm::RbPrmFullBodyPtr_t& body,
     model::ConfigurationIn_t configuration,
     const affMap_t& affordances,
     const std::map<std::string, std::vector<std::string> >& affFilters,
     const fcl::Vec3f& direction, bool& contactMaintained, bool& multipleBreaks,
     const bool allowFailure, const double robustnessTreshold)
    {
      hppDout (info, "compute contacts with previous state");
//static int id = 0;
    const T_Limb& limbs = body->GetLimbs();
    // save old configuration
    core::ConfigurationIn_t save = body->device_->currentConfiguration();
    model::Device::Computation_t flag = body->device_->computationFlag ();
    model::Device::Computation_t newflag = static_cast <model::Device::Computation_t> (model::Device::JOINT_POSITION);
    // load new root position
    body->device_->controlComputation (newflag);
    body->device_->currentConfiguration(configuration);
    body->device_->computeForwardKinematics ();
    // try to maintain previous contacts
    State result = MaintainPreviousContacts(previous,body, body->limbcollisionValidations_, configuration, contactMaintained, multipleBreaks, robustnessTreshold);
    // If more than one are broken, go back to previous state
    // and reposition
    hppDout(notice,"MaintainPreviousCOntact, contact maintained = "<<contactMaintained);
    if(multipleBreaks && !allowFailure)
      {
        fcl::Vec3f normal, position;
        result = previous;
        result.stable = false;
        std::string replaceContact =  result.RemoveFirstContact();
        model::Configuration_t config = previous.configuration_;
        body->device_->currentConfiguration(config);
        body->device_->computeForwardKinematics();
	hpp::model::ObjectVector_t affs = getAffObjectsForLimb (replaceContact,
								affordances, affFilters);
        // if no stable replacement contact found
        // modify contact order to try to replace another contact at the next step
        if(ComputeStableContact(body,result,
				body->limbcollisionValidations_.at(replaceContact),
				replaceContact,body->limbs_.at(replaceContact),
				configuration, config,affs,direction,
				position, normal, robustnessTreshold, true, false,
				body->factory_.heuristics_["random"]) != STABLE_CONTACT)
	  {
            result = previous;
            result.contactOrder_.pop();
            result.contactOrder_.push(replaceContact);
	  }
        body->device_->currentConfiguration(save);
        body->device_->controlComputation (flag);
	//++id;
        // in any case, returns state and raises failure flag
	result = initializeExtraConfigs (body, result);
        return result;
      }
    // no more than one contact was broken
    // we can go on normally
    core::Configuration_t config = result.configuration_;
    bool contactCreated(false);
    // iterate over each const free limb to try to generate contacts
    for(T_Limb::const_iterator lit = limbs.begin(); lit != limbs.end(); ++lit)
      {
        fcl::Vec3f normal, position;
        if(result.contacts_.find(lit->first) == result.contacts_.end()
           && !ContactExistsWithinGroup(lit->second, body->limbGroups_ ,result))
	  {
	    hpp::model::ObjectVector_t affs = getAffObjectsForLimb (lit->first,
								    affordances, affFilters);

            // if the contactMaintained flag remains true,
            // the contacts have not changed, and the state can be merged with the previous one eventually
            contactCreated = ComputeStableContact(body, result,
						  body->limbcollisionValidations_.at(lit->first), lit->first,
						  lit->second, configuration, config, affs, direction, position, normal,
						  robustnessTreshold) != NO_CONTACT || contactCreated;
	  }
      }
    hppDout(notice,"COntacts created = "<<contactCreated);
    contactMaintained = !contactCreated && contactMaintained;
    // reload previous configuration
    // no stable contact was found / limb maintained
    if(!result.stable && !body->noStability_)
      {
        // if no contact changes happened, try to modify one contact
        // existing previously to find a stable value
        if(contactMaintained)
	  {
            hppDout(notice,"result ! stable");
            contactMaintained = false;
            // could not reposition any contact. Planner has failed
            if (!RepositionContacts(result, body, body->collisionValidation_,
				    config, affordances, affFilters, direction, robustnessTreshold))
	      {
                std::cout << "planner is stuck; failure " <<  std::endl;
                body->device_->currentConfiguration(save);
                body->device_->controlComputation (flag);
                result.nbContacts = 0;
		result = initializeExtraConfigs (body, result);
                return result;
	      }
	  }
        // One contact break already happened, the state is invalid
        // if a new contact was created
        else
	  {
            fcl::Vec3f normal, position;
            result = previous;
            result.stable = false;
            std::string replaceContact =  result.RemoveFirstContact();
            if(!replaceContact.empty())
	      {
                model::Configuration_t config = previous.configuration_;
                body->device_->currentConfiguration(config);
                body->device_->computeForwardKinematics();
		hpp::model::ObjectVector_t affs = getAffObjectsForLimb (replaceContact,
									affordances, affFilters);

                // if a contact has already been created this iteration, or new contact is not stable
                // raise failure and switch contact order.
                if(contactCreated || 
		   ComputeStableContact(body,result,
					body->limbcollisionValidations_.at(replaceContact),replaceContact,
					body->limbs_.at(replaceContact), configuration, 
					config, affs, direction,position,
					normal,robustnessTreshold) != STABLE_CONTACT)
		  {
                    multipleBreaks = true;
                    result = previous;
                    result.contactOrder_.pop();
                    result.contactOrder_.push(replaceContact);
		  }
	      }
            body->device_->currentConfiguration(save);
            body->device_->controlComputation (flag);
	    //            ++id;
	    result = initializeExtraConfigs (body, result);
            return result;
	  }
      }
    body->device_->currentConfiguration(save);    
    body->device_->controlComputation (flag);
    result.nbContacts = result.contactNormals_.size();
    result = initializeExtraConfigs (body, result);
    return result;
    }
  } // rbprm
} //hpp
