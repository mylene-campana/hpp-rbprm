// Copyright (c) 2016, LAAS-CNRS
// Authors: Mylene Campana (mcampana@laas.fr), Steve Tonneau
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
#include <hpp/model/joint.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/core/locked-joint.hh>
#include <hpp/core/path.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/rbprm/rbprm-state.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-path.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-interpolation.hh>
#include <hpp/rbprm/stability/stability.hh>
#include <hpp/rbprm/ik-solver.hh>
#include <hpp/rbprm/interpolation/limb-rrt-helper.hh>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::value_type;
    using core::vector_t;
    using core::Configuration_t;

    BallisticInterpolationPtr_t
    BallisticInterpolation::create (const core::ProblemPtr_t &problem,
				    const RbPrmFullBodyPtr_t robot,
				    const State &start,
				    const State &end,
				    const core::PathVectorConstPtr_t path) {
      BallisticInterpolation* rbprmDevice =
	new BallisticInterpolation(problem, robot, start, end, path);
      BallisticInterpolationPtr_t res (rbprmDevice);
      res->init (res);
      return res;
    }

    BallisticInterpolation::BallisticInterpolation
    (const core::ProblemPtr_t& problem, const RbPrmFullBodyPtr_t robot,
     const hpp::rbprm::State &start, const State& end,
     const core::PathVectorConstPtr_t path) : 
      path_(path), start_(start), end_(end), problem_ (problem), robot_(robot) 
	
    {
      // TODO
      lastRootIndex_ = robot_->device_->configSize();
      for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
        hppDout(notice,"LIST OF LIMBS : "<< lit->first);
        if(lit->second->limb_->rankInConfiguration() < lastRootIndex_){
          lastRootIndex_ =lit->second->limb_->rankInConfiguration();
        }
      }
      hppDout(notice,"Last root index = "<<lastRootIndex_);
    }

    BallisticInterpolation::~BallisticInterpolation()
    {
      // NOTHING
    }

    // ========================================================================
    
    fcl::Vec3f BallisticInterpolation::computeDir (const vector_t V0,
						   const vector_t Vimp) {
      fcl::Vec3f dir;
      for (int i = 0; i < 3; i++) {
	//dir [i] = V0 [i] - Vimp [i];
	dir [i] = Vimp [i] - V0 [i]; // inversion for EFORT_normal
      }
      return dir;
    }

    std::size_t BallisticInterpolation::computeLimbLength
    (const model::JointPtr_t limb, const model::JointPtr_t effector)
    {
      std::size_t start = limb->rankInConfiguration();
      std::size_t end = effector->rankInConfiguration()
	+ effector->neutralConfiguration().rows();
      return end - start;
    }

    Configuration_t BallisticInterpolation::fillConfiguration
    (const Configuration_t config, const std::size_t configSize)
    {
      Configuration_t result (configSize);
      const std::size_t trunkSize = config.size ();
      const std::size_t ecsSize = 
	robot_->device_->extraConfigSpace ().dimension ();
      hppDout (info, "original config= " << displayConfig (config));
      for (std::size_t j = 0; j < configSize - ecsSize; j++) {
	if (j < trunkSize)
	  result [j] = config [j];
	else
	  result [j] = 0;
      }
      // copy extra-configs at the end of the config
      hppDout (info, "ecs size in fillConfiguration= " << ecsSize);
      for (std::size_t k = 0; k < ecsSize; k++)
	result [configSize - ecsSize + k] = config [trunkSize - ecsSize + k];
      hppDout (info, "filled config= " << displayConfig (result));
      return result;
    }

    Configuration_t BallisticInterpolation::fillConfiguration
    (const Configuration_t trunkConfig, const Configuration_t refConfig)
    {
      Configuration_t result (refConfig.rows ());
      const std::size_t configSize = refConfig.rows ();
      //const std::size_t trunkSize = trunkConfig.size ();
      const std::size_t trunkSize = 7;
      const std::size_t ecsSize = 
	robot_->device_->extraConfigSpace ().dimension ();
      hppDout (info, "original config= " << displayConfig (trunkConfig));
      for (std::size_t j = 0; j < configSize - ecsSize; j++) {
	if (j < trunkSize)
	  result [j] = trunkConfig [j];
	else
	  result [j] = refConfig [j];
      }
      // copy extra-configs of trunk at the end of the config
      hppDout (info, "ecs size in fillConfiguration= " << ecsSize);
      for (std::size_t k = 0; k < ecsSize; k++)
	result[configSize - ecsSize + k] = trunkConfig[trunkSize - ecsSize + k];
      hppDout (info, "filled config= " << displayConfig (result));
      return result;
    }

    Configuration_t BallisticInterpolation::replaceLimbConfigsInFullConfig
    (const Configuration_t q_full, const Configuration_t q_contact,
     const std::vector <RbPrmLimbPtr_t> limbs) {
      Configuration_t result = q_full;
      hppDout (info, "q_full= " << displayConfig(q_full));
      hppDout (info, "q_contact= " << displayConfig(q_contact));
      const core::DevicePtr_t device = robot_->device_;
      for (std::size_t i = 0; i < limbs.size (); i++ ){
	  const RbPrmLimbPtr_t limb = limbs.at (i);
	  hppDout (info, "limb (in replaceLimb)= " << limb->limb_->name ());
	  model::JointPtr_t effectorClone = device->getJointByName(limb->effector_->name ());
	  const std::size_t startRank (limb->limb_->rankInConfiguration());
	  const std::size_t length (computeLimbLength (limb->limb_, effectorClone));
	  hppDout (info, "startRank= " << startRank);
	  hppDout (info, "result (before replaceLimb)= " << displayConfig(result));
	  for (std::size_t i = 0; i < length; i++) {
	    std::size_t currentRank = startRank + i;
	    result (currentRank) = q_contact (currentRank);
	  }
	  hppDout (info, "result= " << displayConfig(result));
      }
      return result;
    }

    std::vector <RbPrmLimbPtr_t> BallisticInterpolation::activeLimbsFromROMs
    (const ParabolaPathPtr_t pp, const bool startConfig) {
      std::vector <RbPrmLimbPtr_t> activeLimbs;
      const T_Limb& limbs = robot_->GetLimbs();
      std::vector <std::string> ROMnames;
      if (startConfig)
	ROMnames = pp->initialROMnames_;
      else
	ROMnames = pp->endROMnames_;

      hppDout (info, "ROMnames size= " << ROMnames.size ());
      for (T_Limb::const_iterator it = limbs.begin(); it != limbs.end(); it++) {
	const RbPrmLimbPtr_t limb = it->second;
	const std::string effectorName = limb->effector_->name ();
	hppDout (info, "effectorName= " << effectorName);
	for (std::size_t k = 0; k < ROMnames.size (); k++) {
	  hppDout (info, "ROMnames.at(k)= " << ROMnames.at(k));
	  if (effectorName.compare (ROMnames.at(k)) == 0) {
	    hppDout (info, "found in common in name lists= " << effectorName);
	    activeLimbs.push_back (limb);
	  }
	}//for roms
      }//for limbs
      hppDout (info, "nb active limbs= " << activeLimbs.size ());
      return activeLimbs;
    }

    State BallisticInterpolation::computeOffsetContactConfig
    (const BallisticPathPtr_t bp,
     const State& previousState,State& transitionDOFstate, const value_type u_offset,
     const bool increase_u_offset,value_type& lenght,value_type& lenghtTransition,
     const std::size_t maxIter, const value_type alpha) {
      Configuration_t q_trunk_offset, q_contact_offset(previousState.configuration_), q_interp, result;
      q_trunk_offset = previousState.configuration_; 
      hppDout(notice,"q_trunk_offset = "<<displayConfig(q_trunk_offset));
      core::DevicePtr_t robot = robot_->device_;
      bool success, contact_OK = true, multipleBreaks, contactMaintained;
      bool ignore6DOF = false;
      std::size_t iteration = 0;
      fcl::Vec3f dir;
      State state,lastState;
      lastState = previousState;
      value_type currentLenght;
      value_type u;
      if(increase_u_offset){ // forward progression
        u = u_offset;
        currentLenght = 0;
      }else{ // backward progression
        u = - u_offset;
        currentLenght = bp->length();
      }
      std::map<std::string,core::CollisionValidationPtr_t> limbColVal = 
          robot_->getLimbcollisionValidations ();
      
      q_interp = (*bp) (currentLenght, success);
      
      if (u_offset == 0) {
        hppDout (info, "no cushion effect asked, return interpolated config");
        state.configuration_ = q_interp;
        return state;
      }
      
      while (contact_OK && iteration < maxIter && ((((currentLenght<bp->length())/3.)&&increase_u_offset) || ((currentLenght > (bp->length()*2./3.))&&(!increase_u_offset)))){ 
        hppDout (info, "currentLenght= " << currentLenght);        
        currentLenght += u;        
        iteration++;
        hppDout (info, "iteration= " << iteration);
        hppDout (info, "currentLenght = " << currentLenght);
        q_trunk_offset = (*bp) (currentLenght, success);
        dir = bp->evaluateVelocity (currentLenght);
        //state = MaintainPreviousContacts (lastState, limbColVal, q_trunk_offset, contactMaintained, multipleBreaks, successLimbs);
        //state = robot_->MaintainPreviousContacts (lastState, robot_, limbColVal, q_trunk_offset, contactMaintained, multipleBreaks,0.);
        state = rbprm::ComputeContacts(lastState,robot_,q_trunk_offset,problem_->collisionObstacles(),dir,contactMaintained,multipleBreaks,true,0.,ignore6DOF,false);
        /*if(!contactMaintained && !ignore6DOF){ // after the first fail with rotationnal constraint, we relax the problem for longer path
          ignore6DOF = true;
          transitionDOFstate = state;
          if(increase_u_offset)
            lenghtTransition = currentLenght-u;
          else
            lenghtTransition = bp->length()-currentLenght+u;
          state = rbprm::ComputeContacts(lastState,robot_,q_trunk_offset,problem_->collisionObstacles(),dir,contactMaintained,multipleBreaks,true,0.,ignore6DOF,false);
          hppDout(notice,"Relax 6DOF constraints after "<<iteration<<" iterations");
        }*/
        if(contactMaintained){
          contact_OK = true;
          lastState = state;
        }else{
          contact_OK = false;
        }
        hppDout (info, "q_contact_offset= " << displayConfig (q_contact_offset));
        hppDout (info, "## contactMaintained = " << contactMaintained);
        
      }//while
      
      // replace limbs that are accurate:
      // take limb-configs from contact if contact succeeded, 
      // from interp otherwise.
      //result = replaceLimbConfigsInFullConfig (q_interp,                                               q_contact_offset, fixedLimbs);
      
      // replace the limbs not used for contact with their configuration in flexionPose_
      size_t minIndex = robot_->device_->configSize();
      
      for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
        hppDout(notice,"LIST OF LIMBS : "<< lit->first << "contact = "<<lastState.contacts_[lit->first]);
        if(lit->second->limb_->rankInConfiguration() < minIndex){
          minIndex =lit->second->limb_->rankInConfiguration();
        }
        if( ! lastState.contacts_[lit->first]){ // limb is not in contact
          hppDout(notice," Not in contact, index config : "<<lit->second->limb_->rankInConfiguration()<<" -> "<<lit->second->effector_->rankInConfiguration());
          for(size_t i = lit->second->limb_->rankInConfiguration() ; i < lit->second->effector_->rankInConfiguration() ; i++){
            lastState.configuration_[i] = flexionPose_[i];
          }
        } 
      }
      
      // replace the trunkDOF with contactPose value : (we suppose we always work with freeflyer as root ....)
      for (size_t i = 7 ; i < minIndex ; i++){
        lastState.configuration_[i] = flexionPose_[i];
      }
      if(increase_u_offset)
        lenght = currentLenght-u;
      else
        lenght = bp->length()-currentLenght+u;


      lastState.ignore6DOF = ignore6DOF;
      return lastState;
    }

    State BallisticInterpolation::computeTopExtendingPose 
    (const core::PathPtr_t path, const BallisticPathPtr_t bp, value_type& lenght) {
      bool success;
      ParabolaPathPtr_t pp = boost::dynamic_pointer_cast<ParabolaPath>(path);
      State state;
      const vector_t coefs = pp->coefficients ();
      const value_type pathLength = path->length();
      const Configuration_t init = (*path) (0, success);
      const Configuration_t end = (*path) (pathLength, success);
      const value_type theta = coefs (3);
      const value_type x_theta_max = - 0.5 * coefs (1) / coefs (0);
      const value_type x_theta_init = cos(theta)*init (0) + sin(theta)*init (1);
      const value_type x_theta_end = cos(theta)*end (0) + sin(theta)*end (1);
      const value_type u_max = (x_theta_max - x_theta_init)
	/ (x_theta_end - x_theta_init); // in [0,1]
      const value_type z_x_theta_max = coefs (0)*x_theta_max*x_theta_max +
	coefs (1)*x_theta_max + coefs (2);
      core::Configuration_t q_top;
      const std::string robotName = robot_->device_->name ();
      hppDout (info, "robotName= " << robotName);
      const bool parabTallEnough = (robotName.compare ("ant") == 0 &&  z_x_theta_max > 0.44) || (robotName.compare ("spiderman") == 0 &&  z_x_theta_max > 1.1) || (robotName.compare ("frog") == 0 &&  z_x_theta_max > 0.2);
      value_type r = 0.5; // blending coefficient of extending key-frame
      const Configuration_t q_interp_top = (*bp) (u_max*pathLength, success);
      lenght = u_max*pathLength;
      if (parabTallEnough)
	r = 0.8;
      else
	r = 0.4;
      hppDout (info, "r of blending extending pose= " << r);
      
      if (extendingPose_.rows ()) {
	const core::Configuration_t q_trunk_max = (*path) (u_max*pathLength, success);
	const core::Configuration_t q = fillConfiguration (q_trunk_max, extendingPose_); // now q is extending_ at the good top-trunk-configuration
	q_top = blendPoses (q, q_interp_top, r);
      } else {
	q_top = q_interp_top;
	hppDout (info, "no extending pose");
      }
      hppDout (info, "q_top= " << displayConfig(q_top));
      
      state.configuration_ = q_top;
      // test : 
     /* for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
        hppDout(notice,"LIST OF LIMBS : "<< lit->first << "contact = "<<state.contacts_[lit->first]);
      }*/
      
      return state;
    }
    
    
    Configuration_t computeContactPose(const State& state, Configuration_t contactPose, rbprm::RbPrmFullBodyPtr_t robot){
      Configuration_t q = state.configuration_;
      bool useAllLimb = true;
      hppDout(info, "contact configuration = "<<displayConfig(contactPose));
      size_t minIndex = robot->device_->configSize();
      // replace the limbs not used for contact with their configuration in flexionPose_
      for( rbprm::T_Limb::const_iterator lit = robot->GetLimbs().begin();lit != robot->GetLimbs().end(); ++lit){
        hppDout(notice,"Contact Pose, LIST OF LIMBS  : "<< lit->first << "contact = "<<(state.contacts_.find(lit->first) != state.contacts_.end()));
        if(lit->second->limb_->rankInConfiguration() < minIndex){
          hppDout(notice," Min index = "<<lit->second->limb_->rankInConfiguration()) ;         
          minIndex = lit->second->limb_->rankInConfiguration();
        }
        if( (state.contacts_.find(lit->first) == state.contacts_.end())){ // limb is not in contact
          useAllLimb = false;
          hppDout(notice," Not in contact, index config : "<<lit->second->limb_->rankInConfiguration()<<" -> "<<lit->second->effector_->rankInConfiguration());
          for(size_t i = lit->second->limb_->rankInConfiguration() ; i < lit->second->effector_->rankInConfiguration() ; i++){
            q[i] = contactPose[i];
          }
        } 
      }
      
     
      
      // replace the trunkDOF with contactPose value : ( 7 because we suppose we always work with freeflyer as root ....)
      if(!useAllLimb){
        for (size_t i = 7 ; i < minIndex ; i++){
          q[i] = contactPose[i];
        }
      }
      return q;
    }
    
    Configuration_t BallisticInterpolation::computeContactPose(const State& state){
      if(contactPose_.size() == 0)
        return state.configuration_;
      else
        return rbprm::computeContactPose(state,contactPose_,robot_);
    }

    Configuration_t BallisticInterpolation::blendPoses 
    (const Configuration_t q1, const Configuration_t q2, const value_type r)
    {
      Configuration_t result = r*q1 + (1-r)*q2;
      return result;
    }

    // ========================================================================
    
    /*std::vector<Configuration_t>
    BallisticInterpolation::InterpolateConfigs (const double timeStep)
    {
      if(!path_) throw std::runtime_error ("Can not interpolate; not path given to interpolator ");
      std::vector<Configuration_t> configs;
      bool success;
      Configuration_t qStart = start_.configuration_;
      Configuration_t qEnd = end_.configuration_;
      core::DevicePtr_t robot = robot_->device_; // device of fullbody
      const std::size_t subPathNumber = path_->numberPaths ();
      hppDout (info, "number of sub-paths: " << subPathNumber);
      core::PathVectorPtr_t newPath = core::PathVector::create 
	(robot->configSize (), robot->numberDof ());

      for (std::size_t i=0; i<subPathNumber; i++) {
	core::PathPtr_t subpath = path_->pathAtRank (i);
	Configuration_t q1 = subpath->initial ();
	Configuration_t q2 = subpath->end ();
	BallisticPathPtr_t bp = Interpolate (q1, q2, subpath->length (),
					     subpath->coefficients ());
	const core::interval_t& range = path_->timeRange();
	for(value_type t = range.first; t< range.second; t += timeStep)
	  {
	    configs.push_back( bp->operator () (t, success));
	  }
      }
      return configs;
      }*/

    core::PathVectorPtr_t BallisticInterpolation::InterpolateFullPath
    (const core::value_type u_offset) {
      if(!path_) throw std::runtime_error ("Cannot interpolate; not path given to interpolator ");
      Configuration_t qStart = computeContactPose(start_);
      Configuration_t qEnd = computeContactPose(end_);
      core::DevicePtr_t robot = robot_->device_;
      T_StateFrame stateFrames;      
      Configuration_t q2(robot->configSize ()),
	q1contact (robot->configSize ()), q2contact (robot->configSize ());
      BallisticPathPtr_t bp;
      const model::ObjectVector_t &collisionObjects =
	problem_->collisionObstacles();
      hppDout (info, "u_offset= " << u_offset);
      const std::size_t subPathNumber = path_->numberPaths ();
      hppDout (info, "number of sub-paths: " << subPathNumber);
      core::PathVectorPtr_t newPath = core::PathVector::create 
	(robot->configSize (), robot->numberDof ());
      robot_->noStability_ = true; // disable stability for waypoints
      vector_t V0 (3), Vimp (3); fcl::Vec3f dir;
      core::PathPtr_t subpath = path_->pathAtRank (0);
      State state1, state2;
      State contactState1,contactState2,stateTop,contactTransition1,contactTransition2;
      core::value_type lenghtTop,lenghtTakeoff,lenghtLanding,lenghtLanding6DOF,lenghtTakeoff6DOF;


      for (std::size_t i = 0; i < subPathNumber - 1; i++) {
  stateFrames.clear();      
	hppDout (info, "B-interp on path nb: " << i);
	core::PathPtr_t subpath_next = path_->pathAtRank (i+1);
	ParabolaPathPtr_t pp = 
	  boost::dynamic_pointer_cast<ParabolaPath>(subpath);
	ParabolaPathPtr_t pp_next = 
	  boost::dynamic_pointer_cast<ParabolaPath>(subpath_next);
	const value_type pathLength = subpath->length ();
  core::PathPtr_t pathLimb;
	
	if (i == 0) { // keep qStart config which already has contacts
	  hppDout (info, "keep start config");
	  q1contact = qStart;
	  state1 = start_;
	}
	else {
	  q1contact = q2contact;
	  state1 = state2;
	}
	q2 = fillConfiguration (subpath->end (), robot->configSize ());
	hppDout (info, "q2: " << displayConfig(q2));
	V0 = pp_next->V0_; // V0_i+1
	Vimp = pp->Vimp_; // Vimp_i
	dir = computeDir (V0, Vimp);
	hppDout (info, "dir (Vimp-V0)= " << dir);
	state2 = ComputeContacts(robot_, q2, collisionObjects, dir);
	q2contact = computeContactPose(state2);
	// compute average-normal corresponding to new contacts
	std::queue<std::string> contactStack = state2.contactOrder_;
	fcl::Vec3f normalAv = (0,0,0);
	const std::size_t contactNumber = contactStack.size ();
	while(!contactStack.empty())
        {
	  const std::string name = contactStack.front();
	  contactStack.pop();
	  const fcl::Vec3f& normal = state2.contactNormals_.at(name);
	  for (std::size_t j = 0; j < 3; j++)
	    normalAv [j] += normal [j]/contactNumber;
	}
	normalAv.normalize ();
	hppDout (info, "normed normalAv= " << normalAv);
	// If robot has ECS, fill new average-normal in it
	if (robot_->device_->extraConfigSpace ().dimension () > 3 && normalAv.norm () > 0.9) {
	  const std::size_t indexECS = robot_->device_->configSize () - robot_->device_->extraConfigSpace ().dimension ();
	  for (std::size_t i = 0; i < 3; i++)
	    q2contact [indexECS + i] = normalAv [i];
	}
	hppDout (info, "q2contact= " << displayConfig(q2contact));
	hppDout (info, "q1contact= " << displayConfig(q1contact));

	bp = Interpolate (q1contact, q2contact, pathLength, subpath->coefficients ());
  bp->lastRootIndex(lastRootIndex_);
	stateTop = computeTopExtendingPose (subpath, bp,lenghtTop);

	//bp1max = Interpolate (q1contact, q_max,    lenghtTop, subpath->coefficients ());
	//bp2max = Interpolate (q_max, q2contact,  bp->length()-lenghtTop,subpath->coefficients ());

    contactState1 = computeOffsetContactConfig (bp, state1,contactTransition1, u_offset, true,lenghtTakeoff,lenghtTakeoff6DOF);
    contactState2 = computeOffsetContactConfig (bp, state2,contactTransition2, u_offset, false,lenghtLanding,lenghtLanding6DOF);

  /*
	bp1 = Interpolate (q1contact, q_contact_offset1,
			   lenghtTakeoff,
			   subpath->coefficients ());
	bp1max = Interpolate (q_contact_offset1, q_max, 
			      lenghtTop-lenghtTakeoff,
			      subpath->coefficients ());
	bp2max = Interpolate (q_max, q_contact_offset2,
			      bp->length()-lenghtTop-lenghtLanding,
			      subpath->coefficients ());
	bp3 = Interpolate (q_contact_offset2, q2contact,
			   lenghtLanding,
			   subpath->coefficients ());

	newPath->appendPath (bp1);
	newPath->appendPath (bp1max);
	newPath->appendPath (bp2max);
	newPath->appendPath (bp3);
  */
    stateFrames.clear();
    stateFrames.push_back(std::make_pair(0,state1));
   /* if(contactState1.ignore6DOF && (lenghtTakeoff != lenghtTakeoff6DOF))
        stateFrames.push_back(std::make_pair(lenghtTakeoff6DOF,contactTransition1));*/
    stateFrames.push_back(std::make_pair(lenghtTakeoff,contactState1));
    stateFrames.push_back(std::make_pair(lenghtTop,stateTop));
    stateFrames.push_back(std::make_pair(bp->length() - lenghtLanding,contactState2));
   /* if(contactState2.ignore6DOF && (lenghtLanding != lenghtLanding6DOF))
        stateFrames.push_back(std::make_pair(bp->length() -lenghtLanding6DOF,contactTransition2));*/
    stateFrames.push_back(std::make_pair(bp->length(),state2));
  
    
    hppDout(notice, "position initial state frame  = "<<displayConfig(state1.configuration_));
    hppDout(notice, "position initial Contact transition state frame  = "<<displayConfig(contactTransition1.configuration_));
    hppDout(notice, "position initial Contact state frame  = "<<displayConfig(contactState1.configuration_));
    hppDout(notice, "position top frame  = "<<displayConfig(stateTop.configuration_));
    hppDout(notice, "position final contact frame  = "<<displayConfig(contactState2.configuration_));
    hppDout(notice, "position final Contact transition state frame  = "<<displayConfig(contactTransition2.configuration_));
    hppDout(notice, "position final state frame  = "<<displayConfig(state2.configuration_));
    hppDout(notice, "TIME initial state frame  = "<<0);
    hppDout(notice, "TIME initial Contact transition state frame  = "<<lenghtTakeoff6DOF);
    hppDout(notice, "TIME initial Contact state frame  = "<<lenghtTakeoff);
    hppDout(notice, "TIME top frame  = "<<lenghtTop);
    hppDout(notice, "TIME final contact frame  = "<<bp->length() - lenghtLanding);
    hppDout(notice, "TIME final contact transition frame  = "<<bp->length() - lenghtLanding6DOF);
    hppDout(notice, "TIME final state frame  = "<<bp->length());
    hppDout(notice,"test last root index interpolate = "<<bp->lastRootIndex());
    
    
  pathLimb = rbprm::interpolation::interpolateStates(robot_,problem_,bp,stateFrames.begin(),stateFrames.end()-1,2);
  newPath->appendPath(pathLimb);
  
	subpath = subpath_next;
	pp = boost::dynamic_pointer_cast<ParabolaPath>(subpath);
	
	if (i == subPathNumber - 2) { // subpath_next = final subpath
	  q1contact = q2contact;
	  state1 = state2;
	  q2contact = qEnd;
	  state2 = end_;

	  // Compute top-keyframe and adapt ballistic-path with it
	  bp = Interpolate (q1contact, q2contact, subpath->length (),
			    subpath->coefficients ());
    bp->lastRootIndex(lastRootIndex_);
	  stateTop = computeTopExtendingPose (subpath, bp,lenghtTop);
	  //bp1max = Interpolate (q1contact, q_max,	lenghtTop,	subpath->coefficients ());
	  //bp2max = Interpolate (q_max, q2contact,	bp->length()-lenghtTop,	subpath->coefficients ());

      contactState1 = computeOffsetContactConfig (bp, state1,contactTransition1, u_offset, true,lenghtTakeoff,lenghtTakeoff6DOF);
      contactState2 = computeOffsetContactConfig (bp, state2,contactTransition2 ,u_offset, false,lenghtLanding,lenghtLanding6DOF);
/*
	  bp1 = Interpolate (q1contact, q_contact_offset1,
			     lenghtTakeoff,
			     subpath->coefficients ());
	  bp1max = Interpolate (q_contact_offset1, q_max,
				lenghtTop-lenghtTakeoff,
				subpath->coefficients ());
	  bp2max = Interpolate (q_max, q_contact_offset2,
				bp->length()-lenghtTop-lenghtLanding,
				subpath->coefficients ());
	  bp3 = Interpolate (q_contact_offset2, q2contact,
			     lenghtLanding,
			     subpath->coefficients ());

	  newPath->appendPath (bp1);
	  newPath->appendPath (bp1max);
	  newPath->appendPath (bp2max);
	  newPath->appendPath (bp3);*/
      stateFrames.clear();
      stateFrames.push_back(std::make_pair(0,state1));
     /* if(contactState1.ignore6DOF && (lenghtTakeoff != lenghtTakeoff6DOF))
          stateFrames.push_back(std::make_pair(lenghtTakeoff6DOF,contactTransition1));*/
      stateFrames.push_back(std::make_pair(lenghtTakeoff,contactState1));
      stateFrames.push_back(std::make_pair(lenghtTop,stateTop));
      stateFrames.push_back(std::make_pair(bp->length() - lenghtLanding,contactState2));
     /* if(contactState2.ignore6DOF && (lenghtLanding != lenghtLanding6DOF))
          stateFrames.push_back(std::make_pair(bp->length() -lenghtLanding6DOF,contactTransition2));*/
      stateFrames.push_back(std::make_pair(bp->length(),state2));
      
      hppDout(notice, "position initial state frame  = "<<displayConfig(state1.configuration_));
      hppDout(notice, "position initial Contact transition state frame  = "<<displayConfig(contactTransition1.configuration_));
      hppDout(notice, "position initial Contact state frame  = "<<displayConfig(contactState1.configuration_));
      hppDout(notice, "position top frame  = "<<displayConfig(stateTop.configuration_));
      hppDout(notice, "position final contact frame  = "<<displayConfig(contactState2.configuration_));
      hppDout(notice, "position final Contact transition state frame  = "<<displayConfig(contactTransition2.configuration_));
      hppDout(notice, "position final state frame  = "<<displayConfig(state2.configuration_));
      hppDout(notice, "TIME initial state frame  = "<<0);
      hppDout(notice, "TIME initial Contact transition state frame  = "<<lenghtTakeoff6DOF);
      hppDout(notice, "TIME initial Contact state frame  = "<<lenghtTakeoff);
      hppDout(notice, "TIME top frame  = "<<lenghtTop);
      hppDout(notice, "TIME final contact frame  = "<<bp->length() - lenghtLanding);
      hppDout(notice, "TIME final contact transition frame  = "<<bp->length() - lenghtLanding6DOF);
      hppDout(notice, "TIME final state frame  = "<<bp->length());
      hppDout(notice,"test last root index interpolate = "<<bp->lastRootIndex());
    
    pathLimb = rbprm::interpolation::interpolateStates(robot_,problem_,bp,stateFrames.begin(),stateFrames.end()-1,2);
    newPath->appendPath(pathLimb);
    
	}//if final subpath
      }// for subpaths
      return newPath;
    }

    core::PathVectorPtr_t BallisticInterpolation::InterpolateDirectPath
    (const core::value_type u_offset) {
      if(!path_) throw std::runtime_error ("Cannot interpolate; not path given to interpolator ");
      hppDout (info, "direct B-interpolation");
      Configuration_t qStart = computeContactPose(start_);
      Configuration_t qEnd = computeContactPose(end_);
      core::DevicePtr_t robot = robot_->device_;
      BallisticPathPtr_t bp;
      core::PathPtr_t pathLimb;
      core::PathVectorPtr_t newPath = core::PathVector::create 
	(robot->configSize (), robot->numberDof ());
      robot_->noStability_ = true; // disable stability for waypoints
      vector_t V0 (3), Vimp (3); fcl::Vec3f dir;
      State contactState1,contactState2,contactTransition1,contactTransition2,stateTop;
      T_StateFrame stateFrames;
      value_type lenghtTop,lenghtTakeoff,lenghtLanding,lenghtLanding6DOF,lenghtTakeoff6DOF;
      const core::PathPtr_t path = path_->pathAtRank (0);
      const value_type pathLength = path->length ();
      const vector_t pathCoefs = path->coefficients ();

      bp = Interpolate (qStart, qEnd, pathLength, pathCoefs);
      stateTop = computeTopExtendingPose (path, bp,lenghtTop);
      bp->lastRootIndex(lastRootIndex_);
      hppDout(notice, "full lenght = "<<bp->length());
      hppDout(notice,"lenght top = "<<lenghtTop);
      //bp1max = Interpolate (qStart, q_max, lenghtTop, pathCoefs);
      //bp2max = Interpolate (q_max, qEnd, bp->length()-lenghtTop,  pathCoefs);

      contactState1 = computeOffsetContactConfig (bp, start_, contactTransition1,u_offset, true,lenghtTakeoff,lenghtTakeoff6DOF);
      contactState2 = computeOffsetContactConfig (bp, end_, contactTransition2,u_offset, false,lenghtLanding,lenghtLanding6DOF);

      /*
      bp1 = Interpolate (qStart, q_contact_offset1, lenghtTakeoff, pathCoefs);
      hppDout(notice,"BP1 lenght = "<<bp1->length());
      bp1max = Interpolate (q_contact_offset1, q_max, lenghtTop-lenghtTakeoff,pathCoefs);
      bp2max = Interpolate (q_max, q_contact_offset2, bp->length()-lenghtTop-lenghtLanding,	pathCoefs);
      bp3 = Interpolate (q_contact_offset2, qEnd, lenghtLanding, pathCoefs);
      
      hppDout(notice,"lenght bp1max = "<<bp1max->length());
      hppDout(notice,"lenght bp2max = "<<bp2max->length());
      
      hppDout(notice,"BP3 lenght = "<<bp3->length());
      newPath->appendPath (bp1);
      newPath->appendPath (bp1max);
      newPath->appendPath (bp2max);
      newPath->appendPath (bp3);
      */
      stateFrames.push_back(std::make_pair(0,start_));
      /*if(contactState1.ignore6DOF && (lenghtTakeoff != lenghtTakeoff6DOF))
          stateFrames.push_back(std::make_pair(lenghtTakeoff6DOF,contactTransition1));*/
      stateFrames.push_back(std::make_pair(lenghtTakeoff,contactState1));
      stateFrames.push_back(std::make_pair(lenghtTop,stateTop));
      stateFrames.push_back(std::make_pair(bp->length() - lenghtLanding,contactState2));
      /*if(contactState2.ignore6DOF && (lenghtLanding != lenghtLanding6DOF)){
          contactState2.ignore6DOF = false;
          contactTransition2.ignore6DOF = true;
          stateFrames.push_back(std::make_pair(bp->length() -lenghtLanding6DOF,contactTransition2));
      }*/
      stateFrames.push_back(std::make_pair(bp->length(),end_));
      hppDout(notice, "position initial state frame  = "<<displayConfig(start_.configuration_));
      hppDout(notice, "position initial Contact transition state frame  = "<<displayConfig(contactTransition1.configuration_));
      hppDout(notice, "position initial Contact state frame  = "<<displayConfig(contactState1.configuration_));
      hppDout(notice, "position top frame  = "<<displayConfig(stateTop.configuration_));
      hppDout(notice, "position final contact frame  = "<<displayConfig(contactState2.configuration_));
      hppDout(notice, "position final Contact transition state frame  = "<<displayConfig(contactTransition2.configuration_));
      hppDout(notice, "position final state frame  = "<<displayConfig(end_.configuration_));
      hppDout(notice, "TIME initial state frame  = "<<0);
      hppDout(notice, "TIME initial Contact transition state frame  = "<<lenghtTakeoff6DOF);
      hppDout(notice, "TIME initial Contact state frame  = "<<lenghtTakeoff);
      hppDout(notice, "TIME top frame  = "<<lenghtTop);
      hppDout(notice, "TIME final contact frame  = "<<bp->length() - lenghtLanding);
      hppDout(notice, "TIME final contact transition frame  = "<<bp->length() - lenghtLanding6DOF);
      hppDout(notice, "TIME final state frame  = "<<bp->length());
      hppDout(notice,"test last root index interpolate = "<<bp->lastRootIndex());
      pathLimb = rbprm::interpolation::interpolateStates(robot_,problem_,bp,stateFrames.begin(),stateFrames.end()-1,2);
    //  bp->setLimbPath(pathLimb);
      newPath->appendPath(pathLimb);
            
      return newPath;
    }

    BallisticPathPtr_t BallisticInterpolation::Interpolate 
    (const Configuration_t q1, const Configuration_t q2,
     const core::value_type length, const core::vector_t coefficients) {
      BallisticPathPtr_t bp = BallisticPath::create (robot_->device_, q1, q2,
						     length,
						     coefficients);
      return bp;
    }

    // ========================================================================
    
    // assumes unit direction
    std::vector<bool> BallisticInterpolation::setMaintainRotationConstraints(const fcl::Vec3f&) // direction)
    {
      std::vector<bool> res;
      for(std::size_t i =0; i <3; ++i)
        {
	  res.push_back(true);
        }
      return res;
    }

    std::vector<bool> BallisticInterpolation::setTranslationConstraints(const fcl::Vec3f&)// normal)
    {
      std::vector<bool> res;
      for(std::size_t i =0; i <3; ++i)
        {
	  res.push_back(true);
        }
      return res;
    }

    void BallisticInterpolation::LockJointRec(const std::string& limb, const model::JointPtr_t joint, core::ConfigProjectorPtr_t& projector )
    {
        if(joint->name() == limb){      
          hppDout(notice,"LockJointrec : name = "<<joint->name()); 
          return;
        }
        const core::Configuration_t& c = joint->robot()->currentConfiguration();
        core::size_type rankInConfiguration (joint->rankInConfiguration ());
        projector->add(core::LockedJoint::create(joint,c.segment(rankInConfiguration, joint->configSize())));
        for(std::size_t i=0; i< joint->numberChildJoints(); ++i)
        {
            LockJointRec(limb,joint->childJoint(i), projector);
        }
    }

   /* bool BallisticInterpolation::ComputeCollisionFreeConfiguration
    (State& current, core::CollisionValidationPtr_t validation,
     const hpp::rbprm::RbPrmLimbPtr_t& limb,
     model::ConfigurationOut_t configuration,
     const double robustnessTreshold, bool stability)
    {
      for(std::vector<sampling::Sample>::const_iterator cit = limb->sampleContainer_.samples_.begin();
	  cit != limb->sampleContainer_.samples_.end(); ++cit)
        {
	  sampling::Load(*cit, configuration);
	  hpp::core::ValidationReportPtr_t valRep (new hpp::core::CollisionValidationReport);
	  if(validation->validate(configuration, valRep) && (!stability || stability::IsStable(robot_,current) >=robustnessTreshold))
            {
	      current.configuration_ = configuration;
	      return true;
            }
        }
      return false;
    }*/

//    State BallisticInterpolation::MaintainPreviousContacts
//    (const State& previous,
//     std::map<std::string,core::CollisionValidationPtr_t>& limbValidations,
//     model::ConfigurationIn_t configuration, bool& contactMaintained,
//     bool& multipleBreaks, std::vector <RbPrmLimbPtr_t>& successLimbs, const double robustnessTreshold)
//    {
//      hppDout (info, "configuration" <<displayConfig(configuration));
//      multipleBreaks = false;
//      contactMaintained = true;
//      std::vector<std::string> brokenContacts;
//      // iterate over every existing contact and try to maintain them
//      State current;
//      current.configuration_ = configuration;
//      model::Configuration_t config = configuration;
//      core::ConfigurationIn_t save = robot_->device_->currentConfiguration();
//      // iterate over contact filo list
//      std::queue<std::string> previousStack = previous.contactOrder_;
//      while(!previousStack.empty())
//        {
//	  const std::string name = previousStack.front();
//	  previousStack.pop();
//	  const RbPrmLimbPtr_t limb = robot_->GetLimbs().at(name);
//	  // try to maintain contact
//	  const fcl::Vec3f& ppos = previous.contactPositions_.at(name);
//	  const fcl::Vec3f& normal = previous.contactNormals_.at(name);
//    hppDout(info,"Maintain contact, previous contact position : "<<ppos);
//    hppDout(info,"Maintain contact, previous contact normal : "<<normal);
    
//	  core::ConfigProjectorPtr_t proj = core::ConfigProjector::create(robot_->device_,"proj", 5*1e-2, 80);
//	  //BallisticInterpolation::LockJointRec(limb->limb_->name(), robot_->device_->rootJoint(), proj);
    
//    core::size_type rankInConfiguration (robot_->device_->rootJoint()->rankInConfiguration());
//    proj->add(core::LockedJoint::create(robot_->device_->rootJoint(),robot_->device_->currentConfiguration().segment(rankInConfiguration, robot_->device_->rootJoint()->configSize())));
    
//	  const fcl::Vec3f z = limb->effector_->currentTransformation().getRotation() * limb->normal_;
//	  const fcl::Matrix3f& rotation = previous.contactRotation_.at(name);
//	  //proj->add(core::NumericalConstraint::create (constraints::Position::create("",robot_->device_, limb->effector_,fcl::Vec3f(0,0,0), ppos)));
//	  fcl::Transform3f localFrame, globalFrame;
//	  localFrame.setTranslation(ppos);
//	  proj->add(core::NumericalConstraint::create (constraints::Position::create("",robot_->device_, limb->effector_,globalFrame, localFrame, setTranslationConstraints(normal))));
//	  if(limb->contactType_ == hpp::rbprm::_6_DOF)
//            {
//	      proj->add(core::NumericalConstraint::create (constraints::Orientation::create("", robot_->device_, limb->effector_,rotation,setMaintainRotationConstraints(z))));
//            }
//	  if(proj->apply(config))
//            {
//      hppDout(notice,"conf initial = "<<displayConfig(configuration - config ));
//	      hpp::core::ValidationReportPtr_t valRep (new hpp::core::CollisionValidationReport);
//	      if(/*limbValidations.at(name)->validate(config, valRep)*/ true)
//                {
//		  hppDout (info, "contact is valid");
//		  hppDout (info, "contact name= " << name);
//		  current.contacts_[name] = true;
//		  current.contactPositions_[name] = previous.contactPositions_.at(name);
//		  current.contactNormals_[name] = previous.contactNormals_.at(name);
//		  current.contactRotation_[name] = previous.contactRotation_.at(name);
//		  current.contactOrder_.push(name);
//		  current.configuration_ = config;
//		  successLimbs.push_back (limb);
//                }
//	      else
//                {
//		  hppDout (info, "contact is not valid");
//		  hppDout (info, "contact name= " << name);
//		  contactMaintained = false;
//          rbprm::ComputeCollisionFreeConfiguration(current,limbValidations.at(name),limb,current.configuration_,robustnessTreshold,false); // no more necessary since we will take the limb-config from an interpolation
//		  brokenContacts.push_back(name);
//		  hppDout (info, "brokenContacts name= " << name);
//                }
//            }
//	  else
//            {
//	      hppDout (info, "projection cannot be applied");
//	      hppDout (info, "contact name= " << name);
//	      contactMaintained = false;
//          rbprm::ComputeCollisionFreeConfiguration(current,limbValidations.at(name),limb,current.configuration_,robustnessTreshold,false); // no more necessary since we will take the limb-config from an interpolation
//	      brokenContacts.push_back(name);
//	      hppDout (info, "brokenContacts name= " << name);
//            }
//        }
//      hppDout (info, "current config" <<displayConfig(current.configuration_));
//      // reload previous configuration
//      robot_->device_->currentConfiguration(save);
//      hppDout (info, "brokenContacts.size()= " << brokenContacts.size());
//      hppDout (info, "successLimbs.size()= " << successLimbs.size());
//      if(brokenContacts.size() > 1)
//        {
//	  contactMaintained = false;
//	  multipleBreaks = true;
//        }
//      return current;
//      }
  } // model
} //hpp
