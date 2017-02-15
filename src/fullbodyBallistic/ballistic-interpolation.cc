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
#include <hpp/model/body.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/core/basic-configuration-shooter.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/core/config-validations.hh>
#include <hpp/core/locked-joint.hh>
#include <hpp/core/path.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/rbprm/rbprm-state.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-path.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-interpolation.hh>
#include <hpp/rbprm/stability/stability.hh>
#include <hpp/rbprm/ik-solver.hh>
#include <hpp/rbprm/interpolation/limb-rrt-helper.hh>

#include <hpp/rbprm/interpolation/com-rrt.hh>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::value_type;
    using core::vector_t;
    using core::Configuration_t;

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
          lastRootIndex_ = lit->second->limb_->rankInConfiguration();
        }
      }
      hppDout(notice,"Last root index = "<<lastRootIndex_);
    }

    // ========================================================================
    
    fcl::Vec3f BallisticInterpolation::computeDir (const vector_t V0,
						   const vector_t Vimp) {
      fcl::Vec3f dir;
      for (int i = 0; i < 3; i++) {
	//dir [i] = V0 [i] - Vimp [i];
	dir [i] = Vimp [i] - V0 [i]; // inversion for EFORT_normal
      }
      dir.normalize ();
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
    
    Configuration_t replaceFixedLimb(Configuration_t trunk, Configuration_t fixedLimb){
      Configuration_t result(trunk);
      for(size_t i = 0 ; i < trunk.size() ; i++){
        if(fixedLimb[i] != 0){
          result[i] = fixedLimb[i];
          hppDout(notice,"replace limb "<<i<<"result = "<<result[i]);
        }
      }
      return result;
    }

    State BallisticInterpolation::computeOffsetContactConfig
    (const BallisticPathPtr_t bp, const State& previousState,
     T_StateFrame &stateFrames, const value_type u_offset,
     const bool increase_u_offset, value_type& length,
     const std::size_t maxIter) {
      Configuration_t q_trunk_offset, q_contact_offset(previousState.configuration_), q_interp, result;
      q_trunk_offset = previousState.configuration_; 
      hppDout(notice,"q_trunk_offset = "<<displayConfig(q_trunk_offset));
      core::DevicePtr_t robot = robot_->device_;
      bool success, contact_OK = true, multipleBreaks, contactMaintained;
      std::vector<std::string> contactingLimbs;
      Configuration_t lastfixedLimb (robot->configSize());
      for(int i = 0 ; i < lastfixedLimb.size() ; i ++){lastfixedLimb[i] = 0;}
      Configuration_t fixedLimb(lastfixedLimb);
      
      std::size_t iteration = 0;
      fcl::Vec3f dir;
      State state,lastState, lastValidState;
      lastState = previousState;
      value_type currentLength, previousLength;
      value_type u;
      T_StateFrame reverseFrame;
      core::ValidationReportPtr_t validationReport;
      const bool allowFailureMaintainContacts = true;
      const std::size_t nbTries6DOF = 5; // allow to try 5 times the 6DOF constraint before trying 3DOF instead

      hppDout (info, "previousState config for computeOffsetContact= " << displayConfig(previousState.configuration_));
      for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
        if(lastState.contacts_[lit->first]){ // limb is in contact
	  contactingLimbs.push_back(lit->first);
	  hppDout(notice,"contacting limbs : "<<lit->first);
        } else
	  hppDout(info,"NOT contacting limbs: " << lit->first);
      }
      
      if(increase_u_offset){ // forward progression
        u = u_offset;
        currentLength = 0;
      }else{ // backward progression
        u = - u_offset;
        currentLength = bp->length();
      }
      previousLength = currentLength;
     
      q_interp = (*bp) (currentLength, success);
      
      if (u_offset == 0) {
        hppDout (info, "no cushion effect asked, return interpolated config");
        state.configuration_ = q_interp;
        return state;
      }
      
      while (contact_OK && iteration < maxIter && ((((currentLength<bp->length())/3.)&&increase_u_offset) || ((currentLength > (bp->length()*2./3.))&&(!increase_u_offset)))) {
	currentLength += u; // u can be negative
        iteration++;
        hppDout (info, "iteration= " << iteration);
        hppDout (info, "currentLength = " << currentLength);
        q_trunk_offset = (*bp) (currentLength, success);
	hppDout (info, "q_trunk_offset= " << displayConfig(q_trunk_offset));
        dir = bp->evaluateVelocity (currentLength);
	if (increase_u_offset) {
	  // for takeoff, inverse dir  AS DONE IN CORBA (V0) !!
	  dir = - dir;
	}
	hppDout (info, "dir computeOffsetContactConfig= " << dir);
        state = rbprm::ComputeContacts(lastState,robot_,q_trunk_offset, affMap_, affFilters_, dir,contactMaintained,multipleBreaks,allowFailureMaintainContacts,0., nbTries6DOF);
	//state = rbprm::ComputeContacts(previousState,robot_,q_trunk_offset, affMap_, affFilters_, dir,contactMaintained,multipleBreaks,allowFailureMaintainContacts,0.); // ALWAYS USE FIRST STATE AS REF
	hppDout (info, "config after maintainContacts= " << displayConfig(state.configuration_));
	hppDout (info, "## contactMaintained = " << contactMaintained);

	for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
	  if(state.contacts_[lit->first]){ // limb is in contact
	    hppDout(notice,"contacting limbs: "<<lit->first);
	  } else
	    hppDout(info,"NOT contacting limbs: " << lit->first);
	}

        if(contactMaintained){
          contact_OK = true;
        }else {
          hppDout(notice,"contact released, set contact_OK to false");
          contact_OK = false; //NO! other contacts may be maintained
          for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
	    value_type diffNorm = 1;
            if(state.contacts_[lit->first]){ // get limbs still in contact
	      const fcl::Vec3f& scp = state.contactPositions_[lit->first];
	      const fcl::Vec3f& lscp = lastState.contactPositions_[lit->first];
	      diffNorm = sqrt(pow(scp[0] - lscp[0],2) + pow(scp[1] - lscp[1],2) + pow(scp[2] - lscp[2],2));
	      hppDout (info, "state.contactPositions_[lit->first]= " << scp);
	      hppDout (info, "lastState.contactPositions_[lit->first]= " << lscp);
	      hppDout (info, "limb: " << lit->first << " diffNorm= " << diffNorm);
	      if (diffNorm <= 1e-3) {
		// limb contact is OK -> at least one contact OK for this state
		hppDout (info, "set contact_OK to true");
		contact_OK = true;
	      } else
		hppDout (info, "contact could not be maintained (another was created)");
	    }
	    else if (!state.contacts_[lit->first] || diffNorm > 1e-3) { // limb not in contact or contact is NOT OK
              if(lastState.contacts_[lit->first]) { // limb just loose the contact
                for(int i=robot_->GetLimbs().at(lit->first)->limb_->rankInConfiguration() ; i < robot_->GetLimbs().at(lit->first)->effector_->rankInConfiguration() ; i++){
                  fixedLimb[i] = lastState.configuration_[i];
                  hppDout(notice,"fixed limb "<<i<<" = "<<fixedLimb[i]);
                }
              }
            }
          }// for limbs
	}// if contactNotMaintained
	if(contact_OK) { // there is at least one limb with contact maintained
	  hppDout(notice,"another limb in contact, add intermediate state : "<<displayConfig(state.configuration_)); 
	  lastfixedLimb = fixedLimb;

	  //replace limbs that are accurate, trunkDOF with contactPose value, EC
	  state = replaceAccurateTrunkAndLimbs (state, q_trunk_offset, increase_u_offset, contactingLimbs, lastfixedLimb); // NO USE OF LASTSTATE ??

	  bool isValid = problem_->configValidations ()->validate (state.configuration_, validationReport);
	  hppDout (info, "isValid state= " << isValid);
	  hppDout (info, "currentLength= " << currentLength);
	  hppDout (info, "previousLength= " << previousLength);
	  hppDout (info, "fabs(currentLength - previousLength)= " << fabs(currentLength - previousLength));
	  bool stateHasImproved = fabs(currentLength - previousLength) >= u_offset - 1e-6; // to prevent adding the (same) state two consecutive times
	  hppDout (info, "stateHasImproved= " << stateHasImproved);
            
	  if(stateHasImproved && isValid) {
	    if (increase_u_offset) {
	      hppDout (info, "takeoff contact state is pushed in stack");
	      hppDout (info, "config pushed= " << displayConfig(state.configuration_));
	      stateFrames.push_back(std::make_pair(currentLength,state)); // add intermediate state
	    } else {
	      hppDout (info, "landing contact state is pushed in stack");
	      hppDout (info, "config pushed= " << displayConfig(state.configuration_));
	      reverseFrame.push_back(std::make_pair(currentLength,state)); // add intermediate state
	    }
	    lastState = state;
	    previousLength = currentLength; // only updated if contact was OK
	  } else {
	    hppDout (info,"invalid or non-improving state NOT pushed in stack");
	    hppDout (info, "contact_OK set to false to break while");
	    contact_OK = false;
	  }
	} // if contactOK
        hppDout (info, "q_contact_offset= " << displayConfig (q_contact_offset));
      }//while
      hppDout (info, "-- contact broken or max number of iterations reached");
      hppDout (info, "-- or last computed state was invalid (not added)");

      // push reverseStateFrames into stateFrames in correct order
      for (int i = reverseFrame.size () - 1; i >= 0; i--) { // no std::size_t
	hppDout (info, "push reverse state length= " << reverseFrame[i].first);
	hppDout (info, "state= " << displayConfig(reverseFrame[i].second.configuration_));
	stateFrames.push_back(reverseFrame [i]);
      }

      hppDout (info, "lastState= " << displayConfig(lastState.configuration_));
      hppDout (info, "previousLength= " << previousLength);
      length = previousLength; // last valid length
      return lastState;
    }

    State BallisticInterpolation::computeTopExtendingPose 
    (const core::PathPtr_t path, const BallisticPathPtr_t bp,
     value_type& length) {
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
      hppDout (info, "x_theta_init" << x_theta_init);
      hppDout (info, "x_theta_end" << x_theta_end);
      hppDout (info, "x_theta_max" << x_theta_max);
      if (x_theta_max < x_theta_init || x_theta_max > x_theta_end) {
	hppDout (info, "q_top is not relevant");
	return state;
      }
      const value_type u_max = (x_theta_max - x_theta_init)
	/ (x_theta_end - x_theta_init); // in [0,1]
      const value_type z_x_theta_max = coefs (0)*x_theta_max*x_theta_max +
	coefs (1)*x_theta_max + coefs (2);
      core::Configuration_t q_top;
      const std::string robotName = robot_->device_->name ();
      hppDout (info, "robotName= " << robotName);
      hppDout (info, "z_x_theta_max= " << z_x_theta_max);
      const value_type AntThresh = 0.44;
      const value_type JumpermanThresh = 1.1;
      const value_type FrogThresh = 0.2;
      const value_type ArmlessSkeletonThresh = 0.6;
      const value_type SkeletonThresh = 0.6;
      const value_type LampThresh = 0.75;
      const bool parabTallEnough = (robotName.compare ("ant") == 0 &&  z_x_theta_max > AntThresh) || (robotName.compare ("spiderman") == 0 &&  z_x_theta_max > JumpermanThresh) || (robotName.compare ("frog") == 0 &&  z_x_theta_max > FrogThresh) || (robotName.compare ("armlessSkeleton") == 0 &&  z_x_theta_max > ArmlessSkeletonThresh) || (robotName.compare ("skeleton") == 0 &&  z_x_theta_max > SkeletonThresh) || (robotName.compare ("lamp") == 0 &&  z_x_theta_max > LampThresh);
      value_type r = 0.5; // blending coefficient of extending key-frame
      const Configuration_t q_interp_top = (*bp) (u_max*pathLength, success);
      hppDout (info, "q_interp_top= " << displayConfig(q_interp_top));
      core::Configuration_t q;
      length = u_max*pathLength;
      if (parabTallEnough)
	r = 0.8;
      else
	r = 0.4;
      hppDout (info, "r of blending extending pose= " << r);
      
      if (extendingPose_.rows ()) {
	const core::Configuration_t q_trunk_max = (*path) (u_max*pathLength, success);
	hppDout (info, "q_trunk_max= " << displayConfig(q_trunk_max));
	q = fillConfiguration (q_trunk_max, extendingPose_); // now q is extending_ at the good top-trunk-configuration
	q_top = blendPoses (q, q_interp_top, r);
      } else {
	q_top = q_interp_top;
	hppDout (info, "no extending pose");
      }
      hppDout (info, "q_top= " << displayConfig(q_top));

      // If q_top is in collision, try random stuff to find a collision-free one
      core::ValidationReportPtr_t validationReport;
      if (!problem_->configValidations ()->validate (q_top, validationReport)) {
	const core::ConfigurationShooterPtr_t configurationShooter = core::BasicConfigurationShooter::create (robot_->device_);
	hppDout (info, "q_top is NOT valid");
	core::Configuration_t qtmp = q_top;
	std::size_t count = 0;
	rbprm::RbPrmLimbPtr_t limb;
	r = 0.2;
	do{ 
	  count ++;
	  // try to move the limb that is in collision...
	  limb = findLimbInCollision (validationReport);
	  if (!limb) { // weird collision, return empty state
	    hppDout (info, "no limb found in collision, return empty topState");
	    return state;
	  }
	  std::vector <RbPrmLimbPtr_t> limbsvec; // just to use function
	  limbsvec.push_back (limb);
	  core::Configuration_t qrand = *configurationShooter->shoot ();
	  if (q.rows ())
	    q_top = blendPoses (q, qrand, r);
	  qtmp = replaceLimbConfigsInFullConfig (qtmp, qrand, limbsvec);
	  q_top = qtmp;
	} while (!problem_->configValidations ()->validate (qtmp, validationReport) && count < 10000);
	hppDout (info, "valid q_top now: " << displayConfig(qtmp));
      }
      state.configuration_ = q_top;      
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
          for(int i = lit->second->limb_->rankInConfiguration() ; i < lit->second->effector_->rankInConfiguration() ; i++){
            q[i] = contactPose[i];
          }
        } 
      }
      // replace the trunkDOF with contactPose value : ( 7 because we suppose we always work with freeflyer as root ....)
      if(!useAllLimb){
	hppDout (info, "in compute contact-pose, modify also trunk DOF");
	hppDout (info, "minIndex= " << minIndex); // maybe be out of bounds
        for (size_t i = 7 ; i < minIndex; i++){
          q[i] = contactPose[i];
        }
      }
      return q;
    }
 

    Configuration_t BallisticInterpolation::blendPoses 
    (const Configuration_t q1, const Configuration_t q2, const value_type r)
    {
      hppDout (info, "blendPoses q1= " << displayConfig(q1));
      hppDout (info, "blendPoses q2= " << displayConfig(q2));
      Configuration_t result = r*q1 + (1-r)*q2;
      return result;
    }

    rbprm::RbPrmLimbPtr_t BallisticInterpolation::findLimbInCollision
    (const core::ValidationReportPtr_t validationReport) {
      rbprm::RbPrmLimbPtr_t limb, emptyLimb;
      const T_Limb& robotLimbs = robot_->GetLimbs();
      const core::DevicePtr_t device = robot_->device_;
      core::CollisionValidationReportPtr_t colValRep = boost::dynamic_pointer_cast<core::CollisionValidationReport> (validationReport);
	  std::string collisionName = colValRep->object1->name ();
	  hppDout (info, "body in collision= " << collisionName);
	  const std::size_t sz = collisionName.size ();
	  collisionName.resize (sz - 2); // remove "_0" part
	  hppDout (info, "body in collision (resized)= " << collisionName);
	  for(rbprm::T_Limb::const_iterator cit = robotLimbs.begin(); cit != robotLimbs.end(); cit++) {
	    hppDout (info, "limb= " << cit->first);
	    limb = cit->second;
	    model::JointPtr_t effectorClone = device->getJointByName(limb->effector_->name ());
	    const std::size_t startRank (limb->limb_->rankInConfiguration());
	    const std::size_t limbLength (computeLimbLength (limb->limb_, effectorClone));
	    for (std::size_t i = 0; i < limbLength; i++) {
	      std::size_t rank = startRank + i;
	      const std::string bodyName = device->getJointAtConfigRank (rank)->linkedBody ()->name ();
	      hppDout (info, "bodyName= " << bodyName);
	      if (collisionName.compare (bodyName) == 0) { // limb found
		hppDout (info, "found limb in collision");
		return limb;
	      }
	    }
	  }
	  return emptyLimb;
    }

    State BallisticInterpolation::replaceAccurateTrunkAndLimbs
    (const State& refState, const core::Configuration_t refTrunk,
     const bool increase_u_offset,
     const std::vector<std::string> contactingLimbs,
     const core::Configuration_t lastfixedLimb) {
      State state = refState;
      const std::size_t ecsSize = 
	robot_->device_->extraConfigSpace ().dimension ();
      const std::size_t ecIndex = refState.configuration_.size () - ecsSize;

      state.configuration_ = replaceFixedLimb(state.configuration_,
					      lastfixedLimb);

      // take limb-configs from contact if contact succeeded, 
      // from interp otherwise.
      // replace the limbs not used for contact with their configuration in flexionPose_
      size_t minIndex = robot_->device_->configSize();
      core::ValidationReportPtr_t report;
      core::Configuration_t qtmp (state.configuration_);
      for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
        hppDout(notice,"LIST OF LIMBS : "<< lit->first << " contact = "<<state.contacts_[lit->first]);
        if(lit->second->limb_->rankInConfiguration() < minIndex){
          minIndex =lit->second->limb_->rankInConfiguration();
        }
        if(!state.contacts_[lit->first]){ // limb was in contact in refState
          if(std::find(contactingLimbs.begin(),contactingLimbs.end(),lit->first) == contactingLimbs.end()) {
            hppDout(notice," Not in contact, index config : "<<lit->second->limb_->rankInConfiguration()<<" -> "<<lit->second->effector_->rankInConfiguration());
            for(int i = lit->second->limb_->rankInConfiguration() ; i < lit->second->effector_->rankInConfiguration() ; i++){
	      if(increase_u_offset && takeoffContactPose_.size () != 0)
		qtmp[i] = takeoffContactPose_[i];
	      else if (!increase_u_offset && landingContactPose_.size () != 0)
		qtmp[i] = landingContactPose_[i];
            }
            if(robot_->getLimbcollisionValidations().at(lit->first)->validate(qtmp,report)) {
	      hppDout (info, "state with contactPose LIMB DOF is valid, DOF kept");
              state.configuration_ = qtmp;
	    }
          }
        }
      }//for

      // replace the trunkDOF with contactPose value : (we suppose we always work with freeflyer as root ....)
      // (mylene)NO! this is wrong because trunk was planned before
      // trunk DOF should be the value on the parab
      hppDout (info, "using takeoffContactPose and landing for trunk DOF");
      hppDout (info, "minIndex= " << minIndex);
      hppDout (info, "trunkConfigSize_= " << trunkConfigSize_); // no EC
      for (size_t i = 7 ; i < trunkConfigSize_ ; i++) {
	if(increase_u_offset && takeoffContactPose_.size () != 0)
	  qtmp[i] = takeoffContactPose_[i];
	else if (!increase_u_offset && landingContactPose_.size () != 0)
	  qtmp[i] = landingContactPose_[i];
      }
      if(robot_->getCollisionValidation()->validate(qtmp,report)){
	hppDout (info, "state with contactPose TRUNK DOF is valid, DOF kept");
        state.configuration_ = qtmp;
      }

      // set EC of state from trunk config
      if (ecsSize > 0) {
	for (std::size_t k = 0; k < ecsSize; k++) {
	  state.configuration_ [ecIndex + k] = refTrunk [trunkConfigSize_ + k];
	}
      }
      return state;
    }

    State BallisticInterpolation::createStateFromCone
    (const State previousState, library::ContactCones cone) {
      State state = previousState;
      state.nbContacts = cone.coneNumber_;
      for (std::size_t i = 0; i < cone.coneNumber_; i++) {
	for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
	  const hpp::rbprm::RbPrmLimbPtr_t& limb = lit->second;
	  hppDout (info, "limb: " << lit->first << ", cone ROMname: " << cone.ROMnames_ [i]);
	  if (cone.ROMnames_ [i].compare (lit->first) == 0) {
	    // create contact
	    state.contacts_ [lit->first] = true;
	    state.contactNormals_ [lit->first] = cone.directions_ [i];
	    state.contactPositions_ [lit->first] = cone.positions_ [i];
	    // ?? current state of robot important ?
	    robot_->device_->currentConfiguration (previousState.configuration_);
	    const fcl::Vec3f z = limb->effector_->currentTransformation().getRotation() * limb->normal_;
	    const fcl::Matrix3f alignRotation = tools::GetRotationMatrix(z,cone.directions_ [i]);
	    state.contactRotation_ [lit->first] = alignRotation * limb->effector_->currentTransformation().getRotation();
	  }//if
	}//for limbs
      }//for cones
      return state;
    }

    // ========================================================================

     T_PathVectorBP BallisticInterpolation::InterpolateFullPath
     (const core::value_type u_offset, T_StateFrame* stateFramesRef) {
      if(!path_) throw std::runtime_error ("Cannot interpolate; not path given to interpolator ");
      //Configuration_t qStart = computeFlexionContactPose(start_);
      //Configuration_t qEnd = computeFlexionContactPose(end_);
      Configuration_t qStart = start_.configuration_;
      Configuration_t qEnd = end_.configuration_;
      core::DevicePtr_t robot = robot_->device_;
      T_StateFrame stateFrames;
      Configuration_t q2(robot->configSize ()),
	q1contact (robot->configSize ()), q2contact (robot->configSize ());
      BallisticPathPtr_t bp;
      hppDout (info, "u_offset= " << u_offset);
      const std::size_t subPathNumber = path_->numberPaths ();
      hppDout (info, "number of sub-paths: " << subPathNumber);
      std::vector<core::PathVectorPtr_t> newPathVector;
      robot_->noStability_ = true; // disable stability for waypoints
      vector_t V0 (3), Vimp (3); fcl::Vec3f dir;
      core::PathPtr_t subpath = path_->pathAtRank (0);
      State state1, state2;
      State contactState1,contactState2,stateTop;
      core::value_type lengthTop,lengthTakeoff,lengthLanding;//,lengthLanding6DOF,lengthTakeoff6DOF;
      core::ValidationReportPtr_t validationReport;
      T_PathVectorBP resultPairVector;
      trunkConfigSize_ = path_->pathAtRank (0)->initial ().size () - 4;

      for (std::size_t i = 0; i < subPathNumber - 1; i++) {
	stateFrames.clear();      
	hppDout (info, "B-interp on path nb: " << i);
	core::PathPtr_t subpath_next = path_->pathAtRank (i+1);
	ParabolaPathPtr_t pp = 
	  boost::dynamic_pointer_cast<ParabolaPath>(subpath);
	ParabolaPathPtr_t pp_next = 
	  boost::dynamic_pointer_cast<ParabolaPath>(subpath_next);
	const value_type pathLength = subpath->length ();
	core::PathVectorPtr_t pathLimb;
	
	if (i == 0) { // keep qStart config which already has contacts
	  hppDout (info, "keep start config");
	  q1contact = qStart;
	  state1 = start_;
	}
	else {
	  q1contact = q2contact;
	  state1 = state2;
	}

	  const std::string robotName = robot_->device_->name ();
	if (robotName.compare ("spiderman") == 0 && takeoffContactPose_.size () != 0)
	  q2 = fillConfiguration (subpath->end (), takeoffContactPose_);
	else if (flexionPose_.size () != 0)
	  q2 = fillConfiguration (subpath->end (), flexionPose_);
	else
	  q2 = fillConfiguration (subpath->end (), robot->configSize ());
	hppDout (info, "q2: " << displayConfig(q2));
	Configuration_t q2next = subpath_next-> end (); // for theta
	hppDout (info, "q2next: " << displayConfig(q2next));
	V0 = pp_next->V0_; // V0_i+1
	Vimp = pp->Vimp_; // Vimp_i
	dir = computeDir (V0, Vimp);
	hppDout (info, "dir = " << dir);
	hppDout (info, "(Vimp-V0)= " << -computeDir (V0, Vimp));
	robot_->V0dir_ = V0;
	robot_->Vfdir_ = Vimp;
	robot_->thetaBefore_ = atan2 (q2[1] - q1contact[1], q2[0] - q1contact[0]);
	robot_->thetaAfter_ = atan2 (q2next[1] - q2[1], q2next[0] - q2[0]);
	hppDout (info, "robot_->thetaBefore_= " << robot_->thetaBefore_);
	hppDout (info, "robot_->thetaAfter_= " << robot_->thetaAfter_);

	state2 = ComputeContacts(robot_, q2, affMap_, affFilters_, dir);
	hppDout (info, "state2 config= " << displayConfig(state2.configuration_));

	// try to create planning contacts
	// NOT accurate because of trunk-re-orientation after planning
	/*bool contactMaintained, multipleBreaks;
	library::ContactCones cones = pp->contactConesImp_;
	State state2repos = createStateFromCone (state2, cones);
	Configuration_t state2repos_trunkConf = state2repos.configuration_;
	state2repos_trunkConf.resize (trunkConfigSize_);
	const Configuration_t state2repos_trunkConfLong = fillConfiguration(state2repos_trunkConf, robot->configSize ());
	hppDout (info, "state2repos_trunkConfLong= " << displayConfig(state2repos_trunkConfLong));
	State stateTest = rbprm::ComputeContacts(state2repos,robot_, state2repos_trunkConfLong, affMap_, affFilters_, dir,contactMaintained,multipleBreaks,false,0., 5);
	hppDout (info, "stateTest = "<<displayConfig(stateTest.configuration_));
	for (std::size_t i = 0; i < cones.coneNumber_; i++)
	  hppDout (info, "planningContactVerif limb= "  << cones.ROMnames_ [i] << " position= " << cones.positions_ [i] << " direction= " << cones.directions_ [i]);
	*/

	q2contact = computeFlexionContactPose (state2);
	hppDout (info, "q2contact = " << displayConfig(q2contact));
	if (problem_->configValidations ()->validate (q2contact, validationReport))
	  state2.configuration_ = q2contact; // sometimes, state2.configuration_ is "a little" in collision whereas q2contact is not
	else {
	  hppDout (info, "q2contact is abnormally in collision, retry state2");
	  if (!problem_->configValidations ()->validate (state2.configuration_, validationReport))
	    throw std::runtime_error ("q2contact and state2 are abnormally in collision");
	  else
	    q2contact = state2.configuration_;
	}
	// reset params in case of computeContact request somewhere else
	robot_->thetaBefore_ = NULL;
	robot_->thetaAfter_ = NULL;
   
	hppDout (info, "q2contact= " << displayConfig(q2contact));
	hppDout (info, "q1contact= " << displayConfig(q1contact));

	bp = BallisticPath::create (robot_->device_, q1contact, q2contact, pathLength, subpath->coefficients ());
	bp->lastRootIndex(lastRootIndex_);
	stateTop = computeTopExtendingPose (subpath, bp,lengthTop);

	if (stateTop.configuration_.rows())
	  hppDout (info, "topState validity: " << problem_->configValidations ()->validate (stateTop.configuration_, validationReport));
	const bool state1validity = problem_->configValidations ()->validate (state1.configuration_, validationReport);
	const bool state2validity = problem_->configValidations ()->validate (state2.configuration_, validationReport);
	hppDout (info, "state1 validity: " << state1validity);
	hppDout (info, "state2 validity: " << state2validity);

	stateFrames.clear();
	if (state1validity)
	  stateFrames.push_back(std::make_pair(0,state1));
	contactState1 = computeOffsetContactConfig (bp, state1,stateFrames, u_offset, true,lengthTakeoff);
	if (stateTop.configuration_.rows())
	  stateFrames.push_back(std::make_pair(lengthTop,stateTop));
	contactState2 = computeOffsetContactConfig (bp, state2,stateFrames, u_offset, false,lengthLanding);
	if (state2validity)
	  stateFrames.push_back(std::make_pair(bp->length(),state2));
    
	hppDout(notice, "position initial state frame  = "<<displayConfig(state1.configuration_));
	hppDout(notice, "position initial Contact state frame  = "<<displayConfig(contactState1.configuration_));
	hppDout(notice, "position top frame  = "<<displayConfig(stateTop.configuration_));
	hppDout(notice, "position final contact frame  = "<<displayConfig(contactState2.configuration_));
	hppDout(notice, "position final state frame  = "<<displayConfig(state2.configuration_));
	hppDout(notice, "TIME initial state frame  = "<<0);
	//hppDout(notice, "TIME initial Contact transition state frame  = "<<lengthTakeoff6DOF);
	hppDout(notice, "TIME initial Contact state frame  = "<<lengthTakeoff);
	hppDout(notice, "TIME top frame  = "<<lengthTop);
	hppDout(notice, "TIME final contact frame  = "<<lengthLanding);
	//hppDout(notice, "TIME final contact transition frame  = "<<bp->length() - lengthLanding6DOF);
	hppDout(notice, "TIME final state frame  = "<<bp->length());
	hppDout(notice,"test last root index interpolate = "<<bp->lastRootIndex());
    
	hppDout (info, "before interpolateStates in fullpath");
	hppDout (info, "number of states= " << stateFrames.size ());
	pathLimb = rbprm::interpolation::interpolateStatesinPathVector
	  (robot_,problem_,bp,stateFrames.begin(),stateFrames.end()-1,2);
	hppDout (info, "after interpolateStates in fullpath");
	resultPairVector.push_back(std::make_pair(pathLimb,bp));
	subpath = subpath_next;
	pp = boost::dynamic_pointer_cast<ParabolaPath>(subpath);
	
	if (i == subPathNumber - 2) { // subpath_next = final subpath
	  q1contact = q2contact;
	  state1 = state2;
	  q2contact = qEnd;
	  state2 = end_;

	  // Compute top-keyframe and adapt ballistic-path with it
	  bp = BallisticPath::create (robot_->device_, q1contact, q2contact, subpath->length (), subpath->coefficients ());
	  bp->lastRootIndex(lastRootIndex_);
	  stateTop = computeTopExtendingPose (subpath, bp,lengthTop);

	  if (stateTop.configuration_.rows())
	    hppDout (info, "topState validity: " << problem_->configValidations ()->validate (stateTop.configuration_, validationReport));
	  const bool state1validity = problem_->configValidations ()->validate (state1.configuration_, validationReport);
	  const bool state2validity = problem_->configValidations ()->validate (state2.configuration_, validationReport);
	  hppDout (info, "state1 validity: " << state1validity);
	  hppDout (info, "state2 validity: " << state2validity);
    
	  stateFrames.clear();
	  if (state1validity)
	    stateFrames.push_back(std::make_pair(0,state1));
	  contactState1 = computeOffsetContactConfig (bp, state1,stateFrames, u_offset, true,lengthTakeoff);
	  if (stateTop.configuration_.rows())
	    stateFrames.push_back(std::make_pair(lengthTop,stateTop));
	  contactState2 = computeOffsetContactConfig (bp, state2,stateFrames, u_offset, false,lengthLanding);
	  if (state2validity)
	    stateFrames.push_back(std::make_pair(bp->length(),state2));
      
	  hppDout(notice, "position initial state frame  = "<<displayConfig(state1.configuration_));
	  hppDout(notice, "position initial Contact state frame  = "<<displayConfig(contactState1.configuration_));
	  hppDout(notice, "position top frame  = "<<displayConfig(stateTop.configuration_));
	  hppDout(notice, "position final contact frame  = "<<displayConfig(contactState2.configuration_));
	  hppDout(notice, "position final state frame  = "<<displayConfig(state2.configuration_));
	  hppDout(notice, "TIME initial state frame  = "<<0);
	  //hppDout(notice, "TIME initial Contact transition state frame  = "<<lengthTakeoff6DOF);
	  hppDout(notice, "TIME initial Contact state frame  = "<<lengthTakeoff);
	  hppDout(notice, "TIME top frame  = "<<lengthTop);
	  hppDout(notice, "TIME final contact frame  = "<<lengthLanding);
	  //hppDout(notice, "TIME final contact transition frame  = "<<bp->length() - lengthLanding6DOF);
	  hppDout(notice, "TIME final state frame  = "<<bp->length());
	  hppDout(notice,"test last root index interpolate = "<<bp->lastRootIndex());
    
	  *stateFramesRef = stateFrames;
	  hppDout (info, "number of states= " << stateFrames.size ());

	  pathLimb = rbprm::interpolation::interpolateStatesinPathVector (robot_,problem_,bp,stateFrames.begin(),stateFrames.end()-1,2); // interpolateStates in limb-rrt-helper
	  hppDout (info, "after interpolateStates final subpath");
	  resultPairVector.push_back(std::make_pair(pathLimb,bp));
      
	}//if final subpath
      }// for subpaths
      return resultPairVector;
    }

    T_PathVectorBP
    BallisticInterpolation::InterpolateDirectPath
    (const core::value_type u_offset, T_StateFrame* stateFramesRef) {
      if(!path_) throw std::runtime_error ("Cannot interpolate; not path given to interpolator ");
      hppDout (info, "direct B-interpolation");
      //Configuration_t qStart = computeFlexionContactPose (start_);
      //Configuration_t qEnd = computeFlexionContactPose (end_);
      
      hppDout (info, "check contacting limbs");
      State startCopy = start_;
      for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
        if(start_.contacts_[lit->first]){ // limb is in contact
	  hppDout(notice,"contacting limbs: "<<lit->first);
        } else
	  hppDout(info,"NOT contacting limbs: " << lit->first);
      }

      Configuration_t qStart = start_.configuration_;
      Configuration_t qEnd = end_.configuration_;
      core::DevicePtr_t robot = robot_->device_;
      BallisticPathPtr_t bp;
      T_PathVectorBP resultPairVector;
      robot_->noStability_ = true; // disable stability for waypoints
      State contactState1,contactState2,stateTop;
      T_StateFrame stateFrames;
      value_type lengthTop,lengthTakeoff,lengthLanding,lengthLanding6DOF,lengthTakeoff6DOF;
      const core::PathPtr_t path = path_->pathAtRank (0);
      const value_type pathLength = path->length ();
      const vector_t pathCoefs = path->coefficients ();
      core::ValidationReportPtr_t validationReport;
      trunkConfigSize_ = path->initial ().size () - 4;

      bp = BallisticPath::create (robot_->device_, qStart, qEnd, pathLength, pathCoefs);
      stateTop = computeTopExtendingPose (path, bp,lengthTop);
      bp->lastRootIndex(lastRootIndex_);
      hppDout(notice, "full length = "<<bp->length());
      hppDout(notice,"length top = "<<lengthTop);

      if (stateTop.configuration_.rows())
	hppDout (info, "topState validity: " << problem_->configValidations ()->validate (stateTop.configuration_, validationReport));
      hppDout (info, "start_ validity: " << problem_->configValidations ()->validate (start_.configuration_, validationReport));
      hppDout (info, "end_ validity: " << problem_->configValidations ()->validate (end_.configuration_, validationReport));

      stateFrames.clear();
      stateFrames.push_back(std::make_pair(0,start_));
      hppDout (info, "number of states (1)= " << stateFrames.size ());
      contactState1 = computeOffsetContactConfig (bp, start_,stateFrames, u_offset, true,lengthTakeoff);
      hppDout (info, "number of states (1+offsetTakeoff)= " << stateFrames.size ());
      if (stateTop.configuration_.rows())
	stateFrames.push_back(std::make_pair(lengthTop,stateTop));
      hppDout (info, "number of states (2+offsetTakeoff)= " << stateFrames.size ());
      //contactState2 = computeOffsetContactConfig (bp, end_,stateFrames, u_offset, false,lengthLanding); // !!! DEBUG ONLY
      hppDout (info, "number of states (2+offsets)= " << stateFrames.size ());
      stateFrames.push_back(std::make_pair(bp->length(),end_));
      hppDout (info, "number of states (3+offsets)= " << stateFrames.size ());

      hppDout(notice, "position initial state frame  = "<<displayConfig(start_.configuration_));
      hppDout(notice, "position initial Contact state frame  = "<<displayConfig(contactState1.configuration_));
      hppDout(notice, "position top frame  = "<<displayConfig(stateTop.configuration_));
      hppDout(notice, "position final contact frame  = "<<displayConfig(contactState2.configuration_));
      hppDout(notice, "position final state frame  = "<<displayConfig(end_.configuration_));
      hppDout(notice, "TIME initial state frame  = "<<0);
      //hppDout(notice, "TIME initial Contact transition state frame  = "<<lengthTakeoff6DOF);
      hppDout(notice, "TIME initial Contact state frame  = "<<lengthTakeoff);
      hppDout(notice, "TIME top frame  = "<<lengthTop);
      hppDout(notice, "TIME final contact frame  = "<<lengthLanding);
      //hppDout(notice, "TIME final contact transition frame  = "<<bp->length() - lengthLanding6DOF);
      hppDout(notice, "TIME final state frame  = "<<bp->length());
      hppDout(notice,"test last root index interpolate = "<<bp->lastRootIndex());

      for( rbprm::T_Limb::const_iterator lit = robot_->GetLimbs().begin();lit != robot_->GetLimbs().end(); ++lit){
        hppDout(notice,"LIST OF LIMBS END : "<< lit->first << "contact = "<<((end_.contacts_.find(lit->first) != end_.contacts_.end()) && end_.contacts_.at(lit->first)));
      }

      *stateFramesRef = stateFrames;
      hppDout (info, "number of states= " << stateFrames.size ());

      core::PathVectorPtr_t pathLimb = interpolation::interpolateStatesinPathVector (robot_,problem_,bp,stateFrames.begin(),stateFrames.end()-1,2); // limb-rrt-helper
      hppDout (info, "after interpolateStates direct path");
      resultPairVector.push_back(std::make_pair(pathLimb,bp));
      return resultPairVector;
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

    core::Configuration_t BallisticInterpolation::computeFlexionContactPose (const State& state) {
      const std::string robotName = robot_->device_->name ();
      if (robotName.compare ("spiderman") == 0) {
	const std::size_t nbContacts = state.nbContacts;
	hppDout (info, "nbContacts= " << nbContacts);
	/* DEBUG since cannot create contacts with hands !!
	  if (nbContacts <= 2 && flexionPose_.size() != 0) {
	  hppDout (info, "use flexion refConfig");
	  return rbprm::computeContactPose(state, flexionPose_, robot_);
	  }*/
	if (nbContacts > 2 && takeoffContactPose_.size() != 0) {
	  hppDout (info, "use takeoffContact refConfig");
	  return rbprm::computeContactPose(state, takeoffContactPose_, robot_);
	}
	return state.configuration_;
      }// if jumperman robot
      if (flexionPose_.size() != 0)
	return rbprm::computeContactPose(state,flexionPose_,robot_);
      else
	return state.configuration_;
    }

  } // rbprm
} //hpp
