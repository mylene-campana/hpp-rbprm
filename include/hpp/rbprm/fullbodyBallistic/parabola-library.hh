//
// Copyright (c) 2014 CNRS
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

#ifndef HPP_RBPRM_PARABOLA_LIBRARY_HH
# define HPP_RBPRM_PARABOLA_LIBRARY_HH

# include <sstream>
# include <hpp/util/debug.hh>
# include <Eigen/Geometry>
# include <hpp/model/configuration.hh>
# include <hpp/model/device.hh>
# include <hpp/core/problem.hh>
# include <hpp/core/config-validations.hh>
# include <hpp/rbprm/rbprm-validation.hh>
# include <hpp/rbprm/rbprm-validation-report.hh>
# include <polytope/stability_margin.h>
# include "utils/algorithms.h"

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using model::size_type;
    using core::value_type;
    using core::vector_t;


    /// (R2 R1xR2 R1) = A.(r2 r1xr2 r1) where A is the orientation matrix convertible in quaternion
    /// r are expressed in the local frame (character root), R in the global frame
    /// By convention, r2 = (1,0,0), r1 = (0,0,1) so (r2 r1xr2 r1) = Id
    /// n is the extra-configs, i.e. the obstacle normal
    /// Depending on the case, R1 = n so the biped stands on feet
    /// or R2 = -n so the quadruped is on all fours
    /// The remaining R_r is chosen according to theta (last extra-config)
    inline Eigen::Quaternion<value_type> quaternionFromNormalAndTheta
    (const fcl::Vec3f normal, const value_type theta,
     const std::string robotName = "") {
      const value_type nx = normal [0];
      const value_type ny = normal [1];
      const value_type nz = normal [2];

      const int sign = -(theta > M_PI/2) - (theta < -M_PI/2) +
	(int) ((theta > -M_PI/2) && (theta < M_PI/2));
      
      value_type R1x, R1y, R1z, R2x, R2y, R2z; // vector coordinates
      value_type a, b, c; // corresponding to R_r unknown coordinates

      // norm(R_r) = 1 ; R_r perp. to n and R_r following theta
      // See Matlab script or Latex doc for the different cases:
      if (theta != M_PI/2 && theta != -M_PI/2) {
	const value_type tantheta = tan(theta);
	if (nz != 0) {
	  a= sign*sqrt(nz*nz/((1+tantheta*tantheta)
			      *nz*nz+(nx+ny*tantheta)*(nx+ny*tantheta)));
	  b = a*tantheta;
	  c = -a*(nx+ny*tantheta)/nz;
	} else { // nz = 0 case
	  a = 0; b = 0; c = 1;
	}
      } else { // theta = +- pi/2
	a = 0;
	if (ny != 0) {
	  if (theta == M_PI/2)
	    c = sqrt(1/(1+(nz/ny)*(nz/ny)));
	  else // theta = -pi/2
	    c = sqrt(1/(1+(nz/ny)*(nz/ny)));
	  b = -nz/ny*c;
	} else { // ny = 0
	  if (nx == 1) { // no matter b and c)
	    b = 0.5;
	    c = 0.5;
	  } else { // nx = 1
	    b = 1;
	    c = 0;
	  }
	}
      }
      
      hppDout (info, "robotName= " << robotName);
      // quadRobot are 4-legged characters that usually stands on their two feet
      const bool quadRobot = robotName.compare("skeleton_trunk") == 0 || robotName.compare("spiderman_trunk") == 0;
      const bool almostVertical = nz > 0.707; // arbitrary threshold
      const bool xDir = quadRobot && !almostVertical;

      if (!xDir) { // R1 = n, R2 = R_r
	R1x = nx; R1y = ny; R1z = nz;
	R2x = a; R2y = b; R2z = c;
	hppDout (info, "zDir ON");
      }
      else { // R2 = -n, R1 = R_r
	R1x = a; R1y = b; R1z = c;
	R2x = -nx; R2y = -ny; R2z = -nz;
	hppDout (info, "zDir OFF, xDir ON");
      }

      // cross-product R1 x R2
      const value_type Rzx = R1z*R2x - R1x*R2z;
      const value_type Ryz = R1y*R2z - R1z*R2y;
      const value_type Rxy = R1x*R2y - R1y*R2x;
      
      Eigen::Matrix<value_type,3,3> A; //fcl::Matrix3f A;
      A (0,0) = R2x; A (0,1) = Ryz; A (0,2) = R1x;
      A (1,0) = R2y; A (1,1) = Rzx; A (1,2) = R1y;
      A (2,0) = R2z; A (2,1) = Rxy; A (2,2) = R1z;
      //hppDout (info, "A: " << A);
      //hppDout (info, "A.determinant (): " << A.determinant ());

      Eigen::Quaternion<value_type> quat (A);

      return quat;
    }

    // ---------------------------------------------------------------------

    /// Arrange robot orientation according to the surface normal direction
    /// So that the robot is "on" the surface, rotated
    inline core::Configuration_t setOrientation
    (const core::DevicePtr_t& robot, const core::Configuration_t& q) {
      core::Configuration_t qtest = q;
      const core::JointPtr_t jointSO3 = robot->getJointVector () [1];
      const size_type indexSO3 = jointSO3->rankInConfiguration ();
      const size_type index = robot->configSize ()
	- robot->extraConfigSpace ().dimension ();

      const value_type Nnorm = sqrt(pow(q [index],2) + pow(q [index+1],2) + pow(q [index+2],2));
      const value_type nx = q [index]/Nnorm;
      const value_type ny = q [index + 1]/Nnorm;
      const value_type nz = q [index + 2]/Nnorm;
      //const value_type theta = atan2 (2*qw*qz - 2*qx*qy, 1-2*qy*qy-2*qz*qz);
      const value_type theta = q [index + 3];

      fcl::Vec3f normal;
      normal [0] = nx; normal [1] = ny; normal [2] = nz;
      const std::string robotName = robot->name ();
      
      Eigen::Quaternion<value_type> quat = quaternionFromNormalAndTheta
	(normal, theta, robotName);

      const value_type qw = quat.w ();
      const value_type qx = quat.x ();
      const value_type qy = quat.y ();
      const value_type qz = quat.z ();
      const value_type magnitude = sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
      //hppDout (info, "quat: " << " " << qw << " " << qx << " " << qy << " " << qz);
      //hppDout (info, "quaternion magnitude= " << magnitude);
      
      // re-normalize (not needed but sometimes, loss of accuracy...)
      qtest [indexSO3] = qw / magnitude;
      qtest [indexSO3 + 1] = qx / magnitude;
      qtest [indexSO3 + 2] = qy / magnitude;
      qtest [indexSO3 + 3] = qz / magnitude;
      //hppDout (info, "qtest: " << displayConfig (qtest));
      return qtest;
    }

    // ---------------------------------------------------------------------

    /// Compute the GIWC for the given configuration, contactSize containts
    /// two values x y representing the length and width of the robot contact
    /// surfaces (e.g. robot feet).
    inline const polytope::ProjectedCone* computeConfigGIWC
    (core::ProblemPtr_t problem, const core::Configuration_t q,
     const core::vector_t contactSize) {
      hppDout(notice,"## compute GIWC");
      const core::DevicePtr_t robot (problem->robot ());
      const polytope::ProjectedCone* emptyGiwc = NULL;
      //const size_type index = robot->configSize () - robot->extraConfigSpace ().dimension ();
      // fill normal information in node (not sure if needed)
      //const core::vector_t normal (3);
      //normal[0] = q [index]; normal[1] = q [index+1]; normal[2] = q [index+2];

      core::ValidationReportPtr_t report;
      if (!(problem->configValidations())) {
	hppDout(error,"~~ No configs-validation in problem");
      }
      const bool isValid = problem->configValidations()->validate(q,report);
      if (!report) {
	hppDout(error,"~~ Report problem");
      }
      core::RbprmValidationReportPtr_t rbReport =
	boost::dynamic_pointer_cast<core::RbprmValidationReport> (report);
      
      if(!rbReport) {
	  hppDout(error,"~~ Validation Report cannot be cast");
	  return emptyGiwc;
	}
      if(!isValid) {
	  hppDout(warning,"~~ ComputeGIWC : config is not valid");
	  return emptyGiwc;
	}
      
      polytope::T_rotation_t rotContact(3*rbReport->ROMReports.size(),3);
      polytope::vector_t posContact(3*rbReport->ROMReports.size());
      
      // get the 2 object in contact for each ROM :
      hppDout(info,"~~ Number of roms in collision : "<<rbReport->ROMReports.size());
      size_t indexRom = 0 ;
      for(std::map<std::string,core::CollisionValidationReportPtr_t>::const_iterator it = rbReport->ROMReports.begin() ; it != rbReport->ROMReports.end() ; ++it)
	{
	  hppDout(info,"~~ for rom : "<<it->first);
	  core::CollisionObjectPtr_t obj1 = it->second->object1;
	  core::CollisionObjectPtr_t obj2 = it->second->object2;
	  hppDout(notice,"~~ collision between : "<<obj1->name() << " and "<<obj2->name());
	  fcl::CollisionResult result = it->second->result;
        
	  // get intersection between the two objects :
	  obj1->fcl();
	  geom::T_Point vertices1;
	  geom::BVHModelOBConst_Ptr_t model1 =  geom::GetModel(obj1->fcl());
	  hppDout(info,"vertices obj1 : "<<obj1->name()<< " ( "<<model1->num_vertices<<" ) ");
	  for(int i = 0 ; i < model1->num_vertices ; ++i)
	    {
	      vertices1.push_back(Eigen::Vector3d(model1->vertices[i][0], model1->vertices[i][1], model1->vertices[i][2]));
	    }
        
	  obj2->fcl();
	  geom::T_Point vertices2;
	  geom::BVHModelOBConst_Ptr_t model2 =  geom::GetModel(obj2->fcl());
	  hppDout(info,"vertices obj2 : "<<obj2->name()<< " ( "<<model2->num_vertices<<" ) ");
	  for(int i = 0 ; i < model2->num_vertices ; ++i)
	    {
	      vertices2.push_back(Eigen::Vector3d(model2->vertices[i][0], model2->vertices[i][1], model2->vertices[i][2]));
	    }
        
	  /// warning: plane normal harcoded to (0,0,1) here (still true ?)
	  geom::Point pn;
	  geom::T_Point hull = geom::intersectPolygonePlane(model1,model2,pn);

	  if(hull.size() == 0){
	    hppDout(error,"No intersection between rom and environnement");
	    return emptyGiwc;
	  }
        
	  // todo : compute center point of the hull
	  polytope::vector3_t normal,tangent0,tangent1;
	  geom::Point center = geom::center(hull.begin(),hull.end());
	  posContact.segment<3>(indexRom*3) = center;
	  polytope::rotation_t rot; 
	  normal [0] = -result.getContact(0).normal [0]; // of contact surface
	  normal [1] = -result.getContact(0).normal [1];
	  normal [2] = -result.getContact(0).normal [2];
	  hppDout(notice," !!! normal for GIWC : "<<normal.transpose ());
	  // compute tangent vector : 
	  tangent0 = normal.cross(polytope::vector3_t(1,0,0));
	  if(tangent0.dot(tangent0)<0.001)
	    tangent0 = normal.cross(polytope::vector3_t(0,1,0)); 
	  tangent1 = normal.cross(tangent0);
	  rot(0,0) = tangent0(0) ; rot(0,1) = tangent1(0) ; rot(0,2) = normal(0);
	  rot(1,0) = tangent0(1) ; rot(1,1) = tangent1(1) ; rot(1,2) = normal(1);
	  rot(2,0) = tangent0(2) ; rot(2,1) = tangent1(2) ; rot(2,2) = normal(2);
	  rotContact.block<3,3>(indexRom*3,0) = rot;
        
	  indexRom++;
	} // for each ROMS
      
      hppDout (info, "number of contacts: " << rbReport->ROMReports.size());
      polytope::vector_t x(rbReport->ROMReports.size());
      polytope::vector_t y(rbReport->ROMReports.size());
      polytope::vector_t nu(rbReport->ROMReports.size());
      const value_type xContact = contactSize [0];
      const value_type yContact = contactSize [1];
      for(size_t k = 0 ; k<rbReport->ROMReports.size() ; ++k){
        x(k) = xContact; // approx size of foot
        y(k) = yContact; 
        nu(k) = problem->mu_;
      }
      const polytope::ProjectedCone* giwc =
	polytope::U_stance (rotContact, posContact, nu, x, y);
      core::matrix_t Astance = polytope::A_stance(rotContact, posContact);
      hppDout (info, "Astance = " << Astance);
      //hppDout (info, "giwc->v = " << giwc->pImpl_->V_);
      return giwc;
    }

    // ---------------------------------------------------------------------

    /// Return list of contact-cone directions (contacts are in the 
    /// middle of ROM-obstacle intersection, using affordances.
    /// (Old MIG version: return "average" direction from contact-cones)

    namespace library {
      
      typedef struct ContactCones ContactCones;
      struct ContactCones {
	std::size_t coneNumber_;
	std::vector<std::string> ROMnames_;
	std::vector<fcl::Vec3f> directions_;
	std::vector<fcl::Vec3f> positions_;
      };

      inline fcl::Vec3f vectorToVec3f (const polytope::vector3_t vector) {
	fcl::Vec3f result;
	for (std::size_t i = 0; i < 3; i++) result [i] = vector [i];
	return result;
      }

      inline ContactCones computeContactCones 
      (const core::ProblemPtr_t& problem, const core::Configuration_t q) {
	ContactCones contactCones;
	std::vector<fcl::Vec3f> Cones, positions;
	std::vector<std::string> ROMnames;
	core::ValidationReportPtr_t report;
	const core::DevicePtr_t& robot (problem->robot ());
	model::RbPrmDevicePtr_t rbDevice =
	  boost::dynamic_pointer_cast<model::RbPrmDevice> (robot);
	if (!rbDevice) {
	  hppDout(error,"~~ Device cast in RB problem");
	  return contactCones;
	}

	const bool isValid = problem->configValidations()->validate(q,report);
	if(!isValid) {
	  hppDout(warning,"~~ config is not valid");
	  return contactCones;
	}
	if (!report) {
	  hppDout(error,"~~ Report problem");
	  return contactCones;
	}
	core::RbprmValidationReportPtr_t rbReport =
	  boost::dynamic_pointer_cast<core::RbprmValidationReport> (report);
	// checks :
	if(!rbReport)
	  {
	    hppDout(error,"~~ Validation Report cannot be cast");
	    return contactCones;
	  }

	//randomnize the collision pair, in order to get a different surface of contact each time (because only the first one in collision is considered by fcl and put in the report)
	problem->configValidations()->randomnizeCollisionPairs();
      
	const bool isSkelRobot = robot->name ().compare("skeleton_trunk") == 0;

	// get the 2 object in contact for each ROM :
	hppDout(info,"~~ Number of roms in collision : "<<rbReport->ROMReports.size());
	for(std::map<std::string,core::CollisionValidationReportPtr_t>::const_iterator it = rbReport->ROMReports.begin() ; it != rbReport->ROMReports.end() ; ++it)
	  {
	    hppDout(info,"~~ for rom : "<<it->first);
	    std::string ROMnameWithoutSphere = it->first;
	    ROMnameWithoutSphere.resize ((it->first).size () - 6);
	    hppDout(info,"ROMnameWithoutSphere= "<< ROMnameWithoutSphere);
	    ROMnames.push_back (ROMnameWithoutSphere);
	    core::CollisionObjectPtr_t obj1 = it->second->object1;
	    core::CollisionObjectPtr_t obj2 = it->second->object2;
	    hppDout(notice,"~~ collision between : "<<obj1->name() << " and "<<obj2->name());
	    fcl::CollisionResult result = it->second->result;
        
	    // get intersection between the two objects :
	    geom::T_Point vertices1;
	    geom::BVHModelOBConst_Ptr_t model1 =  geom::GetModel(obj1->fcl());
	    for(int i = 0 ; i < model1->num_vertices ; ++i)
	      {
		vertices1.push_back(Eigen::Vector3d(model1->vertices[i][0], model1->vertices[i][1], model1->vertices[i][2]));
	      }
        
	    geom::T_Point vertices2;
	    geom::BVHModelOBConst_Ptr_t model2 =  geom::GetModel(obj2->fcl());
	    for(int i = 0 ; i < model2->num_vertices ; ++i)
	      {
		vertices2.push_back(Eigen::Vector3d(model2->vertices[i][0], model2->vertices[i][1], model2->vertices[i][2]));
	      }
        
	    geom::Point pn; // normal
	    //geom::T_Point hull = geom::intersectPolygonePlane(model1,model2,pn); //my
	    geom::T_Point plane = geom::intersectPolygonePlane(model1,model2,pn); // pierre

	    geom::T_Point hull; // pierre
	    if(plane.size() > 0) // pierre
	      hull = geom::compute3DIntersection(plane,geom::convertBVH(model2)); // pierre
	    else
	      hppDout(error,"No intersection plane between rom and environnement");

	    if(hull.size() == 0){
	      hppDout(error,"No intersection between rom and environnement");
	    }else{
	      geom::Point center = geom::center(hull.begin(),hull.end());
	      hppDout(notice,"Center = "<<center.transpose());
	      hppDout(notice,"Normal : "<<pn.transpose());
	      positions.push_back (vectorToVec3f (center));
	      polytope::vector3_t normal = pn;
	      normal.normalize ();
	      fcl::Vec3f normal_vec3f = vectorToVec3f (normal);
	      // TODO: check normal direction before pushing it!
	      fcl::Vec3f robotCenter;
	      for (int i = 0; i < 3; i++) robotCenter [i] = q [i];
	      if ((robotCenter - vectorToVec3f (center)).dot (normal_vec3f) < 0 && isSkelRobot) {
		hppDout(info,"Normal has been inverted");
		normal_vec3f = -normal_vec3f;
	      }
	      Cones.push_back (normal_vec3f);
	    }
	  } // for each ROMS
	hppDout(notice,"positions.size()= "<< positions.size());

	// DEBUG !!! Fake cone !! for Skeleton-parkour
	/*const size_type extraDim = robot->extraConfigSpace ().dimension ();
	  const size_type ecIndex = robot->configSize() - extraDim;
	  fcl::Vec3f dir, pos;
	  for (int k = 0; k < 3; k++) dir [k] = q [ecIndex + k];
	  for (int k = 0; k < 3; k++) pos [k] = q [k];
	  std::vector<fcl::Vec3f> dirs, poss;
	  dirs.push_back (dir); poss.push_back (pos);
	  std::vector<std::string> ROMs;
	  ROMs.push_back (ROMnames [0]);
	  contactCones.coneNumber_ = 1;
	  contactCones.directions_ = dirs;
	  contactCones.positions_ = poss;
	  contactCones.ROMnames_ = ROMs;*/

	// normal code
	contactCones.coneNumber_ = Cones.size ();
	contactCones.directions_ = Cones;
	contactCones.positions_ = positions;
	contactCones.ROMnames_ = ROMnames;

	hppDout(notice,"positions_.size ()= "<< contactCones.positions_.size ());
	hppDout(notice, "ROMnames_.size ()= "<< contactCones.ROMnames_.size ());
	hppDout(notice, "ROMnames_[0]= "<< contactCones.ROMnames_[0]);

	return contactCones;
      }
    } // namespace library
  } //   namespace rbprm
} // namespace hpp

#endif // HPP_RBPRM_PARABOLA_LIBRARY_HH
