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


    /// Arrange robot orientation according to the surface normal direction
    /// So that the robot is "on" the surface, rotated
    inline core::Configuration_t setOrientation
    (const core::DevicePtr_t& robot, const core::Configuration_t& q) {
      core::Configuration_t qtest = q;
      //hppDout (info, "q: " << displayConfig (q));
      const core::JointPtr_t jointSO3 = robot->getJointVector () [1];
      const size_type indexSO3 = jointSO3->rankInConfiguration ();
      const size_type index = robot->configSize ()
	- robot->extraConfigSpace ().dimension ();

      const value_type nx = q [index];
      const value_type ny = q [index + 1];
      const value_type nz = q [index + 2];
      //const value_type theta = atan2 (2*qw*qz - 2*qx*qy, 1-2*qy*qy-2*qz*qz);
      const value_type theta = q [index + 3];
      //hppDout (info, "theta: " << theta);

      const int sign = -(theta > M_PI/2) - (theta < -M_PI/2) +
	(int) ((theta > -M_PI/2) && (theta < M_PI/2));
      
      value_type x12, y12, z12; // corresponding to R_r coordinates

      // See Matlab script or Latex doc for the different cases:
      if (theta != M_PI/2 && theta != -M_PI/2) {
	const value_type tantheta = tan(theta);
	if (nz != 0) {
	  x12= sign*nz/sqrt((1+tantheta*tantheta)
			    *nz*nz+(nx+ny*tantheta)*(nx+ny*tantheta));
	  y12 = x12*tantheta;
	  z12 = -x12*(nx+ny*tantheta)/nz;
	} else { // nz = 0 case
	  x12 = 0; y12 = 0; z12 = 1;
	}
      } else { // theta = +- pi/2
	x12 = 0;
	if (ny != 0) {
	  if (theta == M_PI/2)
	    z12 = 1/(1+(nz/ny)*(nz/ny));
	  else // theta = -pi/2
	    z12 = 1/(1+(nz/ny)*(nz/ny));
	  y12 = -nz/ny*z12;
	} else { // ny = 0
	  if (nx == 1) { // no matter y12 and z12)
	    y12 = 0.5;
	    z12 = 0.5;
	  } else { // nx = 1
	    y12 = 1;
	    z12 = 0;
	  }
	}
      }

      // cross-product
      const value_type zx = nz*x12-nx*z12;
      const value_type yz = ny*z12-nz*y12;
      const value_type xy = nx*y12-ny*x12;
      
      //fcl::Matrix3f A; // A: rotation matrix expressing R_r in R_0
      Eigen::Matrix<value_type,3,3> A;
      A (0,0) = x12; A (0,1) = yz; A (0,2) = nx;
      A (1,0) = y12; A (1,1) = zx; A (1,2) = ny;
      A (2,0) = z12; A (2,1) = xy; A (2,2) = nz;
      //hppDout (info, "A: " << A);
      //hppDout (info, "A.determinant (): " << A.determinant ());

      Eigen::Quaternion<value_type> quat (A);
      const value_type qw = quat.w ();
      const value_type qx = quat.x ();
      const value_type qy = quat.y ();
      const value_type qz = quat.z ();
      const value_type magnitude = sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
      hppDout (info, "quat: " << " " << qw << " " << qx << " " << qy << " " << qz);
      hppDout (info, "quaternion magnitude= " << magnitude);
      
      // re-normalize (not needed but sometimes, loss of accuracy...)
      qtest [indexSO3] = qw / magnitude;
      qtest [indexSO3 + 1] = qx / magnitude;
      qtest [indexSO3 + 2] = qy / magnitude;
      qtest [indexSO3 + 3] = qz / magnitude;
      hppDout (info, "qtest: " << displayConfig (qtest));
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
        
	  /// warning: plane normal harcoded to (0,0,1) here ->
	  geom::T_Point hull = geom::intersectPolygonePlane(model1,model2,fcl::Vec3f(0,0,1),geom::ZJUMP,result);

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

  } //   namespace rbprm
} // namespace hpp

#endif // HPP_RBPRM_PARABOLA_LIBRARY_HH
