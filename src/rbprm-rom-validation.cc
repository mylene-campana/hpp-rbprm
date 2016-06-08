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

#include "hpp/rbprm/rbprm-rom-validation.hh"
#include <hpp/fcl/collision.h>
#include <hpp/fcl/BVH/BVH_model.h>
#include <hpp/rbprm/rbprm-validation-report.hh>
#include "utils/algorithms.h"

namespace hpp {
  using namespace core;
  namespace rbprm {

    RbPrmRomValidationPtr_t RbPrmRomValidation::create
    (const model::DevicePtr_t& robot, const NormalFilter& normalFilter)
    {
      RbPrmRomValidation* ptr = new RbPrmRomValidation (robot, normalFilter);
      return RbPrmRomValidationPtr_t (ptr);
    }

    RbPrmRomValidation::RbPrmRomValidation (const model::DevicePtr_t& robot
                                           ,const NormalFilter& normalFilter)
        : hpp::core::CollisionValidation(robot)
        , filter_(normalFilter)
        , unusedReport_(new CollisionValidationReport)
    {
        if(!normalFilter.unConstrained_)
        {
            collisionRequest_.enable_contact = true;
            collisionRequest_.num_max_contacts = 30;
        }
    }


    bool RbPrmRomValidation::validate (const Configuration_t& config)
    {
        return validate(config, unusedReport_);
    }

    bool RbPrmRomValidation::validate (const Configuration_t& config,
                    ValidationReportPtr_t& validationReport)
    {        
        bool collision = !hpp::core::CollisionValidation::validate(config, validationReport);
        if(collision && !filter_.unConstrained_)
        {
            collision = false;
            CollisionValidationReport* report =
                    static_cast <CollisionValidationReport*> (validationReport.get());
            for(std::size_t i = 0; i< report->result.numContacts() && !collision; ++i)
            {
                // retrieve triangle
                const fcl::Contact& contact =  report->result.getContact(i);
                assert(contact.o2->getObjectType() == fcl::OT_BVH); // only works with meshes
                const fcl::BVHModel<fcl::OBBRSS>* surface = static_cast<const fcl::BVHModel<fcl::OBBRSS>*> (contact.o2);
                const fcl::Triangle& tr = surface->tri_indices[contact.b2];
                const fcl::Vec3f& v1 = surface->vertices[tr[0]];
                const fcl::Vec3f& v2 = surface->vertices[tr[1]];
                const fcl::Vec3f& v3 = surface->vertices[tr[2]];
                fcl::Vec3f normal = (v2 - v1).cross(v3 - v1);
                normal.normalize();
                if(normal.dot(filter_.normal_)>=filter_.range_){
                    collision = true;
                }
            }
	    /* Pierre work
            // test contact area size :  
            //TODO remove it after the precomputation of the meshs is working
           if(collision){
              geom::BVHModelOBConst_Ptr_t model1 =  geom::GetModel(colReport->object1->fcl());
              geom::BVHModelOBConst_Ptr_t model2 =  geom::GetModel(colReport->object2->fcl());
              geom::T_Point intersection = geom::intersectPolygonePlane(model1,model2,filter_.normal_,0,colReport->result,false,filter_.range_);
              
              double a = geom::area(intersection.begin(),intersection.end());
              if(a <= 0.1)
                collision = false;
		}*/
        }
        return collision;
    }
  }// namespace rbprm
}// namespace hpp
