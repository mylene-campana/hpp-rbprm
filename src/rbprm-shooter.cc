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
#include <hpp/core/collision-validation-report.hh>
#include <hpp/rbprm/rbprm-shooter.hh>
#include <hpp/model/collision-object.hh>
#include <hpp/model/joint.hh>
#include <hpp/fcl/collision_object.h>
#include <hpp/fcl/BVH/BVH_model.h>
#include <hpp/core/collision-validation.hh>
#include <Eigen/Geometry>
#include <hpp/model/configuration.hh>
#include <hpp/rbprm/fullbodyBallistic/parabola-library.hh>

namespace hpp {
using namespace core;
using namespace fcl;
namespace
{
    typedef fcl::BVHModel<OBBRSS> BVHModelOB;
    typedef boost::shared_ptr<const BVHModelOB> BVHModelOBConst_Ptr_t;

    BVHModelOBConst_Ptr_t GetModel(const fcl::CollisionObjectConstPtr_t object)
    {
        assert(object->collisionGeometry()->getNodeType() == BV_OBBRSS);
        const BVHModelOBConst_Ptr_t model = boost::static_pointer_cast<const BVHModelOB>(object->collisionGeometry());
        assert(model->getModelType() == BVH_MODEL_TRIANGLES);
        return model;
    }

    double TriangleArea(rbprm::TrianglePoints& tri)
    {
        double a, b, c;
        a = (tri.p1 - tri.p2).norm();
        b = (tri.p2 - tri.p3).norm();
        c = (tri.p3 - tri.p1).norm();
        double s = 0.5 * (a + b + c);
        return sqrt(s * (s-a) * (s-b) * (s-c));
    }

    std::vector<double> getTranslationBounds(const model::RbPrmDevicePtr_t robot)
    {
        const JointPtr_t root = robot->rootJoint();
        std::vector<double> res;
        for(std::size_t i =0; i<3; ++i)
        {
            if(root->isBounded(i))
            {
                res.push_back(root->lowerBound(i));
                res.push_back(root->upperBound(i));
            }
            else            {

                res.push_back(-std::numeric_limits<double>::max());
                res.push_back(std::numeric_limits<double>::max());
            }
        }
        return res;
    }

    void SetConfigTranslation(const model::RbPrmDevicePtr_t robot, ConfigurationPtr_t config, const Vec3f& translation)
    {
        std::vector<double> bounds = getTranslationBounds(robot);
        for(std::size_t i =0; i<3; ++i)
        {
            (*config)(i)= std::min(bounds[2*i+1], std::max(bounds[2*i], translation[i]));
        }
    }

    void Translate(const model::RbPrmDevicePtr_t robot, ConfigurationPtr_t config, const Vec3f& translation)
    {
        // bound to positions limits
        std::vector<double> bounds = getTranslationBounds(robot);
        for(int i =0; i<3; ++i)
        {
            (*config)(i)=std::min(bounds[2*i+1], std::max(bounds[2*i], (*config)(i)+translation[i]));
        }
    }

    void SampleRotationRec(ConfigurationPtr_t config, JointVector_t& jv, std::size_t& current)
    {
        JointPtr_t joint = jv[current++];
        std::size_t rank = joint->rankInConfiguration ();
        joint->configuration ()->uniformlySample (rank, *config);
        if(current<jv.size())
            SampleRotationRec(config,jv,current);
    }

    void SampleRotation(model::DevicePtr_t so3, ConfigurationPtr_t config, JointVector_t& /*jv*/)
    {
      size_t id = 0;
      if(so3->rootJoint())
      {
          std::size_t rank = 3;
          //so3->rootJoint()->configuration()->uniformlySample(rank,*config);
          Eigen::Matrix <value_type, 3, 1> confso3;
        
          id+=1;
          model::JointPtr_t joint = so3->rootJoint();
          for(int i =0; i <3; ++i)
          {
             //hppDout(notice,"JOINT : name = "<<joint->name()<<" size = "<<joint->configSize());
             //hppDout(notice,"JOINT : lower bound = "<<joint->lowerBound(0)<<"upper bound = "<<joint->upperBound(0));
              joint->configuration()->uniformlySample (i, confso3);
              if(i<2)
                joint = joint->childJoint(0);
          }
          //std::cout<<confso3<<std::endl;
          Eigen::Quaterniond qt = Eigen::AngleAxisd(confso3(0), Eigen::Vector3d::UnitZ())
            * Eigen::AngleAxisd(confso3(1), Eigen::Vector3d::UnitY())
            * Eigen::AngleAxisd(confso3(2), Eigen::Vector3d::UnitX());
          //std::cout<<"quat = "<<qt.w()<<" , "<<qt.x()<<" , "<<qt.y()<<" , "<<qt.z()<<" , "<<std::endl;
          (*config)(rank+0) = qt.w();
          (*config)(rank+1) = qt.x();
          (*config)(rank+2) = qt.y();
          (*config)(rank+3) = qt.z();
          //std::cout<<model::displayConfig(*config)<<std::endl;
          
      }
     /* if(id < jv.size()){
          std::cout<<"SampleRotationRec"<<std::endl;
          SampleRotationRec(config,jv,id);
      }   */
        
    }

    model::DevicePtr_t initSo3()
    {
        DevicePtr_t so3Robot = model::Device::create("so3Robot");
        /*model::JointPtr_t res = new model::JointSO3(fcl::Transform3f());
        res->name("defaultSO3");*/
        //so3Robot->rootJoint(0);
        return so3Robot;
    }

    void seRotationtLimits (model::DevicePtr_t so3Robot, const std::vector<double>& limitszyx)
    {
        model::Joint* previous = so3Robot->rootJoint();
        if(previous == 0)
        {
            // init joints
            previous = new model::jointRotation::Bounded(fcl::Transform3f());
            model::Joint * jy = new model::jointRotation::Bounded(fcl::Transform3f());
            model::Joint * jx = new model::jointRotation::Bounded(fcl::Transform3f());
            so3Robot->rootJoint(previous);
            previous->addChildJoint (jy);
            jy->addChildJoint (jx);
            previous->name("so3z");
            jy->name("so3y");
            jx->name("so3x");
        }
        // set limits model::JointPtr_t
        assert(limitszyx.size() == 6);
        int i = 0;
        model::JointPtr_t current = previous;
        for(std::vector<double>::const_iterator cit = limitszyx.begin();
            cit != limitszyx.end(); ++cit, ++i)
        {
            if(i % 2 == 0)
            {
                current->lowerBound(0, *cit);
            }
            else
            {
                current->upperBound(0, *cit);
                
                current = current->numberChildJoints() > 0 ? current->childJoint(0) : 0;
            }
        }
    }

} // namespace

  namespace rbprm {

    RbPrmShooterPtr_t RbPrmShooter::create (const model::RbPrmDevicePtr_t& robot,
                                            const ObjectVector_t& geometries,
                                            const std::vector<std::string>& filter,
                                            const std::map<std::string, rbprm::NormalFilter>& normalFilters,
                                            const std::size_t shootLimit, const std::size_t displacementLimit,
					    const std::size_t nbFilterMatch)
    {
       /* unsigned int seed = (unsigned int)(time(NULL));
        srand (seed);
        hppDout(notice,"&&&&&& SEED = "<<seed);
*/
     // srand (0);

        RbPrmShooter* ptr = new RbPrmShooter (robot, geometries, filter, normalFilters, shootLimit, displacementLimit, nbFilterMatch);
        RbPrmShooterPtr_t shPtr (ptr);
        ptr->init (shPtr);
        return shPtr;
    }

    void RbPrmShooter::init (const RbPrmShooterPtr_t& self)
    {
        ConfigurationShooter::init (self);
        weak_ = self;
    }


    void RbPrmShooter::BoundSO3(const std::vector<double>& limitszyx)
    {
        seRotationtLimits(eulerSo3_, limitszyx);
    }

    // TODO: outward

    RbPrmShooter::RbPrmShooter (const model::RbPrmDevicePtr_t& robot,
				const ObjectVector_t& geometries,
				const std::vector<std::string>& filter,
				const std::map<std::string, rbprm::NormalFilter>& normalFilters,
				const std::size_t shootLimit,
				const std::size_t displacementLimit,
				const std::size_t nbFilterMatch)
      : shootLimit_(shootLimit)
      , displacementLimit_(displacementLimit)
      , filter_(filter)
      , robot_ (robot)
      , validator_(rbprm::RbPrmValidation::create(robot_, filter, normalFilters, nbFilterMatch))
      , eulerSo3_(initSo3())
      , fullOrientationMode_ (false)
    {
      hppDout (info, "constructor RbPrmShooter");
      for (std::size_t i = 0; i < filter_.size (); i++) {
      hppDout (info, "rbShooter filter= " << filter_[i]);
    }
        for(hpp::core::ObjectVector_t::const_iterator cit = geometries.begin();
            cit != geometries.end(); ++cit)
        {
            validator_->addObstacle(*cit);
        }
        this->InitWeightedTriangles(geometries);
    }

    void RbPrmShooter::InitWeightedTriangles(const model::ObjectVector_t& geometries)
    {
        double sum = 0;
        for(model::ObjectVector_t::const_iterator objit = geometries.begin();
          objit != geometries.end(); ++objit)
        {
            const  fcl::CollisionObjectPtr_t& colObj = (*objit)->fcl();
            BVHModelOBConst_Ptr_t model =  GetModel(colObj); // TODO NOT TRIANGLES
            for(int i =0; i < model->num_tris; ++i)
            {
                TrianglePoints tri;
                Triangle fcltri = model->tri_indices[i];
                tri.p1 = colObj->getRotation() * model->vertices[fcltri[0]] + colObj->getTranslation();
                tri.p2 = colObj->getRotation() * model->vertices[fcltri[1]] + colObj->getTranslation();
                tri.p3 = colObj->getRotation() * model->vertices[fcltri[2]] + colObj->getTranslation();;
                double weight = TriangleArea(tri);
                sum += weight;
                weights_.push_back(weight);
                fcl::Vec3f normal = (tri.p2 - tri.p1).cross(tri.p3 - tri.p1);
                normal.normalize();
                triangles_.push_back(std::make_pair(normal,tri));
            }
            double previousWeight = 0;
            for(std::vector<double>::iterator wit = weights_.begin();
                wit != weights_.end(); ++wit)
            {
                previousWeight += (*wit) / sum;
                (*wit) = previousWeight;
            }
        }
    }


  const RbPrmShooter::T_TriangleNormal &RbPrmShooter::RandomPointIntriangle() const
  {
      return triangles_[rand() % triangles_.size()];
  }

  const RbPrmShooter::T_TriangleNormal& RbPrmShooter::WeightedTriangle() const
  {
      double r = ((double) rand() / (RAND_MAX));
      hppDout (info, "r= " << r);
      std::vector<T_TriangleNormal>::const_iterator trit = triangles_.begin();
      for(std::vector<double>::const_iterator wit = weights_.begin();
          wit != weights_.end();
          ++wit, ++trit)
      {
          if(*wit >= r) return *trit;
      }
      return triangles_[triangles_.size()-1]; // not supposed to happen
  }

hpp::core::ConfigurationPtr_t RbPrmShooter::shoot () const
{
    const size_type extraDim = robot_->extraConfigSpace ().dimension ();
    const size_type index = robot_->configSize() - extraDim;
    const bool hasECS = extraDim >= 3;
    JointVector_t jv = robot_->getJointVector ();
    ConfigurationPtr_t config (new Configuration_t (robot_->Device::currentConfiguration()));
    std::size_t limit = shootLimit_;
    bool found(false);
    CollisionValidationReport* report;
    ValidationReportPtr_t reportShPtr(new CollisionValidationReport);
    value_type thetaSample = 0;
    Vec3f normal;
    hppDout (info, "fullOrientationMode in shooter= " << fullOrientationMode_);
    bool valid = false;
    while(limit >0 && (!found || !valid))
    {
        // pick one triangle randomly
        const T_TriangleNormal* sampled(0);
        double r = ((double) rand() / (RAND_MAX));
	hppDout (info, "r= " << r);
        if(r > 0.5)
            sampled = &RandomPointIntriangle();
        else
            sampled = &WeightedTriangle();
        const TrianglePoints& tri = sampled->second;
        //http://stackoverflow.com/questions/4778147/sample-random-point-in-triangle
        double r1, r2;
        r1 = ((double) rand() / (RAND_MAX)); r2 = ((double) rand() / (RAND_MAX));
        Vec3f p = (1 - sqrt(r1)) * tri.p1 + (sqrt(r1) * (1 - r2)) * tri.p2
                + (sqrt(r1) * r2) * tri.p3;
	hppDout (info, "r1= " << r1);
	hppDout (info, "r2= " << r2);
	hppDout (info, "p= " << p);
	// get normal even if not in collision...
	normal = (tri.p2 - tri.p1).cross(tri.p3 - tri.p1);
	normal.normalize();
	hppDout (info, "normal= " << normal);
        //set configuration position to sampled point
        SetConfigTranslation(robot_,config, p);
	//hppDout (info, "config= " << displayConfig (*config));
	if (hasECS && fullOrientationMode_) {
	  if (!validator_->trunkValidation_->validate(*config, reportShPtr)) {
	    report = static_cast<CollisionValidationReport*>(reportShPtr.get());
	    Vec3f normal_test = triangles_[report->result.getContact(0).b2].first;
	    hppDout (info, "normal_test= " << normal_test);
	  }
	  for (size_type i=0; i<3; i++) (*config) [index + i] = normal [i];
	  thetaSample = 2 * M_PI * rand ()/RAND_MAX - M_PI;
	  (*config) [index + 3] = thetaSample;
	  *config = setOrientation (robot_, *config);
	}
	else {
	  SampleRotation(eulerSo3_, config, jv);
	  hppDout (info, "random rotation was sampled");
	}
	//hppDout (info, "config= " << displayConfig (*config));
        // rotate and translate randomly until valid configuration found or
        // no obstacle is reachable
        std::size_t limitDis = displacementLimit_;
        Vec3f lastDirection = normal; // will be updated if necessary
        while(!found && limitDis >0)
        {
            valid = validator_->trunkValidation_->validate(*config, reportShPtr);
            report = static_cast<CollisionValidationReport*>(reportShPtr.get());
            found = valid && validator_->validateRoms(*config, filter_);
            if(valid &!found)
            {
                // try to rotate to reach rom
                for(; limitDis>0 && !found; --limitDis)
                {
		  //SampleRotation(eulerSo3_, config, jv);
		  if (hasECS && fullOrientationMode_) {
		    thetaSample = 2 * M_PI * rand ()/RAND_MAX - M_PI;
		    hppDout (info, "thetaSample= " << thetaSample);
		    for (size_type i=0; i<3; ++i)
		      (*config) [index + i] = normal [i];
		    (*config) [index + 3] = thetaSample;
		    *config = setOrientation (robot_, *config);
		    hppDout (info, "normal= " << normal);
		  } else {
		    SampleRotation(eulerSo3_, config, jv);
		    hppDout (info, "random rotation was sampled");
		  }
		  //hppDout (info, "config= " << displayConfig(*config));
                    found = validator_->validate(*config, filter_);
                    if(!found)
                    {
                        Translate(robot_, config, -lastDirection *
                                  1 * ((double) rand() / (RAND_MAX)));
                    }
                    found = validator_->validate(*config, filter_);
                }
                if(!found) break;
            }
            else if (!valid)// move out of collision
            {
                // retrieve Contact information
                //lastDirection = -report.result.getContact(0).normal;
                // mouve out by penetration depth
                // v0 move away from normal
                //get normal from collision tri
                lastDirection = triangles_[report->result.getContact(0).b2].first;
                Translate(robot_,config, lastDirection *
                          (std::abs(report->result.getContact(0).penetration_depth) +0.03));
                 limitDis--;
            }
        }
        // Shoot extra configuration variables
        /*for (size_type i=0; i<extraDim; ++i)
        {
            value_type lower = robot_->extraConfigSpace ().lower (i);
            value_type upper = robot_->extraConfigSpace ().upper (i);
            value_type range = upper - lower;
            if ((range < 0) ||
              (range == std::numeric_limits<double>::infinity()))
            {
                std::ostringstream oss
                  ("Cannot uniformy sample extra config variable ");
                oss << i << ". min = " << ", max = " << upper << std::endl;
                throw std::runtime_error (oss.str ());
            }
            (*config) [index + i] = (upper - lower) * rand ()/RAND_MAX + lower;
        }
	*/
        // save the normal (code from Mylène)
        if(extraDim >= 3 ){
	  hppDout (info, "lastDirection= " << lastDirection);
          for (size_type i=0; i<3; ++i)
            (*config) [index + i] = lastDirection [i];
	  if (fullOrientationMode_)
	    *config = setOrientation (robot_, *config);
	  }
	valid = validator_->trunkValidation_->validate(*config, reportShPtr);
        limit--;
    }
    if (!found) std::cout << "no config found" << std::endl;
    return config;
}


  }// namespace rbprm
}// namespace hpp
