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

#include <hpp/util/debug.hh>
#include <hpp/model/configuration.hh>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/joint-configuration.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/rbprm/fullbodyBallistic/timed-ballistic-path.hh>
#include <hpp/core/straight-path.hh>
#include <Eigen/Geometry>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::value_type;
    using core::vector_t;
    using core::interval_t;
    using core::size_type;
    
    TimedBallisticPath::TimedBallisticPath (const core::DevicePtr_t& device,
                                            core::ConfigurationIn_t init,
                                            core::ConfigurationIn_t end,
                                            value_type length,
                                            value_type alpha, value_type theta, value_type v0,BallisticPathPtr_t bp) :
      parent_t (interval_t (0, length), device->configSize (),
                device->numberDof ()), device_ (device), initial_ (init),
      end_ (end), length_ (length),alpha_(alpha),theta_(theta),v0_(v0),g_(9.81),bp_(bp)
    {
      assert (device);
      xTheta0_ = initial_[0]*cos(theta_) + initial_[1]*sin(theta_);
    }
    
    TimedBallisticPath::TimedBallisticPath (const rbprm::BallisticPathPtr_t ballisticPath) :
      parent_t (interval_t(0,ballisticPath->length()),ballisticPath->device()->configSize(),ballisticPath->device()->numberDof ()),device_(ballisticPath->device()),initial_(ballisticPath->initial()),end_(ballisticPath->end()),length_(ballisticPath->length()),g_(9.81),lastRootIndex_(ballisticPath->lastRootIndex()),bp_(ballisticPath)
    {
      // TODO : get coeffs
      alpha_ = ballisticPath->coefficients()[4];
      theta_ = ballisticPath->coefficients()[3];
      const value_type x_theta_0_dot = ballisticPath->coefficients()[5];
      v0_ = sqrt((1 + tan(alpha_)*tan(alpha_))) * x_theta_0_dot;
      hppDout(notice,"Create timed path : alpha = "<<alpha_<<"  theta = "<<theta_<<"  v0 = "<<v0_);
      xTheta0_ = initial_[0]*cos(theta_) + initial_[1]*sin(theta_);
      length_ = computeLength(initial_,end_);
    }
    
    TimedBallisticPath::TimedBallisticPath (const BallisticPathPtr_t bp1,const BallisticPathPtr_t bp1Max,const BallisticPathPtr_t bp2Max,const BallisticPathPtr_t bp2) :
      parent_t (interval_t(0,bp1->length()),bp1->device()->configSize(),bp1->device()->numberDof ()),device_(bp1->device()),initial_(bp1->initial()),end_(bp2->end()),length_(bp1->length()+bp1Max->length()+bp2Max->length()+bp2->length()),g_(9.81),bp1_(bp1),bp1Max_(bp1Max),bp2Max_(bp2Max),bp2_(bp2){
      // TODO : get coeffs
      alpha_ = bp1->coefficients()[4];
      theta_ = bp1->coefficients()[3];
      const value_type x_theta_0_dot = bp1->coefficients()[5];
      v0_ = sqrt((1 + tan(alpha_)*tan(alpha_))) * x_theta_0_dot;
      hppDout(notice,"Create timed path : alpha = "<<alpha_<<"  theta = "<<theta_<<"  v0 = "<<v0_);
      xTheta0_ = initial_[0]*cos(theta_) + initial_[1]*sin(theta_);
      length_ = computeLength(initial_,end_);
      
      // compute lenght of each subpath:
      value_type xTheta1 = cos(theta_)*bp1->end()[0] + sin(theta_)*bp1->end()[1] - xTheta0_;
      value_type xTheta1Max = cos(theta_)*bp1Max->end()[0] + sin(theta_)*bp1Max->end()[1] - xTheta0_;
      value_type xTheta2Max = cos(theta_)*bp2Max->end()[0] + sin(theta_)*bp2Max->end()[1] - xTheta0_;
      value_type xTheta2 = cos(theta_)*bp2->end()[0] + sin(theta_)*bp2->end()[1] - xTheta0_;
      hppDout(notice,"Bounds on xTheta for each subpath : "<<xTheta1<< " , "<<xTheta1Max<< " , "<<xTheta2Max<< " , "<<xTheta2<< " , ");
      t1_ = xTheta1/(v0_*cos(alpha_));
      t1Max_ = xTheta1Max/(v0_*cos(alpha_));
      t2Max_ = xTheta2Max/(v0_*cos(alpha_));
      t2_ = xTheta2/(v0_*cos(alpha_));
      hppDout(notice,"Bounds on t for each subpath : "<<t1_<< " , "<<t1Max_<< " , "<<t2Max_<< " , "<<t2_<< " , ");
    }
    
    TimedBallisticPath::TimedBallisticPath (const TimedBallisticPath& path) :
      parent_t (path), device_ (path.device_), initial_ (path.initial_),
      end_ (path.end_),length_ (path.length_),alpha_(path.alpha_),theta_(path.theta_),v0_(path.v0_),lastRootIndex_(path.lastRootIndex_)
    {
      xTheta0_ = initial_[0]*cos(theta_) + initial_[1]*sin(theta_);
    }
    
    bool TimedBallisticPath::impl_compute (core::ConfigurationOut_t result,
                                           value_type t) const
    {
        if (t == 0 || initial_(0) == end_(0)) {
          for(size_t i = 0 ; i < device_->configSize(); i++){
            result[i] = initial_[i];
          }
          return true;
        }
        if (t >= length_) {
          for(size_t i = 0 ; i < device_->configSize(); i++){
            result[i] = end_[i] ;
          }
          return true;
        }
      value_type xTheta = v0_*cos(alpha_)*t ;
      //hppDout(info,"xTheta = "<<xTheta);
      value_type u; // param in BallisticPath representation (need convertion from param t)
      // interpolation for the articulation are done by ballisticPath (one of the 4 subpath)
    /*  if(t < t1_){
        u = (t/t1_) * bp1_->length();
        (*bp1_)(result,u);
        //hppDout(notice,"1 : u ="<<u<<"  lenght = "<<bp1_->length());
      }else if(t < t1Max_){
        u = ((t-t1_)/(t1Max_- t1_)) * bp1Max_->length();
        (*bp1Max_)(result,u);
        //hppDout(notice,"1Max : u ="<<u<<"  lenght = "<<bp1Max_->length());
      }else if(t < t2Max_){
        u = ((t-t1Max_)/(t2Max_-t1Max_)) * bp2Max_->length();
        (*bp2Max_)(result,u);
        //hppDout(notice,"2Max : u ="<<u<<"  lenght = "<<bp2Max_->length());
      }else{
        u = ((t-t2Max_)/(t2_-t2Max_)) * bp2_->length();
        (*bp2_)(result,u);
        //hppDout(notice,"2 : u ="<<u<<"  lenght = "<<bp2_->length());        
      }*/
      
      
        
      // replace with the correct position / orientation for the center

      result[0] = xTheta*cos(theta_) + initial_[0];
      result[1] = xTheta*sin(theta_) + initial_[1];     
      result[2] = -0.5*g_*t*t + v0_*sin(alpha_)*t + initial_[2];
      
      
    // std::cout<<"t = "<<t;
     // std::cout<<"   ; z0 = "<<initial_[2]<< "   ; z = "<<result[2]<<"  ; alpha0 = "<<alpha_<<std::endl;
      
      
      /* Quaternions interpolation */
      u = t/length_;
      const core::JointPtr_t SO3joint = device_->getJointByName ("base_joint_SO3");
      const std::size_t rank = SO3joint->rankInConfiguration ();
      const core::size_type dimSO3 = SO3joint->configSize ();
      SO3joint->configuration ()->interpolate
	(initial_, end_, u, rank, result);

      const size_type nbConfig = device_->configSize();
      const size_type ecsDim = device_->extraConfigSpace ().dimension ();
      
      for(size_t i = 7 ; i < nbConfig-ecsDim ; i++){
          result (i) = (1 - u) * initial_ (i) + u * end_ (i);
      }
      return true;
    }
    
    
    core::PathPtr_t TimedBallisticPath::extract (const interval_t& subInterval) const throw (hpp::core::projection_error)
    {
      bool success;
      hppDout(notice,"extract timed path : "<<subInterval.second - subInterval.first);
      core::Configuration_t q1 ((*this) (subInterval.first, success)); // straight
      core::Configuration_t q2 ((*this) (subInterval.second, success)); // straight
      interval_t bpInterval(bp_->computeLength(bp_->initial(),q1),bp_->computeLength(bp_->initial(),q2));
      core::PathPtr_t extractBp = bp_->extract(bpInterval);
      core::PathPtr_t result = rbprm::TimedBallisticPath::create(device_,q1,q2,computeLength(q1,q2),alphaAtParam(subInterval.first),theta_,velocityAtParam(subInterval.first),
                                                                 boost::dynamic_pointer_cast<BallisticPath>(extractBp),lastRootIndex_);
      return result;
    }
    
    core::PathPtr_t TimedBallisticPath::reverse () const{
      hppDout(notice, "reverse path parabola !!!!!!!!!!!!!!!!!!!!!!!!");
      bool success;
      core::Configuration_t q1 ((*this) (length_, success));
      core::Configuration_t q2 ((*this) (0, success));
      core::PathPtr_t result = TimedBallisticPath::create (device_, q1, q2, length_,
                                                           alpha_,theta_,v0_,bp_);
      return result;
    }
    
    core::DevicePtr_t TimedBallisticPath::device () const
    {
      return device_;
    }

    value_type TimedBallisticPath::velocityAtParam(double t) const{
      value_type xTheta = v0_*cos(alpha_)*t;
      value_type v = sqrt(v0_*v0_-2.*g_*xTheta*tan(alpha_) + (g_*xTheta)/(v0_*cos(alpha_))*(g_*xTheta)/(v0_*cos(alpha_)));
      return v;
    }
    
    value_type TimedBallisticPath::alphaAtParam(value_type t) const{
      Eigen::Vector3d dir = evaluateVelocity(t);
      Eigen::Vector3d xTheta(cos(theta_),sin(theta_),0);
      dir.normalize();
      xTheta.normalize();
      hppDout(notice,"direction for alpha = "<<dir);
      hppDout(notice,"xTheta for alpha = "<<xTheta);
      value_type alphaT = acos(xTheta.dot(dir));
      if(dir[2] < 0 )
        alphaT = -alphaT;
      hppDout(notice,"Current alpha = "<<alphaT<< " with alpha0 = "<<alpha_<<"  ; at t ="<<t);
      return alphaT;
    }
    
    value_type TimedBallisticPath::computeLength
    (const core::ConfigurationIn_t q1, const core::ConfigurationIn_t q2) const {
      value_type z = q1[2] - q2[2]; // difference of height
      value_type xTheta1 = q1[0]*cos(theta_) + q1[1]*sin(theta_);
      value_type xTheta2 = q2[0]*cos(theta_) + q2[1]*sin(theta_);
      value_type t1 = xTheta1 / (v0_*cos(alpha_));
      value_type t2 = xTheta2 / (v0_*cos(alpha_));

      value_type lenght = std::fabs(t2 - t1);
      hppDout(info,"difference of height = "<<z);
      hppDout(notice, "total flying time = "<<lenght);
      hppDout(notice,"xTheta1 = "<<xTheta1<< "   , xTheta2 = "<<xTheta2);
      return lenght;
    }
    
    

    vector_t TimedBallisticPath::evaluateVelocity (const value_type t) const {
      vector_t vel (3);
      bool success;
      const value_type theta = bp_->coefficients()[3];
      const value_type alpha = bp_->coefficients()[4];;
      const value_type x_theta_0_dot = bp_->coefficients()[5];;
      const value_type inv_x_theta_0_dot_sq = 1/(x_theta_0_dot*x_theta_0_dot);
      const value_type x_theta_0 = bp_->coefficients()[6];
      const core::Configuration_t q = (*this) (t, success);
      const value_type x_theta = q [0]*cos(theta) + q [1]*sin(theta);
      vel [0] = x_theta_0_dot * cos(theta);
      vel [1] = x_theta_0_dot * sin(theta);
      vel [2] = x_theta_0_dot * (-9.81 * (x_theta - x_theta_0)*inv_x_theta_0_dot_sq + tan(alpha));
      return vel;
    }
    

    
    
  } //   namespace rbprm
} // namespace hpp

