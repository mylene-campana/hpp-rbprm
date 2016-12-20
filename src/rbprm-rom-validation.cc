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
    (const model::DevicePtr_t& robot, const std::vector<std::string>& affFilters)
    {
      RbPrmRomValidation* ptr = new RbPrmRomValidation (robot, affFilters);
      return RbPrmRomValidationPtr_t (ptr);
    }

    RbPrmRomValidation::RbPrmRomValidation (const model::DevicePtr_t& robot
					    ,const std::vector<std::string>& affFilters)
      : hpp::core::CollisionValidation(robot)
      , filter_(affFilters)
      , unusedReport_(new CollisionValidationReport) {}

    bool RbPrmRomValidation::validate (const Configuration_t& config)
    {
      return validate(config, unusedReport_);
    }

    bool RbPrmRomValidation::validate (const Configuration_t& config,
				       ValidationReportPtr_t& validationReport)
    {
      ValidationReportPtr_t romReport;
      bool collision = !hpp::core::CollisionValidation::validate(config, romReport);
      CollisionValidationReportPtr_t reportCast = boost::dynamic_pointer_cast<CollisionValidationReport>(romReport);
      RbprmValidationReportPtr_t rbprmReport =boost::dynamic_pointer_cast<RbprmValidationReport>(validationReport);
      if(!rbprmReport)
        hppDout(notice,"problem cast rbprm-validation-report instance");

      if(collision)
        {
	  if(rbprmReport){  // if the report is a correct rbprm report, we add the rom information
	    rbprmReport->ROMReports.insert(std::make_pair(robot_->name(),boost::dynamic_pointer_cast<CollisionValidationReport>(romReport)));
	    rbprmReport->ROMFilters.insert(std::make_pair(robot_->name(),collision));
	  }else{
	    validationReport = romReport;
	  }
	}
      rbprmReport->ROMFilters.insert(std::make_pair(robot_->name(),collision)); // here report = 0;;
      return collision;
    }

  }// namespace rbprm
}// namespace hpp
