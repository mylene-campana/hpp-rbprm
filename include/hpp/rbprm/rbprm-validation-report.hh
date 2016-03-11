//
// Copyright (c) 2016 CNRS
// Authors: Pierre Fernbach
//
// This file is part of hpp-rbprm
// hpp-rbprm is free software: you can redistribute it
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

#ifndef HPP_RBPRM_VALIDATION_REPORT_HH
# define HPP_RBPRM_VALIDATION_REPORT_HH

# include <hpp/core/validation-report.hh>
# include <hpp/core/collision-validation-report.hh>

namespace hpp {
  namespace core {
    /// \addtogroup validation
    /// \{

    /// Validate a configuration with respect to collision
    ///
    struct HPP_CORE_DLLAPI RbprmValidationReport : public CollisionValidationReport
    {
      /// Directing vector between collision point of geometries
      fcl::Vec3f outwardCOllisionDirection;
    }; // class RbprmValidationReport
    /// \}
  } // namespace core
} // namespace hpp

#endif // HPP_RBPRM__VALIDATION_REPORT_HH
