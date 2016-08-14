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

#include <hpp/rbprm/interpolation/limb-rrt.hh>
#include <hpp/rbprm/interpolation/time-constraint-utils.hh>
#include <hpp/core/bi-rrt-planner.hh>

namespace hpp {
using namespace core;
  namespace rbprm {
  namespace interpolation {
    void SetLimbRRTConstraints::operator ()(LimbRRTHelper& helper, const State& from, const State& to)
    {
        helper.SetContactConstraints(from, to);
    }

    core::PathPtr_t limbRRT(RbPrmFullBodyPtr_t fullbody, core::ProblemPtr_t referenceProblem,
                 const rbprm::CIT_State &startState, const rbprm::CIT_State &endState, const std::size_t numOptimizations)
    {
        LimbRRTShooterFactory factory;
        return interpolateStates<LimbRRTHelper, LimbRRTShooterFactory, CIT_State >
                (fullbody, referenceProblem, factory, startState, endState, numOptimizations);
    }

    core::PathPtr_t limbRRTFromPath(RbPrmFullBodyPtr_t fullbody, core::ProblemPtr_t referenceProblem, const PathPtr_t refPath,
                         const CIT_StateFrame &startState, const CIT_StateFrame &endState, const  std::size_t numOptimizations)
    {
        LimbRRTShooterFactory factory;
        return interpolateStatesFromPath<LimbRRTHelper, LimbRRTShooterFactory>
                (fullbody, referenceProblem, factory, refPath, startState, endState, numOptimizations);
    }

  }// namespace interpolation
  }// namespace rbprm
}// namespace hpp