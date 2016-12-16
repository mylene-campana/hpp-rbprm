//
// Copyright (c) 2016 CNRS
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
#include <hpp/core/basic-configuration-shooter.hh>
#include <hpp/core/connected-component.hh>
#include <hpp/core/config-validations.hh>
#include <hpp/core/path-validation-report.hh>
#include <hpp/core/node.hh>
#include <hpp/core/edge.hh>
#include <hpp/core/path.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/roadmap.hh>
#include <hpp/rbprm/rbprm-path-validation.hh>
#include <hpp/rbprm/rbprm-validation-report.hh>
#include <hpp/rbprm/rbprm-state.hh>
#include <hpp/rbprm/fullbodyBallistic/ballistic-planner.hh>
#include <polytope/stability_margin.h>
#include "utils/algorithms.h"
#include <hpp/rbprm/fullbodyBallistic/parabola-library.hh>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::value_type;
    using core::vector_t;
    using model::size_type;

    BallisticPlanner::BallisticPlanner (const core::Problem& problem):
      PathPlanner (problem), problem_ (core::ProblemPtr_t(&problem)),
      configurationShooter_ (problem.configurationShooter()),
      smParabola_(rbprm::SteeringMethodParabola::create((core::ProblemPtr_t(&problem)))),
      roadmap_(boost::dynamic_pointer_cast<core::Roadmap>(core::RbprmRoadmap::create (problem.distance (),problem.robot()))),
      fullRobot_ (RbPrmFullBody::create(problem.robot ())),
      contactSize_ (core::vector_t(2))
    {
      hppDout(notice,"Constructor ballistic-planner");
    }

    BallisticPlanner::BallisticPlanner (const core::Problem& problem,
					const core::RoadmapPtr_t& roadmap) :
      PathPlanner (problem, roadmap), problem_ (core::ProblemPtr_t(&problem)),
      configurationShooter_ (problem.configurationShooter()),
      smParabola_(rbprm::SteeringMethodParabola::create((core::ProblemPtr_t(&problem)))),
      roadmap_(roadmap),
      fullRobot_ (RbPrmFullBody::create(problem.robot ())),
      contactSize_ (core::vector_t(2))
    {
      hppDout(notice,"Constructor ballistic-planner with Roadmap");
      hppDout(info,"contactSize_= " << contactSize_);
    }

    void BallisticPlanner::oneStep ()
    {
      core::DevicePtr_t robot (problem ().robot ());
      core::PathPtr_t fwdPath, bwdPath;
      DelayedEdge_t fwdDelayedEdge, bwdDelayedEdge;
      DelayedEdges_t delayedEdges; // store valid forward and backward
      //const size_type indexECS = robot->configSize () - robot->extraConfigSpace ().dimension ();

      // shoot a RB-valid random configuration using rbprm-shooter
      core::ConfigurationPtr_t q_rand;
      core::Configuration_t q_tmp;
      core::ValidationReportPtr_t report;
      bool valid = false;
      std::vector<fcl::Vec3f> contactCones;
      
      hppDout(notice,"# oneStep BEGIN");
      while (!valid) {
	q_rand = configurationShooter_->shoot ();
	valid = problem ().configValidations()->validate(*q_rand,report);
      }
      hppDout (info, "q_rand: " << displayConfig (*q_rand));

      // Add q_rand as a new node: here for the parabola, as the impact node
      core::NodePtr_t impactNode = roadmap ()->addNode (q_rand);
      impactNode->indexInRM (roadmap ()->nodeIndex_);
      roadmap ()->nodeIndex_++;
      core::RbprmNodePtr_t impactNodeRb = static_cast<core::RbprmNodePtr_t>(impactNode);
      if(impactNodeRb) hppDout(notice, "Node correctly casted to rbprmNode");
      else hppDout(error, "Impossible to cast node to rbprmNode");
      contactCones = library::computeContactCones (problem_, *q_rand);
      // Note: orientation not updated with 2D-CC direction since not computed
      impactNodeRb->contactCones (&contactCones);

      // try to connect the random configuration to each connected component
      // of the roadmap.
      for (core::ConnectedComponents_t::const_iterator itcc =
	     roadmap ()->connectedComponents ().begin ();
	   itcc != roadmap ()->connectedComponents ().end (); ++itcc) {
	core::ConnectedComponentPtr_t cc = *itcc;
	// except its own connected component of course
	if (cc != impactNode->connectedComponent ()) {

	  // iteration on each node of the current connected-component
	  for (core::NodeVector_t::const_iterator n_it = cc->nodes ().begin (); 
	       n_it != cc->nodes ().end (); ++n_it){
	    core::RbprmNodePtr_t n_itRb = static_cast<core::RbprmNodePtr_t>(*n_it);
	    if(!impactNodeRb) hppDout(error, "Impossible to cast node to rbprmNode");
	    core::ConfigurationPtr_t qCC = (*n_it)->configuration ();
	    hppDout (info, "qCC: " << displayConfig (*qCC));

	    // NOT TIME-OPTIMAL (I can store the cones in the roadmap... so I computed them each time)
	    contactCones = library::computeContactCones (problem_, *qCC);
	    n_itRb->contactCones (&contactCones);

	    // Create forward and backward paths
	    //fwdPath = (*smParabola_) (*qCC, *q_rand);
	    fwdPath = (*smParabola_) (n_itRb, impactNodeRb);
	    //bwdPath = (*smParabola_) (*q_rand, *qCC);
	    bwdPath = (*smParabola_) (impactNodeRb, n_itRb);
	    // if a path is returned (i.e. not null), then it is valid

	    if (fwdPath) {
	      hppDout (info, "forward path is valid");       
	      fwdDelayedEdge = DelayedEdge_t (*n_it, impactNode, fwdPath);
	      delayedEdges.push_back (fwdDelayedEdge);
	    }

	    if (bwdPath) {
	      hppDout (info, "backward path is valid");          
	      bwdDelayedEdge = DelayedEdge_t (impactNode, *n_it, bwdPath);
	      delayedEdges.push_back (bwdDelayedEdge);
	    }
	  }//for nodes in cc
	}//avoid impactNode cc
      }//for cc in roadmap

      // Insert in roadmap all forward delayed edges (DE)
      bool avoidAddIdenticalEdge = true;
      for (DelayedEdges_t::const_iterator itEdge = delayedEdges.begin ();
	   itEdge != delayedEdges.end (); ++itEdge) {
	const core::NodePtr_t& nodeDE = itEdge-> get <0> ();
	const core::NodePtr_t& node2DE = itEdge-> get <1> ();
	const core::PathPtr_t& pathDE = itEdge-> get <2> ();
	core::EdgePtr_t edge = roadmap ()->addEdge (nodeDE, node2DE, pathDE);
	hppDout(info, "connection between q1: " 
		<< displayConfig (*(nodeDE->configuration ()))
		<< "and q2: "
		<< displayConfig (*(node2DE->configuration ())));
	edge->indexInRM (roadmap ()->edgeIndex_);
	// assure that forward and backward edges have same edgeIndex
	// still kept for roadmap display (even if non-symmetric)
	if (!avoidAddIdenticalEdge) {
	  roadmap ()->edgeIndex_++;
	  avoidAddIdenticalEdge = true;
	} else
	  avoidAddIdenticalEdge = false;
      }
    }

    void BallisticPlanner::tryDirectPath ()
    {
        // call steering method here to build a direct conexion
      hppDout (info, "Try direct path -------------");
        core::PathPtr_t path;
        core::PathPtr_t fwdPath, bwdPath;
        const core::ConfigurationPtr_t q_init = roadmap ()->initNode()->configuration ();
        const core::RbprmNodePtr_t initNode = static_cast<core::RbprmNodePtr_t>(roadmap ()->initNode());
	if (!initNode) hppDout (error, "problem to cast RbprmNode");
	std::vector<fcl::Vec3f> cones;
	cones = library::computeContactCones (problem_, *q_init);
	initNode->contactCones (&cones);

      for (core::Nodes_t::const_iterator itn = roadmap ()->goalNodes ().begin();itn != roadmap ()->goalNodes ().end (); ++itn) {
	if (!*itn) hppDout (error, "problem to get goal node");
        const core::RbprmNodePtr_t goalNode = static_cast<core::RbprmNodePtr_t>(*itn);
	if (!goalNode) hppDout (error, "problem to cast RbprmNode");
	const core::ConfigurationPtr_t q_goal = (*itn)->configuration ();
	cones = library::computeContactCones (problem_, *q_goal);
	goalNode->contactCones (&cones);
        assert (*q_init != *q_goal);

        // Create forward and backward paths
	fwdPath = (*smParabola_) (initNode, goalNode);
	bwdPath = (*smParabola_) (goalNode, initNode);
        // if a path is returned (i.e. not null), then it is valid
        if (fwdPath) {
          hppDout (info, "forward path is valid");
          roadmap ()->addEdge(roadmap ()->initNode(), *itn, fwdPath);
        }
        if (bwdPath) {
          hppDout (info, "backward path is valid");
          roadmap ()->addEdge(*itn, roadmap ()->initNode(), bwdPath);
        }
      } //for qgoals
    }

  } // namespace core
} // namespace hpp
