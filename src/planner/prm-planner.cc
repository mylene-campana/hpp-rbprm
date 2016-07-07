//
// Copyright (c) 2015-2016 CNRS
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
#include <hpp/core/node.hh>
#include <hpp/core/edge.hh>
#include <hpp/core/path.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/roadmap.hh>
#include <hpp/rbprm/planner/prm-planner.hh>

namespace hpp {
  namespace rbprm {
    using model::displayConfig;
    using core::value_type;
    using core::vector_t;
    using model::size_type;

    PrmPlanner::PrmPlanner (const core::Problem& problem):
      PathPlanner (problem),
      configurationShooter_ (problem.configurationShooter()),
      qProj_ (problem.robot ()->configSize ()),
      smParabola_(rbprm::SteeringMethodParabola::create((core::ProblemPtr_t(&problem)))),
      rbRoadmap_(core::RbprmRoadmap::create (problem.distance (),problem.robot())), roadmap_(boost::dynamic_pointer_cast<core::Roadmap>(rbRoadmap_))
    {
      hppDout(notice,"Constructor prm-planner");
    }

    PrmPlanner::PrmPlanner (const core::Problem& problem,
			    const core::RoadmapPtr_t& roadmap) :
      PathPlanner (problem, roadmap),
      configurationShooter_ (problem.configurationShooter()),
      qProj_ (problem.robot ()->configSize ()),
      smParabola_(rbprm::SteeringMethodParabola::create((core::ProblemPtr_t(&problem)))),
      rbRoadmap_(core::RbprmRoadmap::create (problem.distance (),problem.robot())), roadmap_(roadmap)
      //roadmap_(boost::dynamic_pointer_cast<core::Roadmap>(rbRoadmap_))
    {
      hppDout(notice,"Constructor prm-planner with Roadmap");
    }

    void PrmPlanner::oneStep ()
    {
      core::DevicePtr_t robot (problem ().robot ());
      core::PathPtr_t localPath;
      DelayedEdge_t fwdDelayedEdge, bwdDelayedEdge;
      DelayedEdges_t fwdDelayedEdges;
      
      // shoot a RB-valid random configuration using *normally* rbprm-shooter
      core::ConfigurationPtr_t q_rand;
      hppDout(notice,"# oneStep BEGIN");
      q_rand = configurationShooter_->shoot ();
      hppDout (info, "q_rand: " << displayConfig (*q_rand));

      // Add q_rand as a new node: here for the parabola, as the impact node
      core::NodePtr_t impactNode = roadmap ()->addNode (q_rand);
      impactNode->indexInRM (roadmap ()->nodeIndex_);
      roadmap ()->nodeIndex_++;

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
	    core::ConfigurationPtr_t qCC = (*n_it)->configuration ();

	    // Create forward local path from qCC to q_rand
	    localPath = (*smParabola_) (*qCC, *q_rand);

	    // if a forward path is returned, it is valid
	    if (localPath) {
	      // Save forward & backward delayed edges
	      fwdDelayedEdge = DelayedEdge_t (*n_it, impactNode, localPath);
	      fwdDelayedEdges.push_back (fwdDelayedEdge);
		
	      // Assuming that SM is symmetric (V0max = Vfmax)
	      // WARN: I had to reverse *n_it - impNode HERE
	      // to add edges consecutively to same vector fwdDelayedEdges
	      bwdDelayedEdge = DelayedEdge_t (impactNode, *n_it,
					      localPath->reverse ());
	      fwdDelayedEdges.push_back (bwdDelayedEdge);
	    } //if SM has returned a non-empty path
	  }//for nodes in cc
	}//avoid impactNode cc
      }//for cc in roadmap

      // Insert in roadmap all forward delayed edges (DE)
      bool avoidAddIdenticalEdge = true;
      for (DelayedEdges_t::const_iterator itEdge = fwdDelayedEdges.begin ();
	   itEdge != fwdDelayedEdges.end (); ++itEdge) {
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
	if (!avoidAddIdenticalEdge) {
	  roadmap ()->edgeIndex_++;
	  avoidAddIdenticalEdge = true;
	} else
	  avoidAddIdenticalEdge = false;
      }
    }

    void PrmPlanner::tryDirectPath ()
    {
      // call steering method here to build a direct conexion
      core::PathPtr_t path;
      std::vector<std::string> filter;
      core::NodePtr_t initNode = roadmap ()->initNode();
      for (core::Nodes_t::const_iterator itn = roadmap ()->goalNodes ().begin();itn != roadmap ()->goalNodes ().end (); ++itn) {
        core::ConfigurationPtr_t q1 ((initNode)->configuration ());
        core::ConfigurationPtr_t q2 ((*itn)->configuration ());
        assert (*q1 != *q2);
        path = (*smParabola_) (*q1, *q2);
        if (path) { // has no collision
	  hppDout(notice, "#### direct parabola path is valid !");
	  roadmap ()->addEdge (initNode, *itn, path);
	  roadmap ()->addEdge (*itn, initNode, path->reverse());
	} else {
	  hppDout(notice, "#### direct parabola path not valid !");
	}
      } //for qgoals
    }

  } // namespace core
} // namespace hpp
