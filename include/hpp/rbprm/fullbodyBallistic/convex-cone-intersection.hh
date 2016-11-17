
#ifndef CONVEX_CONE_INTERSECTION
#define CONVEX_CONE_INTERSECTION


#include <iostream>
#include <Eigen/Dense>
#include <Eigen/src/Core/util/Macros.h>
#include <vector>

#include <hpp/model/configuration.hh>
#include <hpp/fcl/collision_data.h>
#include <hpp/core/fwd.hh>

namespace convexCone
{
  using hpp::model::displayConfig;
  using hpp::core::value_type;
  using hpp::core::vector_t;
  typedef std::vector<fcl::Vec3f> Cones; // list of cone normals

  /// Return true if this is a force-closure case
  bool force_closure (const Cones cones, const value_type mu);

  /// Compute and return true if intersection exists
  bool intersection_existence_zone (const fcl::Vec3f n, const fcl::Vec3f n2,
				    const value_type theta,const value_type mu);

  /// Return true if t=n1-n2 is parallel to plane_theta
  bool is_parallel_plane_theta (const fcl::Vec3f n, const fcl::Vec3f n2,
			       const value_type theta);

  /// Compute the normalized bissectrice direction between the directions
  /// [O-Mplus] and [O-Mminus]
  fcl::Vec3f compute_bissectrice_dir (const fcl::Vec3f Mplus,
				      const fcl::Vec3f Mminus);

  /// Compute points resulting from the intersection of the convex-cone
  /// and the plane_theta.
  void compute_M_points (fcl::Vec3f* Mplus, fcl::Vec3f* Mminus,
			 const fcl::Vec3f n, const fcl::Vec3f n2,
			 const value_type theta, const value_type mu);

  /// Compute if the intersection between a cone and plane_theta is empty.
  /// Return false if intersection is a point (= "empty").
  /// IROS VERSION
  bool cone_plane_inter (const fcl::Vec3f n, const value_type theta,
			 const value_type mu);

  /// Compute the intersection between a cone-circle and plane_theta.
  /// TODO !!!!!
  /// SEE IF RESULT OF EXISTENCE IF COHERENT WITH IROS
  bool cone_circle_plane_inter (fcl::Vec3f* Mplus, fcl::Vec3f* Mminus,
				const fcl::Vec3f n, const value_type theta,
				const value_type mu);

  /// Compute opening angle from M points
  // computeMaxRange
  value_type compute_angle (std::vector<fcl::Vec3f> M_vec,
			    const value_type theta);

  /// Full algorithm computing the angle of the intersection between the
  /// convex cone and the plane_theta
  vector_t compute_convex_cone_inter (const Cones cones,
					const value_type theta,
					const value_type mu);

} //namespace convexCone

#endif //CONVEX_CONE_INTERSECTION
