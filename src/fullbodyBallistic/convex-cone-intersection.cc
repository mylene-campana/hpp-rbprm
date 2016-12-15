//
// Authors: Mylene Campana (2016)
//
// Math tools to compute a convex cone intersection
// The following functions may be integrated in HPP later
// (once all cases and formula are completed).
//
// This file is part of hpp-rbprm
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
// hpp-rbprm. If not, see
// <http://www.gnu.org/licenses/>.

#include "hpp/rbprm/fullbodyBallistic/convex-cone-intersection.hh"
#include <hpp/util/debug.hh>

namespace convexCone
{

  bool force_closure (Cones cones, const value_type mu) {
    hppDout (info, "Force-closure case detection ----");
    const std::size_t N = cones.size ();
    const value_type phi = atan (mu);
    hppDout (info, "phi= " << phi);
    if (N == 2) {
      	hppDout (info, "2 cones");
      const value_type scalProd = cones [0][0]*cones [1][0] + cones [0][1]*cones [1][1] + cones [0][2]*cones [1][2];
      hppDout (info, "scalProd = " << scalProd);
      const value_type psi = acos (-scalProd);
      hppDout (info, "psi = " << psi);
      hppDout (info, "2*phi = " << 2*phi);
      //hppDout (info, "fabs(cos(2*phi)) = " << fabs(cos(2*phi)));
      //if (scalProd < 0 && fabs(cos(2*phi)) < fabs (scalProd)) {
      if (psi > -2*phi && psi < 2*phi) {
	hppDout (info, "2 cones and force closure case DETECTED");
	return true;
      }
    }
    if (N > 2) {
      hppDout (info, "3+ cones, force closure computation NOT IMPLEMENTED");
    }
    return false;
  }

  // -----------------------------------------------------------------------

  bool intersection_existence_zone (const fcl::Vec3f n, const fcl::Vec3f n2,
				    const value_type theta,const value_type mu){
    bool result = false;
    const value_type Ctheta = cos(theta); // n_theta = [Stheta,-Ctheta,0]]
    const value_type Stheta = sin(theta);
    const value_type scalProd1 = n [0]*Stheta - n [1]*Ctheta;
    const value_type scalProd2 = n2 [0]*Stheta - n2 [1]*Ctheta;
    
    // test intersection of each cone with plane_theta  and intersection of the cones' middle-zone with plane_theta:
    //result = cone_plane_inter (n, theta, mu) || cone_plane_inter (n2, theta, mu)|| ((scalProd1 > 0 && scalProd2 < 0) || (scalProd1 < 0 && scalProd2 > 0));
    // test only intersection of the cones' middle-zone with plane_theta:
    result = ((scalProd1 > 0 && scalProd2 < 0) || (scalProd1 < 0 && scalProd2 > 0));
    hppDout (info, "intersection_zone existence = " << result);
    return result;
  }

  // -----------------------------------------------------------------------

  bool is_parallel_plane_theta (const fcl::Vec3f n, const fcl::Vec3f n2,
				const value_type theta) {
    const value_type tx = n[0] - n2[0];
    const value_type ty = n[1] - n2[1];
    hppDout (info, "tx = " << tx << ", ty = " << ty);
    if (theta != M_PI/2 && theta != -M_PI/2) {
      const value_type Ttheta = tan(theta);
      if (ty - tx*Ttheta == 0) {
	hppDout (info, "plane_theta is parallel to t_dir");
	return true;
      }
    } else { // theta = pi/2 or theta = -pi/2
      if (tx == 0) {
	hppDout (info, "plane_theta is parallel to t_dir");
	return true;
      }
    }
    hppDout (info, "plane_theta is NOT parallel to t_dir");
    return false;
  }

  // -----------------------------------------------------------------------

  fcl::Vec3f compute_bissectrice_dir (const fcl::Vec3f Mplus,
				      const fcl::Vec3f Mminus) {
    fcl::Vec3f CC_dir;
    value_type Mp_norm, Mm_norm, xm, ym, zm, xp, yp, zp;
    // normalize:
    Mp_norm = Mplus.norm ();
    Mm_norm = Mminus.norm ();
    xp = Mplus [0] / Mp_norm; yp = Mplus [1] / Mp_norm;
    zp = Mplus [2] / Mp_norm; xm = Mminus [0] / Mm_norm;
    ym = Mminus [1] / Mm_norm; zm = Mminus [2] / Mm_norm;
    CC_dir [0] = 0.5*(xp + xm); 
    CC_dir [1] = 0.5*(yp + ym);
    CC_dir [2] = 0.5*(zp + zm);
    CC_dir.normalize ();
    return CC_dir;
  }


  // -----------------------------------------------------------------------

  void compute_M_points (fcl::Vec3f* Mplus, fcl::Vec3f* Mminus,
			 const fcl::Vec3f n, const fcl::Vec3f n2,
			 const value_type theta, const value_type mu) {
    
    // constants
    const value_type nx = n[0]; const value_type ny = n[1];
    const value_type nz = n[2]; const value_type n2x = n2[0];
    const value_type n2y = n2[1]; const value_type n2z = n2[2];
    const value_type Cyz = ny*n2z - nz*n2y; // cross-product n x n2
    const value_type Czx = nz*n2x - nx*n2z;
    const value_type Cxy = nx*n2y - ny*n2x;
    const value_type K = nx*n2x + ny*n2y + nz*n2z;
    const value_type mu12 = 1 + mu*mu;

    const value_type tx = nx - n2x;
    const value_type ty = ny - n2y;
    const value_type tz = nz - n2z;

    value_type xPp, yPp, zPp, xPm, yPm, zPm; // P_plus and P_minus coordinates
    value_type discr;
    value_type xMp, yMp, zMp, xMm, yMm, zMm; // M_plus and M_minus coordinates

    /// If t is not parallel to plan_theta, allow to compute P points
    
    hppDout (info, "compute P points");
    hppDout (info, "Cyz = " << Cyz << ", Czx = " << Czx << ", Cxy = " << Cxy);
    value_type A, B, C;
    if (nz != 0) {
      if (Cyz != 0) {
	// case I
	hppDout (info, "case I");
	A = 1 + Czx*Czx/(Cyz*Cyz) + pow((Czx*ny + nx*Cyz)/(nz*Cyz), 2);
	B = -2*(Czx*(nz*K-n2z)/(Cyz*Cyz) + (Czx*ny + Cyz*nx)*(Cyz + ny*nz*K - ny*n2z)/(nz*nz*Cyz*Cyz));
	C = pow((nz*K - n2z)/Cyz,2) + pow((Cyz + ny*nz*K - ny*n2z)/(nz*Cyz),2) - mu12;
	discr = B*B - 4*A*C;
	hppDout (info, "discr= " << discr);
	if (discr < 0)
	  hppDout (error, "case I, problem with negative discr");
	xPp = 0.5*(-B + sqrt(discr))/A;
	xPm = 0.5*(-B - sqrt(discr))/A;
	yPp = Czx/Cyz*xPp - (nz*K-n2z)/Cyz;
	yPm = Czx/Cyz*xPm - (nz*K-n2z)/Cyz;
	zPp = (1 - xPp*nx - yPp*ny)/nz;
	zPm = (1 - xPm*nx - yPm*ny)/nz;
      } else { // (Cyz == 0)
	if (Czx != 0) {
	  // case II
	  hppDout (info, "case II");
	  xPp = (nz*K - n2z)/Czx; xPm = xPp;
	  A = 1 + ny*ny/(nz*nz);
	  B = -2*ny*(1 - xPp*nx)/(nz*nz);
	  C = pow((1 - xPp*nx)/nz,2) + xPp*xPp - mu12;
	  discr = B*B - 4*A*C;
	  hppDout (info, "discr= " << discr);
	  if (discr < 0)
	    hppDout (error, "case II, problem with negative discr");
	  yPp = 0.5*(-B + sqrt(discr))/A;
	  yPm = 0.5*(-B - sqrt(discr))/A;
	  zPp = (1 - xPp*nx - yPp*ny)/nz; zPm = (1 - xPm*nx - yPm*ny)/nz;
	} else { // (Czx == 0)
	  // impossible
	  hppDout (error, "should not happen, n // n2" );
	}
      }
    } else { // nz == 0
      if (ny != 0) {
	if (n2z != 0) {
	  // case III
	  hppDout (info, "case III");
	  A = 1 + pow (nx/ny,2) + pow (Cxy/(ny*n2z),2);
	  B = -2*nx/(ny*ny) + 2*Cxy*(ny*K-n2y)/pow(ny*n2z,2);
	  C = 1/(ny*ny) + pow((ny*K-n2y)/(ny*n2z),2) - mu12;
	  discr = B*B - 4*A*C;
	  hppDout (info, "discr= " << discr);
	  if (discr < 0)
	    hppDout (error, "case III, problem with negative discr");
	  xPp = 0.5*(-B + sqrt(discr))/A;
	  yPp = (1-nx*xPp)/ny;
	  zPp = (Cxy*xPp + ny*K - n2y)/(ny*n2z);
	  xPm = 0.5*(-B - sqrt(discr))/A;
	  yPm = (1-nx*xPm)/ny;
	  zPm = (Cxy*xPm + ny*K - n2y)/(ny*n2z);
	} else { // n2z == 0
	  if (Cxy != 0) {
	    // case IV
	    hppDout (info, "case IV");
	    xPp = (n2y - ny*K)/Cxy;
	    yPp = (1-nx*xPp)/ny;
	    zPp = sqrt(mu*mu - xPp*xPp - yPp*yPp);
	    xPm = xPp;
	    yPm = yPp;
	    zPm = -sqrt(mu*mu - xPm*xPm - yPm*yPm);
	  } else { // Cxy == 0
	    // case V
	    hppDout (info, "case V");
	    xPp = nx;
	    yPp = (1-nx*xPp)/ny;
	    zPp = sqrt(mu*mu-nx*nx-pow(nx,4)/(ny*ny));
	    xPm = xPp;
	    yPm = yPp;
	    zPm = -sqrt(mu*mu-nx*nx-pow(nx,4)/(ny*ny));
	  }
	}
      } else { // ny == 0
	if (n2z != 0) {
	  // case VI
	  hppDout (info, "case VI");
	  xPp = 1;
	  yPp = mu*sqrt(n2z*n2z/(n2z*n2z+n2y*n2y));
	  zPp = -n2y*yPp/n2z;
	  xPm = 1;
	  yPm = -mu*sqrt(n2z*n2z/(n2z*n2z+n2y*n2y));
	  zPm = -n2y*n2z/yPm;
	} else { // n2z == 0
	  if (n2y != 0) {
	    // case VII
	    hppDout (info, "case VII");
	    xPp = 0; yPp = 0; zPp = mu;
	    xPm = 0; yPm = 0; zPm = -mu;
	  } else { // n2y == 0
	    // impossible
	    hppDout (error, "should not happen, n // n2");
	  }
	}
      }
    }
    hppDout (info, "P_plus: " << xPp << ", " << yPp << ", " << zPp);
    hppDout (info, "P_minus: " << xPm << ", " << yPm << ", " << zPm);

    /// Compute M points from P
    value_type alpha;
    if (theta != M_PI/2 && theta != -M_PI/2) {
      const value_type Ttheta = tan(theta);
      const value_type tytxtTheta = ty - tx*Ttheta; // should not be = 0 here
      hppDout (info, "tz= " << tz);
      hppDout (info, "tytxtTheta= " << tytxtTheta);
      // case I_m and II_m
      alpha = -(yPp - Ttheta*xPp)/tytxtTheta;
      xMp = alpha*tx + xPp; yMp = alpha*ty + yPp; zMp = alpha*tz + zPp;
      alpha = -(yPm - Ttheta*xPm)/tytxtTheta;
      xMm = alpha*tx + xPm; yMm = alpha*ty + yPm; zMm = alpha*tz + zPm;
    } else { // theta = pi/2 or theta = -pi/2
      hppDout (info, "tx= " << tx);
      // Case III_m
      alpha = -xPp/tx;
      xMp = 0; yMp = alpha*ty + yPp; zMp = alpha*tz + zPp;
      alpha = -xPm/tx;
      xMm = 0; yMm = alpha*ty + yPp; zMm = alpha*tz + zPm;
    }
    hppDout (info, "M_plus: " << xMp << ", " << yMp << ", " << zMp);
    hppDout (info, "M_minus: " << xMm << ", " << yMm << ", " << zMm);
    (*Mplus) [0] = xMp; (*Mplus) [1] = yMp; (*Mplus) [2] = zMp;
    (*Mminus) [0] = xMm; (*Mminus) [1] = yMm; (*Mminus) [2] = zMm;

  } // compute_M_points

  // -----------------------------------------------------------------------

  bool cone_plane_inter (const fcl::Vec3f n, const value_type theta,
			 const value_type mu) {
    hppDout (info, "IROS version of cone-plane intersection existence ----");
    const value_type U = n [0]; // n_x
    const value_type V = n [1]; // n_y
    const value_type W = n [2]; // n_z
    hppDout (info, "U= " << U << ", V= " << V << ", W= " << W);
    const value_type denomK = U*U + V*V - W*W*mu*mu;
    const bool tanThetaDef = theta != M_PI /2 && theta != -M_PI /2;
    hppDout (info, "tanThetaDef = " << tanThetaDef);

    if (denomK > -1e-6 && denomK < 1e-6) { // denomK (or 'A') = 0
      hppDout (info, "denomK = 0 case");
      if (tanThetaDef) {
	hppDout (info, "cone-plane intersection OK");
	return true;
      } else { // theta = +-pi/2
	if (V != 0) {
	  hppDout (info, "cone-plane intersection OK");
	  return true;
	} else { // V = 0
	  hppDout (info, "cone-plane intersection is a line");
	  return true; // accepted for convex-cone
	}
      }
    } // if denomK = 0
      
    if (tanThetaDef) {
      const value_type tantheta = tan(theta);
      const value_type discr = (U*U+W*W)*mu*mu - V*V - U*U*tantheta*tantheta + (V*V + W*W)*mu*mu*tantheta*tantheta + 2*(1+mu*mu)*U*V*tantheta;
      hppDout (info, "discr: " << discr);
      if (discr < 0) {
	hppDout (info, "cone-plane intersection empty");
	return false;
      }
      /*if (discr < 5e-2) {
	hppDout (info, "cone-plane intersection too small");
	return false;
	}*/ // accepted here (line)
      hppDout (info, "cone-plane intersection OK");
      return true;
    }
    else { // theta = +-pi/2
      const value_type discr =  -U*U+(V*V + W*W)*(mu*mu);
      hppDout (info, "discr: " << discr);
      if (discr < 0) {
	hppDout (info, "cone-plane intersection empty");
	return false;
      }
      /*if (discr < 5e-2) {
	hppDout (info, "cone-plane intersection too small");
	return false;
      }*/ // accepted here (line)
      hppDout (info, "cone-plane intersection OK");
      return true;
    }
  }

  // -----------------------------------------------------------------------

  bool cone_circle_plane_inter (fcl::Vec3f* Mplus, fcl::Vec3f* Mminus,
				const fcl::Vec3f n, const value_type theta,
				const value_type mu) {
    hppDout (info, "cone-circle plane intersection existence -------");
    const value_type nx = n[0]; const value_type ny = n[1];
    const value_type nz = n[2];
    const value_type mu12 = 1 + mu*mu;
    value_type x, y;
    const bool tanThetaDef = theta != M_PI /2 && theta != -M_PI /2;
    value_type A, B, C;
    hppDout (info, "tanThetaDef = " << tanThetaDef);
    if (tanThetaDef) {
      const value_type tanTheta = tan(theta);
      const value_type nxytTheta = nx + tanTheta*ny;
      if (nz != 0) {
	A = 1 + tanTheta*tanTheta + pow((nxytTheta/nz),2);
	B = -2*nxytTheta/(nz*nz);
	C = 1/(nz*nz) - mu12;
	const value_type delta = B*B - 4*A*C;
	hppDout (info, "delta = " << delta);
	if (delta > 0) {
	  // Case I_MC
	  hppDout (info, "case I_MC");
	  x = 0.5*(-B + sqrt(delta))/A;
	  (*Mplus) [0] = x;
	  (*Mplus) [1] = x * tanTheta;
	  (*Mplus) [2] = (1 - nx*x - ny*x*tanTheta)/nz;
	  x = 0.5*(-B - sqrt(delta))/A;
	  (*Mminus) [0] = x;
	  (*Mminus) [1] = x * tanTheta;
	  (*Mminus) [2] = (1 - nx*x - ny*x*tanTheta)/nz;
	  hppDout (info, "cone-circle-plane intersection OK");
	  hppDout (info, "M_plus: " << *Mplus);
	  hppDout (info, "M_minus: " << *Mminus);
	  return true;
	}
	if (delta == 0) {
	  // Case II_MC
	  hppDout (info, "case II_MC");
	  x = -B*0.5/A;
	  (*Mplus) [0] = x;
	  (*Mplus) [1] = x * tanTheta;
	  (*Mplus) [2] = (1 - nx*x - ny*x*tanTheta)/nz;
	  *Mminus = *Mplus;
	  hppDout (info, "cone-circle-plane intersection OK");
	  hppDout (info, "M_plus: " << *Mplus);
	  return true;
	}
	if (delta < 0) {
	  hppDout (info, "intersection cone-circle-plane is empty");
	  return false;
	}
      } // nz != 0
      else { // nz == 0
	if (ny != 0) {
	  if (nxytTheta != 0) {
	    // Case V_MC
	    hppDout (info, "case V_MC");
	    x = 1/nxytTheta;
	    (*Mplus) [0] = x;
	    (*Mplus) [1] = (1 - nx*x)/ny;
	    (*Mplus) [2] = sqrt(mu12 - (1+tanTheta*tanTheta)*x*x);
	    (*Mminus) [0] = x;
	    (*Mminus) [1] = (1 - nx*x)/ny;
	    (*Mminus) [2] = -sqrt(mu12 - (1+tanTheta*tanTheta)*x*x);
	    hppDout (info, "M_plus: " << *Mplus);
	    hppDout (info, "M_minus: " << *Mminus);
	    hppDout (info, "cone-circle-plane intersection OK");
	    return true;
	  } else { // nxytTheta == 0
	    hppDout (info, "cone-circle-plane intersection is empty");
	    return false;
	  }
	} // ny != 0
	else { // ny == 0
	  // Case VI_MC
	  hppDout (info, "case VI_MC");
	  (*Mplus) [0] = 1;
	  (*Mplus) [1] = tanTheta;
	  (*Mplus) [2] = sqrt(mu12 - 1 - tanTheta*tanTheta);
	  (*Mminus) [0] = 1;
	  (*Mminus) [1] = tanTheta;
	  (*Mminus) [2] = -sqrt(mu12 - 1 - tanTheta*tanTheta);
	  hppDout (info, "M_plus: " << *Mplus);
	  hppDout (info, "M_minus: " << *Mminus);
	  hppDout (info, "cone-circle-plane intersection OK");
	  return true;
	}
      }
    } else { // theta = +- pi/2
      if (nz != 0) {
	A = 1 + pow(ny/nz,2);
	B = -2*ny/(nz*nz);
	C = 1/(nz*nz) - mu12;
	const value_type delta = B*B - 4*A*C;
	hppDout (info, "delta = " << delta);
	if (delta > 0) {
	  // Case III_MC
	  hppDout (info, "case III_MC");
	  y = 0.5*(-B + sqrt(delta))/A;
	  (*Mplus) [0] = 0;
	  (*Mplus) [1] = y;
	  (*Mplus) [2] = (1 - ny*y)/nz;
	  y = 0.5*(-B - sqrt(delta))/A;
	  (*Mminus) [0] = 0;
	  (*Mminus) [1] = y;
	  (*Mminus) [2] = (1 - ny*y)/nz;
	  hppDout (info, "M_plus: " << *Mplus);
	  hppDout (info, "M_minus: " << *Mminus);
	  hppDout (info, "cone-circle-plane intersection OK");
	  return true;
	}
	if (delta == 0) { // Case IV_MC
	  hppDout (info, "case IV_MC");
	  y = -B*0.5/A;
	  (*Mplus) [0] = 0;
	  (*Mplus) [1] = y;
	  (*Mplus) [2] = (1 - ny*y)/nz;
	  (*Mminus) = *Mplus;
	  hppDout (info, "M_plus: " << *Mplus);
	  hppDout (info, "cone-circle-plane intersection OK");
	  return true;
	}
	if (delta < 0) {
	  hppDout (info, "intersection cone-circle-plane is empty");
	  return false;
	}
      } // if nz != 0
      else { // nz == 0
	if (ny != 0) {
	  // Case VII_MC
	  hppDout (info, "case VII_MC");
	  (*Mplus) [0] = 0;
	  (*Mplus) [1] = 1/ny;
	  (*Mplus) [2] = sqrt(mu12 - 1/(ny*ny));
	  (*Mminus) [0] = 0;
	  (*Mminus) [1] = 1/ny;
	  (*Mminus) [2] = -sqrt(mu12 - 1/(ny*ny));
	  hppDout (info, "M_plus: " << *Mplus);
	  hppDout (info, "M_minus: " << *Mminus);
	  hppDout (info, "cone-circle-plane intersection OK");
	  return true;
	} // ny != 0
	else {
	  hppDout (info, "cone-circle-plane intersection is empty");
	  return false;
	}
      } // nz == 0
    } // theta = +- pi/2
    return false;
  }

  // -----------------------------------------------------------------------

  // computeMaxRange
  value_type compute_angle (std::vector<fcl::Vec3f> M_vec,
			    const value_type theta, fcl::Vec3f* Mplus_border,
			    fcl::Vec3f* Mminus_border) {
    value_type phi_cc = 0;
    const std::size_t N_points = M_vec.size ();
    const value_type cTheta = cos(theta);
    const value_type sTheta = sin(theta);
    value_type z_i, x_theta_i,  psi_i;
    fcl::Vec3f Mplus = M_vec [0];// border max angle
    fcl::Vec3f Mminus = M_vec [0]; // border min angle
    value_type minPsi = std::numeric_limits <value_type>::infinity();
    value_type maxPsi = 0;
    for (std::size_t i = 0; i < N_points; i++) {
      fcl::Vec3f M_i = M_vec [i];
      hppDout (info, "M_i" << M_i); // for VIEWER-PLOT
      z_i = M_i[2];
      x_theta_i = M_i[0]*cTheta + M_i[1]*sTheta;
      psi_i = atan2 (z_i, x_theta_i);
      if (psi_i > maxPsi) {
	Mplus = M_i;
	maxPsi = psi_i;
      }
      if (psi_i < minPsi) {
	Mminus = M_i;
	minPsi = psi_i;
      }
    }
    hppDout (info, "END of intersection points ---"); // for VIEWER-PLOT
    const value_type crossScal = Mplus[0]*Mminus[0] + Mplus[1]*Mminus[1] + Mplus[2]*Mminus[2];
    const value_type scal_p = Mplus[0]*Mplus[0] + Mplus[1]*Mplus[1] + Mplus[2]*Mplus[2];
    const value_type scal_n = Mminus[0]*Mminus[0] + Mminus[1]*Mminus[1] + Mminus[2]*Mminus[2];
    hppDout (info, "crossScal = " << crossScal);
    hppDout (info, "scal_n = " << scal_n );
    hppDout (info, "scal_p = " << scal_p );
    phi_cc = 0.5*acos (crossScal/(scal_p*scal_n));
    *Mplus_border = Mplus;
    *Mminus_border = Mminus;
    hppDout (info, "phi_cc = " << phi_cc);
    hppDout (info, "M_border = " << Mplus); hppDout (info, "M_border = " << Mminus); // same line for VIEWER-PLOT
    hppDout (info, "END of border points ---"); // for VIEWER-PLOT
    return phi_cc;
  }

  // -----------------------------------------------------------------------

  // General algorithm
  vector_t compute_convex_cone_inter (const Cones cones,
				      const value_type theta,
				      const value_type mu) {
    value_type angle = 0;
    vector_t result (4);
    const std::size_t Ncones = cones.size ();
    hppDout (info, "Number of contact cones= " << Ncones);
    std::vector<fcl::Vec3f> M_vec;
    if (force_closure (cones, mu)) {
      // force-closure case, return max angle as if no sliding constraint
      angle = M_PI;
      result [0] = angle;
      result [1] = 0;
      result [2] = 0;
      result [3] = 0;
      return result;
    }

    if (Ncones == 1) {
      fcl::Vec3f n = cones [0];
      fcl::Vec3f Mplus; fcl::Vec3f Mminus;
      if (cone_plane_inter (n, theta, mu)) {
	cone_circle_plane_inter (&Mplus, &Mminus, n, theta, mu);
	M_vec.push_back (Mplus);
	if (Mplus != Mminus)
	  M_vec.push_back (Mminus);
      }
    }
    else {
      for (std::size_t i = 0; i < Ncones; i++) {
	for (std::size_t j = i + 1; j < Ncones; j++) {
	  hppDout (info, "i= " << i << ", j= " << j << " **********");
	  fcl::Vec3f n = cones [i];
	  fcl::Vec3f n2 = cones [j];
	  fcl::Vec3f Mplus; fcl::Vec3f Mminus;
	  if (cone_plane_inter (n, theta, mu)) {
	    cone_circle_plane_inter (&Mplus, &Mminus, n, theta, mu);
	    M_vec.push_back (Mplus);
	    if (Mplus != Mminus)
	      M_vec.push_back (Mminus);
	  }
	  if (cone_plane_inter (n2, theta, mu)) {
	    cone_circle_plane_inter (&Mplus, &Mminus, n2, theta, mu);
	    M_vec.push_back (Mplus);
	    if (Mplus != Mminus)
	      M_vec.push_back (Mminus);
	  }
	  if (!is_parallel_plane_theta (cones [i], cones [j], theta)) {
	    //intersection_existence_zone (n, n2, theta, mu); // NOT USED
	    // If t is not parallel to plan_theta, allow to compute P->M points
	    compute_M_points (&Mplus, &Mminus, n, n2, theta, mu);
	    M_vec.push_back (Mplus);
	    M_vec.push_back (Mminus);
	  }
	}
      }
    }
    hppDout (info, "M_vec.size () = " << M_vec.size ());
    if (!M_vec.size ()) {
      // plane_theta - convex-cone intersection is empty
      hppDout (info, "CC-zone plane_theta intersection EMPTY -----------");
      result [0] = 0;
      result [1] = 0;
      result [2] = 0;
      result [3] = 0;
      return result;
    }
    fcl::Vec3f Mplus; fcl::Vec3f Mminus; // borders
    angle = compute_angle (M_vec, theta, &Mplus, &Mminus);
    const fcl::Vec3f CC_dir = compute_bissectrice_dir (Mplus, Mminus);
    hppDout (info, "CC_dir = " << CC_dir);
    result [0] = angle;
    result [1] = CC_dir [0]; // direction of 2D-convex-cone
    result [2] = CC_dir [1];
    result [3] = CC_dir [2];
    return result;
  }

} // namespace convexCone
