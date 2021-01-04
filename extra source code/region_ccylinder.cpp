/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "region_ccylinder.h"
#include <cmath>
#include <cstring>
#include "update.h"
#include "domain.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20
enum{CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

RegCcylinder::RegCcylinder(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  options(narg-5,&arg[5]);

  // check open face settings

  if (openflag && (open_faces[3] || open_faces[4] || open_faces[5]))
    error->all(FLERR,"Invalid region cylinder open setting");

  Radius = xscale * force->numeric(FLERR, arg[2]);
  theta = force->numeric(FLERR, arg[3]);
  radius = xscale * force->numeric(FLERR, arg[4]);

  n1x = 0;
  n1y = 1;
  n2x = std::sin(theta);
  n2y = -std::cos(theta);
  half_flag = (theta > M_PI);

  

  // error check

  if (radius <= 0.0 || Radius <= 0.0) error->all(FLERR,"Illegal region cylinder command");

  // extent of ccylinder
 

  if (interior) {
    bboxflag = 1;
	double rSum = radius + Radius;
      extent_xlo = -rSum;
      extent_xhi = rSum;
      extent_ylo = -rSum;
      extent_yhi = rSum;
	  extent_zlo = -radius;
	  extent_zhi = radius;
    
  } else bboxflag = 0;

  // particle could be close to cylinder surface and 2 ends
  // particle can only touch surface and 1 end

  cmax = 3;
  contact = new Contact[cmax];
  if (interior) tmax = 2;
  else tmax = 1;
}

/* ---------------------------------------------------------------------- */

RegCcylinder::~RegCcylinder()
{
  delete [] contact;
}

/* ---------------------------------------------------------------------- */

void RegCcylinder::init()
{
  Region::init();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegCcylinder::inside(double x, double y, double z)
{
  if (half_flag)
  {
	  if (n1x*x + n1y * y < 0.0 && n2x*x + n2y * y < 0.0)
		  return 0;
  }
  else
  {
	  if (!(n1x*x + n1y * y > 0.0 && n2x*x + n2y * y > 0.0))
		  return 0;
  }

  double dist = sqrt(x * x + y * y) - Radius;
  if (dist * dist + z * z <= radius * radius)
	  return 1;
  else
	  return 0;
}

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of cylinder
   can be one contact for each of 3 cylinder surfaces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on cylinder to x
   special case: no contact with curved surf if x is on center axis
------------------------------------------------------------------------- */

int RegCcylinder::surface_interior(double *x, double cutoff)
{
  
	
  double del0 = x[0], del1 = x[1], del2 = x[2];
  double dot1 = n1x * del0 + n1y * del1;
  double dot2 = n2x * del0 + n2y * del1;

	  // x is exterior to torus

  if (half_flag)
  {
	  if (dot1 < 0.0 && dot2 < 0.0)
		  return 0;
  }
  else
  {
	  if (!(dot1 > 0.0 && dot2 > 0.0))
		  return 0;
  }

  int n = 0;
  double vir_x = sqrt(del0 * del0 + del1 * del1);
  double dist = vir_x - Radius;
  double r = sqrt(dist * dist + del2 * del2);

  // x is exterior to torus
  if (r > radius)
	  return 0;

 
  double r_diff = radius - r; 
  
  // x is on the surface or interior
  if (r_diff < cutoff && !open_faces[2])
  {
	  double p = 1.0 - radius / r;
	  double del_dist = dist * p;
	  double temp = dist * p;
	  contact[n].r = r_diff;
	  contact[n].delx = temp / vir_x * del0;
	  contact[n].dely = temp / vir_x * del1;
	  contact[n].delz = del2 * p;
	  contact[n].radius = -2.0*radius;
	  contact[n].iwall = 2;
	  contact[n].varflag = 0;
	  n++;
  }
  if (dot1 < cutoff && del0 > 0 && !open_faces[0]) {
	  contact[n].r = dot1;
	  contact[n].delx = contact[n].delz = 0;
	  contact[n].dely = dot1;
	  contact[n].radius = 0;
	  contact[n].iwall = 0;
	  contact[n].varflag = 0;
	  n++;
  }
  if (dot2 < cutoff && del0 * n2y - del1 * n1x < 0 && !open_faces[1]) {
	  contact[n].r = dot2;
	  contact[n].delx = dot2 * n2x;
	  contact[n].dely = dot2 * n2y;
	  contact[n].delz = 0.0;
	  contact[n].radius = 0;
	  contact[n].iwall = 1;
	  contact[n].varflag = 0;
	  n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of cylinder
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on cylinder to x
------------------------------------------------------------------------- */

int RegCcylinder::surface_exterior(double *x, double cutoff)
{
	error->one(FLERR, "surface_exterior function called");
	return 0;
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegCcylinder::shape_update()
{
      error->one(FLERR,"shape_update is called");
 
}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegCcylinder::variable_check()
{
	error->one(FLERR, "variable_check is called");

}


/* ----------------------------------------------------------------------
   Set values needed to calculate velocity due to shape changes.
   These values do not depend on the contact, so this function is
   called once per timestep by fix/wall/gran/region.

------------------------------------------------------------------------- */

void RegCcylinder::set_velocity_shape()
{
	error->one(FLERR, "set_velocity_shape is called");

}



/* ----------------------------------------------------------------------
   add velocity due to shape change to wall velocity
------------------------------------------------------------------------- */

void RegCcylinder::velocity_contact_shape(double *vwall, double *xc)
{
	error->one(FLERR, "velocity_contact_shape is called");
}

