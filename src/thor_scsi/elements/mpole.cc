#include <thor_scsi/elements/mpole.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/elements_enums.h>
#include <thor_scsi/elements/utils.h>
#include <tps/math_pass.h>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

/* ==========================================================================
   Support functions
   required by other elements too ?

   If so: should go to element_helpers.h
   ========================================================================== */

namespace thor_scsi::elements {
	template<typename T>
	void edge_focus(const tsc::ConfigType &conf, const double irho, const double phi /* in degrees */,
			    const double gap, ss_vect<T> &ps);
	template<typename T>
	void p_rot(const tsc::ConfigType &conf, double phi, ss_vect<T> &ps);

	template<typename T>
	void bend_fringe(const tsc::ConfigType &conf, const double hb, ss_vect<T> &ps);

	template<typename T>
	void quad_fringe(const tsc::ConfigType &conf, const double b2, ss_vect<T> &ps);

}
/**
 * Todo:
 *     check phi
 */
template<typename T>
void tse::edge_focus(const tsc::ConfigType &conf, const double irho, const double phi /* in degrees */,
		    const double gap, ss_vect<T> &ps)
{
  ps[px_] += irho*tan(degtorad(phi))*ps[x_];
  if (!conf.dip_edge_fudge) {
    // Remark: Leads to a diverging Taylor map (see SSC-141).
    // ps[py_] -=
    //   irho*tan(degtorad(phi)-get_psi(irho, phi, gap))
    //   *ps[y_]/(1e0+ps[delta_]);
    // Leading order correction.
    ps[py_] -=
      irho*tan(degtorad(phi)-get_psi(irho, phi, gap))*ps[y_]*(1e0-ps[delta_]);
  } else
    ps[py_] -= irho*tan(degtorad(phi)-get_psi(irho, phi, gap))*ps[y_];
}

// how does it differ from PRotTransform? Should it be implemented there?
template<typename T>
void tse::p_rot(const tsc::ConfigType &conf, double phi, ss_vect<T> &ps)
{
  T  c, s, t, pz, val; // p
  ss_vect<T> ps1;

  c = cos(degtorad(phi));
  s = sin(degtorad(phi));
  t = tan(degtorad(phi));
  pz = get_p_s(conf, ps);

  if (!conf.H_exact && !conf.Cart_Bend) {
     ps[px_] = s*pz + c*ps[px_];
  } else {
    // ps1 = ps; p = c*pz - s*ps1[px_];
    // px[x_]   = ps1[x_]*pz/p; px[px_] = s*pz + c*ps1[px_];
    // px[y_]  += ps1[x_]*ps1[py_]*s/p;
    // px[ct_] += (1e0+ps1[delta_])*ps1[x_]*s/p;

    ps1 = ps; val = 1e0 - ps1[px_]*t/pz;
    ps[x_]  = ps1[x_]/(c*val);
    ps[px_] = ps1[px_]*c + s*pz;
    ps[y_]  = ps1[y_] + t*ps1[x_]*ps1[py_]/(pz*val);
    ps[ct_] = ps1[ct_] + ps1[x_]*(1e0+ps1[delta_])*t/(pz*val);
  }
}

template<typename T>
void tse::bend_fringe(const tsc::ConfigType &conf, const double hb, ss_vect<T> &ps)
{
  T          coeff, u, pz, pz2, pz3;
  ss_vect<T> ps1;

  coeff = -hb/2e0;
  ps1 = ps;
  pz = get_p_s(conf, ps);
  pz2 = sqr(pz); pz3 = pz*pz2;
  u = 1e0 + 4e0*coeff*ps1[px_]*ps1[y_]*ps1[py_]/pz3;
  if (u >= 0e0) {
    ps[y_]  = 2e0*ps1[y_]/(1e0+sqrt(u));
    ps[x_]  = ps1[x_] - coeff*sqr(ps[y_])*(pz2+sqr(ps1[px_]))/pz3;
    ps[py_] = ps1[py_] + 2e0*coeff*ps1[px_]*ps[y_]/pz;
    ps[ct_] = ps1[ct_] - coeff*ps1[px_]*sqr(ps[y_])*(1e0+ps1[delta_])/pz3;
  } else {
    printf("bend_fringe: *** Speed of light exceeded!\n");
    ps[x_] = NAN; ps[px_] = NAN; ps[y_] = NAN; ps[py_] = NAN;
    ps[delta_] = NAN; ps[ct_] = NAN;
  }
}


template<typename T>
void tse::quad_fringe(const tsc::ConfigType &conf, const double b2, ss_vect<T> &ps)
{
  T u, p_s;

  u = b2/(12e0*(1e0+ps[delta_])); p_s = u/(1e0+ps[delta_]);
  ps[py_] /= 1e0 - 3e0*u*sqr(ps[y_]); ps[y_] -= u*tse::cube(ps[y_]);
  if (conf.Cavity_on) ps[ct_] -= p_s*tse::cube(ps[y_])*ps[py_];
  ps[px_] /= 1e0 + 3e0*u*sqr(ps[x_]);
  if (conf.Cavity_on) ps[ct_] += p_s*tse::cube(ps[x_])*ps[px_];
  ps[x_] += u*tse::cube(ps[x_]); u = u*3e0; p_s = p_s*3e0;
  ps[y_] = exp(-u*sqr(ps[x_]))*ps[y_]; ps[py_] = exp(u*sqr(ps[x_]))*ps[py_];
  ps[px_] += 2e0*u*ps[x_]*ps[y_]*ps[py_];
  if (conf.Cavity_on) ps[ct_] -= p_s*sqr(ps[x_])*ps[y_]*ps[py_];
  ps[x_] = exp(u*sqr(ps[y_]))*ps[x_]; ps[px_] = exp(-u*sqr(ps[y_]))*ps[px_];
  ps[py_] -= 2e0*u*ps[y_]*ps[x_]*ps[px_];
  if (conf.Cavity_on) ps[ct_] += p_s*sqr(ps[y_])*ps[x_]*ps[px_];
}



/* ==========================================================================
   End support functions
   ========================================================================== */


/*
 * Review if mpole should be simplified....
 */
template<typename T>
void tse::MpoleType::_localPass(tsc::ConfigType &conf, ss_vect<T> &ps)
{
  int          seg = 0, i;
  double       dL = 0e0, dL1 = 0e0, dL2 = 0e0,
               dkL1 = 0e0, dkL2 = 0e0, h_ref = 0e0;

  // now handled by LocalCoordinateElement
  //GtoL(ps, dS, dT, Pc0, Pc1, Ps1);

  if (conf.emittance && !conf.Cavity_on) {
    // Needs A^-1.
    curly_dH_x = 0e0;
    for (i = 0; i <= 5; i++)
      /* Synchrotron integrals */ dI[i] = 0e0;
  }

  switch (Pmethod) {

  case Meth_Fourth:
    if (conf.mat_meth && (Porder <= Quad)) {
      ps = mat_pass(M_elem, ps);

      // if (conf.emittance && !conf.Cavity_on) &&
      // 	(PL != 0e0) && (Pirho != 0e0)) get_dI_eta_5(this);
    } else {
      // Fringe fields.
      if (conf.quad_fringe && (PB[Quad+HOMmax] != 0e0))
	tse::quad_fringe(conf, PB[Quad+HOMmax], ps);
      if (!conf.Cart_Bend) {
	if (Pirho != 0e0)
	  tse::edge_focus(conf, Pirho, PTx1, Pgap, ps);
      } else {
	tse::p_rot(conf, PTx1, ps); bend_fringe(conf, Pirho, ps);
      }

      if (Pthick == thick) {
	if (!conf.Cart_Bend) {
	  // Polar coordinates.
	  h_ref = Pirho; dL = PL/PN;
	} else {
	  // Cartesian coordinates.
	  h_ref = 0e0;
	  if (Pirho == 0e0)
	    dL = PL/PN;
	  else
	    dL = 2e0/Pirho*sin(PL*Pirho/2e0)/PN;
	}


#warning "Enable dl1, dL2, dkL1, dkL2 calculation"
#if 0
	dL1 = c_1*dL; dL2 = c_2*dL; dkL1 = d_1*dL; dkL2 = d_2*dL;
#endif

	for (seg = 1; seg <= PN; seg++) {
	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }

	  drift_pass(conf, dL1, ps);
	  tse::thin_kick(conf, Porder, *this->intp, dkL1, Pirho, h_ref, ps);
	  drift_pass(conf, dL2, ps);
	  tse::thin_kick(conf, Porder, *this->intp, dkL2, Pirho, h_ref, ps);

	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += 4e0*is_tps<tps>::get_curly_H(ps);
	    dI[4] += 4e0*is_tps<tps>::get_dI_eta(ps);
	  }

	  drift_pass(conf, dL2, ps);
	  thin_kick(conf, Porder, *this->intp, dkL1, Pirho, h_ref, ps);
	  drift_pass(conf, dL1, ps);

	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }
	}

	if (conf.emittance && !conf.Cavity_on) {
	  // Needs A^-1.
	  curly_dH_x /= 6e0*PN;
	  dI[1] += PL*is_tps<tps>::get_dI_eta(ps)*Pirho;
	  dI[2] += PL*sqr(Pirho);
	  dI[3] += PL*fabs(tse::cube(Pirho));
	  dI[4] *=
	    PL*Pirho*(sqr(Pirho)+2e0*PBpar[Quad+HOMmax])
	    /(6e0*PN);
	  dI[5] += PL*fabs(tse::cube(Pirho))*curly_dH_x;
	}
      } else
	thin_kick(conf, Porder, *this->intp, 1e0, 0e0, 0e0, ps);

      // Fringe fields.
      if (!conf.Cart_Bend) {
	if (Pirho != 0e0)
	  tse::edge_focus(conf, Pirho, PTx2, Pgap, ps);
      } else {
	bend_fringe(conf, -Pirho, ps); p_rot(conf, PTx2, ps);
      }
#warning "Deactivated quad fringe"
#if 0
      if (conf.quad_fringe && (PB[Quad+HOMmax] != 0e0))
	quad_fringe(conf, -PB[Quad+HOMmax], ps);
#endif
    }
    break;

  default:
    std::cerr <<  "Mpole_Pass: Method not supported " << this->name
	      <<  " method " << Pmethod <<  std::endl;
    throw ts::NotImplemented();
    break;
  }

  // now handled by LocalCoordinateElement
  // LtoG(ps, dS, dT, Pc0, Pc1, Ps1);
}
