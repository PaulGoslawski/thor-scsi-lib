"""
"""
import thor_scsi.lib as tslib
from .accelerator import instrument_with_radiators
import numpy as np
from scipy.constants import c
from dataclasses import dataclass

from thor_scsi.utils.closed_orbit import compute_closed_orbit
from thor_scsi.utils.linear_optics import compute_M_diag
from thor_scsi.utils.phase_space_vector import vec_mat2ss_vect_tps
from thor_scsi.utils.output import mat2txt, vec2txt

import logging
logger = logging.getLogger("thor-scsi")


X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]

x_, px_, y_, py_, ct_, delta_ = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta
]


@dataclass
class RadiationResult:
    #: the different taus
    relaxation_constants: np.ndarray
    emittance: np.ndarray
    #: for consistency checks
    fractional_tunes: np.ndarray


def compute_circ(acc):
    return np.sum([elem.getLength() for elem in acc])


def compute_diffusion_coefficients(rad_del_kicks):
    dD_rad = \
        np.array([rk.getDiffusionCoefficientsIncrements()
                  for rk in rad_del_kicks])
    D_rad = np.sum(dD_rad, axis=0)
    return D_rad


def compute_rad_prop(acc, calc_config, x0, dE, alpha_rad, D_rad):
    prt_deb = not False

    dof = 3
    J = np.zeros(dof)
    tau = np.zeros(dof)
    eps = np.zeros(dof)
    C = compute_circ(acc)
    print("\nC = ", C)
    U_0 = calc_config.Energy*dE
    for k in range(dof):
        J[k] = 2e0*(1e0+x0[delta_])*alpha_rad[k]/dE
        tau[k] = -C/(c*alpha_rad[k])
        eps[k] = -D_rad[k]/(2e0*alpha_rad[k])

    if prt_deb:
        print("\nE [GeV]     = {:3.1f}".format(1e-9*calc_config.Energy))
        print("U0 [keV]    = {:3.1f}".format(1e-3*U_0))
        print("eps         = {:12.6e} {:12.6e} {:12.6e}".
              format(eps[X_], eps[Y_], eps[Z_]))
        print("tau [msec]  = {:8.6f} {:8.6f} {:8.6f}".
              format(1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]))
        print("J           = {:8.6f} {:8.6f} {:8.6f}".
              format(J[X_], J[Y_], J[Z_]))
        print("alpha_rad   = {:13.6e} {:13.6e} {:13.6e}".
                format(alpha_rad[X_], alpha_rad[Y_], alpha_rad[Z_]))

        print("D_rad       = {:12.6e} {:12.6e} {:12.6e}".
              format(D_rad[X_], D_rad[Y_], D_rad[Z_]))



def compute_radiation(
    acc: tslib.Accelerator,
    calc_config: tslib.ConfigType,
    E,
    eps
):

    dof = 3

    calc_config.Energy    = E
    calc_config.radiation = True
    calc_config.emittance = False
    calc_config.Cavity_on = True

    # Install radiators that radiation is calculated
    rad_del_kicks = instrument_with_radiators(acc, energy=E)

    r = compute_closed_orbit(acc, calc_config, delta=0e0, eps=eps)
    #print("M:\n" + mat2txt(map2numpy(r.one_turn_map)))
    M = r.one_turn_map[:6, :6]

    logger.info(
        "M:\n" + mat2txt(M)
        + "\nx0 =", vec2txt(r.x0)
    )

    calc_config.dE = 0e0
    ps = r.x0
    acc.propagate(calc_config, ps)
    dE = calc_config.dE

    A, A_inv, alpha_rad = compute_M_diag(dof, M)

    calc_config.emittance = True

    A_cpy = vec_mat2ss_vect_tps(r.x0, A)
    acc.propagate(calc_config, A_cpy)

    D_rad = compute_diffusion_coefficients(rad_del_kicks)

    compute_rad_prop(acc, calc_config, r.x0, dE, alpha_rad, D_rad)


# def calculate_radiation(
#     acc: tslib.Accelerator,
#     *,
#     energy,
#     calc_config: tslib.ConfigType = None,
#     install_radiators: bool = True,
#     dof=3
# ):
#     """

#     Todo:
#         Rename function
#         Inspect if radiators are installed?
#         check in calc config if radation is requested
#         in this case install the radiators

#     1. calculate fix point and Poincarè Map M with damped system
#        (i.e. radiation on and cavity on (without dispersion in a second case)
#     2. diagonalise M = A $\Gamma$ A$^{-1}$
#     3. eigenvalues:

#         - complex part: tunes,
#         - real part: damping times  (refer equation)

#        use eigen values of symplectic matrix to identify the planes

#     4. propagate A, thin kick will create diffusion coeffs (don't forget to
#        zero them before calculation starts (sum it up afterwards

#     """

#     if install_radiators:
#         logger.debug("Installing radiators")
#         # keep variable as long as you need to do calculations
#         rad_del = instrument_with_radiators(acc, energy=energy)

#     if calc_config is None:
#         raise AssertionError
#         calc_config = tslib.ConfigType()
#         calc_config.radiation = True
#         # is this used anywhere?
#         calc_config.emittance = False
#         calc_config.Cavity_on = True

#     calc_config.Energy = energy
#     logger.debug(
#         f"calc_config radiation { calc_config.radiation} emmittance {calc_config.emittance} Cavity on {calc_config.Cavity_on}"
#     )


#     # Diffusion coefficients
#     acc.propagate(calc_config, Ap)
#     print("Ap after calc")
#     print(Ap)

#     r = RadiationResult(relaxation_constants=w.real, fractional_tunes=w.imag)
