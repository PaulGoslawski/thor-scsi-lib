import os
import numpy as np

from thor_scsi.lib import (
    ss_vect_double,
    ss_vect_tps,
    ConfigType,
    RadiationDelegate,
    RadiationDelegateKick,
    ObservedState
)

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics

from thor_scsi.utils.phase_space_vector import map2numpy

from thor_scsi.utils.phase_space_vector import ss_vect_tps2ps_jac, \
    array2ss_vect_tps

X_ = 0
Y_ = 1
Z_ = 2


c0      = 2.99792458e8
m_e     = 0.51099906e6
# q_e = 1.60217663e-19
q_e     = 1.602e-19
mu_0    = 4e0*np.pi*1e-7
eps_0   = 1e0/(c0**2*mu_0)
r_e     = q_e/(4e0*np.pi*eps_0*m_e)
C_u     = 55e0/(24e0*np.sqrt(3e0))
C_gamma = 4e0*np.pi*r_e/(3e0*(1e-9*m_e)**3)
#cl_rad  = C_gamma*(conf.Energy)**3/(2e0*np.pi)
h_bar   = 6.58211899e-16
C_q     = 3e0*C_u*h_bar*c0/(4e0*m_e)
#q_fluct = C_q*C_gamma/(np.pi*sqr(1e-9*m_e))*pow(conf.Energy, 5e0)


def comp_hor_emit(lat, conf, prt):
    conf.Energy = 2.5 # [GeV].
    C = np.sum([elem.getLength() for elem in lat])

    C_q_scl = 1e18*C_q/(m_e)**2
    E_0     = 1e9*conf.Energy
    T_0     = C/c0

    conf.emittance = False

    M_tmp = linear_optics.compute_map(lat, conf)
    M = map2numpy(M_tmp)
    print(M)
    A_j = linear_optics.find_phase_space_origin(M)
    A_tmp = np.zeros([7, 7], dtype=np.float)
    A_tmp[:6, :6] = A_j
    A = array2ss_vect_tps(A_tmp)
    print(A)

    conf.Cavity_on = False
    conf.radiation = True
    conf.emittance = True

    lat.propagate(conf, A)

    type_name = "Marker"
    ps_zero = ss_vect_double()
    rad_del = [RadiationDelegate() for elem in
               lat.elementsWithNameType(type_name)]
    for a_del, elem in zip(rad_del, lat.elementsWithNameType(type_name)):
        elem.setRadiationDelegate(a_del)
        # Just use that that the marker knows who is calling him
        a_del.view(elem, ps_zero, ObservedState.start, 0)

    # Not used any more better to clean up the name space
    del ps_zero


    # Radiation delegate for field kick
    type_name = "Bending"
    radiators = [elem for elem in lat.elementsWithNameType(type_name)]
    rad_del_kick = [
        # RK()
        RadiationDelegateKick()
        for elem in radiators
    ]
    for a_del, elem in zip(rad_del_kick, radiators):
        # Should be set from accelerator
        a_del.setEnergy(1e9*conf.Energy)
        elem.setRadiationDelegate(a_del)

    for e in radiators:
        print(repr(e))
        break

    lat.propagate(conf, A)

    # Inspect curly_H in
    for a_del in rad_del:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x:5f}"
        print(txt)

    for a_del in rad_del_kick:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        dI = a_del.getSynchrotronIntegralsIncrements()
        D_rad = a_del.getDiffusionCoefficientsIncrements()
        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x: 10.6e}"
        txt += "    dI " + ",".join(["{: 10.6e}".format(v) for v in dI])
        txt += "   "
        txt += "    D_rad" + ",".join(["{: 10.6e}".format(v) for v in D_rad])
        txt += "   "
        print(txt)

    I_tmp = np.sum([a_del.getSynchrotronIntegralsIncrements() for a_del in
                    rad_del_kick], axis=0)
    I = np.zeros(6, float)
    I[1:] = I_tmp
    D_rad = np.sum([a_del.getDiffusionCoefficientsIncrements() for a_del in
                    rad_del_kick], axis=0)

    print(f"{I=}")
    print(f"{D_rad=}")

    U_0         = 1e9*C_gamma*conf.Energy**4*I[2]/(2e0*np.pi)
    eps_x       = C_q_scl*conf.Energy**2*I[5]/(I[2]-I[4])
    sigma_delta = np.sqrt(C_q_scl*conf.Energy**2*I[3]/(2e0*I[2]+I[4]))
    J = np.zeros(3, float)
    J[X_]       = 1e0 - I[4]/I[2]
    J[Z_]       = 2e0 + I[4]/I[2]
    J[Y_]       = 4e0 - J[X_] - J[Z_]

    for k in range(0, 3):
#        tau[k] = 4e0*np.pi*T_0/(C_gamma*cube(1e-9*E_0)*J[k]*I[2])

        if True:
            print("\n  I[1..5]:")
            for k in range(1, 6):
                print(" {:10.3e}".format(I[k]))
            print("\n")

            print("\n  U_0   [keV]    = {:5.1f}\n".format(1e-3*U_0))
            print("  eps_x [nm.rad] = {:6.4f}\n".format(1e9*eps_x))
            print("  sigma_delta    = {:9.3e}\n".format(sigma_delta))
            print("  J              = [{:5.3f}, {:5.3f}, {:5.3f}]\n".
                  format(J[X_], J[Y_], J[Z_]))
            print("  tau   [msec]   = [{:9.3e}, {:9.3e}, {:9.3e}]\n".
                  format(1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]))

    return [eps_x, sigma_delta, U_0, J, tau, I]


t_dir  = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

lat  = accelerator_from_config(t_file)
conf = ConfigType()


#[eps_x, sigma_delta, U_0, J, tau, I] = comp_hor_emit(lat, conf, True)
comp_hor_emit(lat, conf, True)
