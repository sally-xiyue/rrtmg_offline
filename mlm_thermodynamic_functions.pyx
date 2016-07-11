import numpy as np
cimport numpy as np
# import matplotlib.pyplot as plt
include 'parameters.pxi'
from libc.math cimport exp, cos

# constants
cdef double cp = 1004.0
# g = 9.80665
# Rd = 287.1
# Rv = 461.5
# eps_v = Rd / Rv
# eps_vi = 1 / eps_v
# p_tilde = 1.0e5
#
# f1 = 15.0
# f0 = 72.0
# kap = 170.0


cdef get_pressure(double [:] z, double p_surface, double rho0):
    # thetav = theta * (1.0 + 0.61 * qv - ql)

    cdef double [:] p_profile = np.zeros(z.shape[0])
    p_profile[0] = p_surface
    cdef double dz = z[1] - z[0]
    for i in np.arange(1, z.shape[0]):
        # rho_ = p_profile[i-1] / (Rd * (thetav * (p_profile[i-1]/p_tilde)**(Rd/cp)))
        p_profile[i] = p_profile[i - 1] - g * rho0 * dz

    return p_profile


cdef double saturation_vapor_pressure(double temp_val) nogil:
    # Twarm = 273.0
    # Tcold = 235.0
    # pow_n = 0.1

    # if temp_val >= Tcold and temp_val <= Twarm:
    #     lam = ((temp_val - Tcold)/(Twarm - Tcold))**pow_n
    # else:
    #     lam = 1.0

    cdef double es0 = 611.21
    cdef double t0 = 273.16

    cdef double a3l = 17.502
    cdef double a4l = 32.19
    cdef double ewl = es0 * exp(a3l * (temp_val - t0) / (temp_val - a4l))

    # a3i = 22.587
    # a4i = -0.7
    # ti = t0 - 23.0
    # ewi = es0 * np.exp(a3i * (temp_val - ti)/(temp_val - a4i))

    return ewl


cdef qv_star(double p0, double qt, double pv):
    return pv * (Rd / Rv) * (1 - qt) / (p0 - pv)


cdef thetal_isdac(double p_, double t_, double ql_, double qt_):
    cdef double rl_ = ql_ / (1 - qt_)

    return t_ * (p_tilde / p_) ** (Rd / cp) * exp(-(rl_ * 2.501e6) / (t_ * cp))


cdef get_theta_v(double p_, double t_, double ql_, double qt_):
    cdef double rl_ = ql_ / (1 - qt_)
    cdef double rv_ = (qt_ - ql_) / (1 - qt_)

    return t_ * (p_tilde / p_) ** (Rd / cp) * (1.0 + 0.61 * rv_ - rl_)


cdef sat_adjst(double p_, double thetal_, double qt_):
    """
    Use saturation adjustment scheme to compute temperature and ql given thetal and qt.
    :param p_: pressure [Pa]
    :param thetal_: liquid water potential temperature  [K]
    :param qt_:  total water specific humidity
    :return: t_2, ql_2
    """

    # Compute temperature
    cdef double t_1 = thetal_ * (p_ / p_tilde) ** (Rd / cp)
    # Compute saturation vapor pressure
    cdef double pv_star_1 = saturation_vapor_pressure(t_1)
    # Compute saturation mixing ratio
    cdef double qs_1 = qv_star(p_, qt_, pv_star_1)

    cdef:
        double ql_1
        double f_1
        double t_2
        double pv_star_2
        double qs_2
        double ql_2
        double t_n

    if qt_ <= qs_1:
        # If not saturated return temperature and ql = 0.0
        return t_1, 0.0
    else:
        ql_1 = qt_ - qs_1
        f_1 = thetal_ - thetal_isdac(p_, t_1, ql_1, qt_)
        t_2 = t_1 + 2.501e6 * ql_1 / cp
        pv_star_2 = saturation_vapor_pressure(t_2)
        qs_2 = qv_star(p_, qt_, pv_star_2)
        ql_2 = qt_ - qs_2

        while np.fabs(t_2 - t_1) >= 1e-9:
            pv_star_2 = saturation_vapor_pressure(t_2)
            qs_2 = qv_star(p_, qt_, pv_star_2)
            ql_2 = qt_ - qs_2
            f_2 = thetal_ - thetal_isdac(p_, t_2, ql_2, qt_)
            t_n = t_2 - f_2 * (t_2 - t_1) / (f_2 - f_1)
            t_1 = t_2
            t_2 = t_n
            f_1 = f_2

    return t_2, ql_2


cdef qv_star_rh(double p0, double rh, double pv):
    cdef double val = eps_v * pv / (p0 - pv) / (1 + rh * eps_v * pv / (p0 - pv))
    return val


cdef qv_unsat(double p0, double pv):
    cdef double val = 1.0 / (eps_vi * (p0 - pv) / pv + 1.0)
    return val

#
# def get_radiation_isdac(ql, z, rho0):
#     dz = np.mean(np.diff(z))
#     nz = ql.shape[0]
#     f_rad = np.zeros(nz + 1)
#     f_heat = np.zeros(nz)
#     q_1 = 0.0
#     f_rad[0] = f1 * np.exp(-q_1)
#     for k in xrange(1, nz + 1):
#         q_1 += kap * rho0 * ql[k - 1] * dz
#         f_rad[k] += f1 * np.exp(-q_1)
#
#     q_0 = 0.0
#     f_rad[nz] += f0 * np.exp(-q_0)
#     for k in xrange(nz - 1, -1, -1):
#         q_0 += kap * rho0 * ql[k] * dz
#         f_rad[k] += f0 * np.exp(-q_0)
#
#     for k in xrange(nz):
#         f_heat[k] = -(f_rad[k + 1] - f_rad[k]) / dz / rho0
#
#     return f_rad, f_heat


# large-scale subsidence
cdef get_ls_subsidence(double [:] z, double zi_i, double div_frac):
    cdef double divergence = 5.0e-6 * div_frac
    cdef double [:] w_ls = np.zeros_like(z)
    cdef Py_ssize_t k
    for k in xrange(len(z)):
        if z[k] <= zi_i:
            w_ls[k] = -divergence * zi_i * cos(pi * 0.5 * (zi_i - z[k]) / zi_i)
        else:
            w_ls[k] = -divergence * zi_i * cos(pi * 0.5 * (z[k] - zi_i) / (10000. - zi_i))

    return w_ls
