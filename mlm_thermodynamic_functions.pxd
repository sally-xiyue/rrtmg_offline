cdef get_pressure(double [:] z, double p_surface, double rho0)

cdef double saturation_vapor_pressure(double temp_val) nogil

cdef qv_star(double p0, double qt, double pv)

cdef thetal_isdac(double p_, double t_, double ql_, double qt_)

cdef get_theta_v(double p_, double t_, double ql_, double qt_)

cdef sat_adjst(double p_, double thetal_, double qt_)

cdef qv_star_rh(double p0, double rh, double pv)

cdef qv_unsat(double p0, double pv)

cdef double get_alpha(double p0, double T, double qt, double qv)

cdef get_ls_subsidence(double [:] z, double zi_i, double div_frac)