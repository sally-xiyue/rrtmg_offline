cimport Radiation
cimport timestepping
from NetCDFIO cimport NetCDFIO_Stats

cdef class MixedLayerModel:
    cdef:
        double gamma_thetal
        double thetal_i
        double thetal_ft
        double rh_ft
        double div_frac
        double dz
        Py_ssize_t nz
        double [:] z
        double p_surface
        double zi_i
        double t_surface
        double qt_surface
        double rho0
        double rh_i
        double [:] pressure
        double [:] temperature
        double [:] qv
        double [:] z_interface
        double [:] pressure_i
        double [:] ql
        double [:] qi
        double [:] thetal
        double [:] qt
        double [:] values
        double [:] tendencies
        double efficiency

    # cdef:
    #     # Radiation.Radiation radiation
    #     int count
    #     double radiation_frequency
    #     double next_radiation_calculate

    cpdef initialize(self, NetCDFIO_Stats NS)
    cpdef update(self, timestepping.TimeStepping TS, Radiation.Radiation Ra)
    cpdef stats_io(self, NetCDFIO_Stats NS)