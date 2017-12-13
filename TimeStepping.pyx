#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: initializedcheck=False
#cython: cdivision=True

import numpy as np
cimport numpy as np
cimport MixedLayerModel

cdef class TimeStepping:

    def __init__(self):
        return

    cpdef initialize(self, namelist, MixedLayerModel.MixedLayerModel mlm):

        try:
            self.dt = namelist['time_stepping']['dt_initial']
        except:
            self.dt = 600.0

        try:
            self.dt_max = namelist['time_stepping']['dt_max']
        except:
            self.dt_max = 600.0

        try:
            self.t = namelist['time_stepping']['t']
        except:
            self.t = 0.0

        try:
            self.t_max = namelist['time_stepping']['t_max']
        except:
            self.t_max = 12 * 3600.0 # default 12 hours of simulation

        self.initialize_second(mlm)

        return


    cpdef update_second(self, MixedLayerModel.MixedLayerModel mlm):

        cdef:
            Py_ssize_t i

        # with nogil:
        if self.rk_step == 0:
            for i in xrange(3):
                self.value_copies[i] = mlm.values[i]
                mlm.values[i] += mlm.tendencies[i]*self.dt
                mlm.tendencies[i] = 0.0
        else:
            for i in xrange(3):
                mlm.values[i] = 0.5 * (self.value_copies[i] + mlm.values[i] + mlm.tendencies[i] * self.dt)
                mlm.tendencies[i] = 0.0
            self.t += self.dt

        if self.t + self.dt > self.t_max:
            self.dt = self.t_max - self.t

        return

    cdef void initialize_second(self, MixedLayerModel.MixedLayerModel mlm):

        self.rk_step = 0
        self.n_rk_steps = 2

        self.value_copies = np.zeros((3,), dtype=np.double, order='c')

        return