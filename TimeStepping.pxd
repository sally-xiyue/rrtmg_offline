cimport MixedLayerModel

cdef class TimeStepping:
    cdef:
        public double dt
        public double dt_max
        public double t
        public double t_max
        double [:] value_copies
        public Py_ssize_t rk_step
        public Py_ssize_t n_rk_steps
        void initialize_second(self, MixedLayerModel.MixedLayerModel mlm)

    cpdef initialize(self, case_ic, MixedLayerModel.MixedLayerModel mlm)
    cpdef update_second(self, MixedLayerModel.MixedLayerModel mlm)