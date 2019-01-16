cimport MixedLayerModel

cdef class NetCDFIO_Stats:
    cdef:
        object root_grp
        object profiles_grp
        object ts_grp

        str stats_path
        str path_plus_file

        public double last_output_time
        public double frequency

    cpdef initialize(self, dict namelist, MixedLayerModel.MixedLayerModel mlm)
    cpdef open_files(self)
    cpdef close_files(self)
    cpdef setup_stats_file(self, MixedLayerModel.MixedLayerModel mlm)
    cpdef add_profile(self, var_name)
    cpdef add_ts(self, var_name)
    cpdef write_profile(self, var_name, double[:] data)
    cpdef write_ts(self, var_name, double data)
    cpdef write_simulation_time(self, double t)
