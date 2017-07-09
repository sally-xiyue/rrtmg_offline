#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: initializedcheck=False
#cython: cdivision=True

import netCDF4 as nc
import os
import shutil
import numpy as np
cimport numpy as np
import cython
cimport MixedLayerModel

cdef class NetCDFIO_Stats:

    def __init__(self):
        self.root_grp = None
        self.profiles_grp = None
        self.ts_grp = None

        return

    # @cython.wraparound(True)
    cpdef initialize(self, dict namelist, MixedLayerModel.MixedLayerModel mlm):

        self.last_output_time = 0.0
        # self.uuid = str(namelist['meta']['uuid'])
        self.frequency = namelist['stats_io']['frequency']

        self.stats_path = str(namelist['stats_io']['output_root'])
        try:
            os.mkdir(self.stats_path)
        except:
            pass

        self.path_plus_file = str(self.stats_path + 'Stats.' + namelist['meta']['simname'] + '.nc')

        self.setup_stats_file(mlm)
        return

    cpdef open_files(self):
        self.root_grp = nc.Dataset(self.path_plus_file, 'r+', format='NETCDF4')
        self.profiles_grp = self.root_grp.groups['profiles']
        self.ts_grp = self.root_grp.groups['timeseries']

        return

    cpdef close_files(self):
        self.root_grp.close()
        return

    cpdef setup_stats_file(self, MixedLayerModel.MixedLayerModel mlm):

        root_grp = nc.Dataset(self.path_plus_file, 'w', format='NETCDF4')

        profile_grp = root_grp.createGroup('profiles')
        profile_grp.createDimension('z', mlm.nz)
        profile_grp.createDimension('t', None)
        z = profile_grp.createVariable('z', 'f8', ('z'))
        z[:] = np.array(mlm.z)
        profile_grp.createVariable('t', 'f8', ('t'))
        del z

        ts_grp = root_grp.createGroup('timeseries')
        ts_grp.createDimension('t', None)
        ts_grp.createVariable('t', 'f8', ('t'))

        root_grp.close()
        return

    cpdef add_profile(self, var_name):
        root_grp = nc.Dataset(self.path_plus_file, 'r+', format='NETCDF4')
        profile_grp = root_grp.groups['profiles']
        new_var = profile_grp.createVariable(var_name, 'f8', ('t', 'z'))
        root_grp.close()
        return

    cpdef add_ts(self, var_name):
        root_grp = nc.Dataset(self.path_plus_file, 'r+', format='NETCDF4')
        ts_grp = root_grp.groups['timeseries']
        new_var = ts_grp.createVariable(var_name, 'f8', ('t',))
        root_grp.close()
        return

    cpdef write_profile(self, var_name, double[:] data):
        var = self.profiles_grp.variables[var_name]
        var[-1, :] = np.array(data)
        return

    @cython.wraparound(True)
    cpdef write_ts(self, var_name, double data):
        var = self.ts_grp.variables[var_name]
        var[-1] = data
        return

    cpdef write_simulation_time(self, double t):
        profile_t = self.profiles_grp.variables['t']
        profile_t[profile_t.shape[0]] = t

        ts_t = self.ts_grp.variables['t']
        ts_t[ts_t.shape[0]] = t
        return


