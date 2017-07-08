cimport MixedLayerModel
cimport TimeStepping
cimport NetCDFIO
cimport Radiation
import time
import numpy as np


class Simulation:

    def __init__(self):
        return

    def initialize(self, namelist):
        self.mlm = MixedLayerModel.MixedLayerModel(namelist)
        self.TS = TimeStepping.TimeStepping()
        self.StatsIO = NetCDFIO.NetCDFIO_Stats()
        self.Ra = Radiation.Radiation(namelist)

        self.StatsIO.initialize(namelist, self.mlm)
        self.TS.initialize(namelist, self.mlm)
        self.Ra.initialize(self.mlm, self.StatsIO)
        self.mlm.initialize(self.StatsIO)

        return

    def run(self):

        self.mlm.update(self.TS, self.Ra)
        self.Ra.initialize_profiles(self.mlm)

        while self.TS.t < self.TS.t_max:
            time1 = time.time()
            for self.TS.rk_step in xrange(self.TS.n_rk_steps):
                self.Ra.update(self.mlm, self.TS)
                self.mlm.update(self.TS, self.Ra)
                self.TS.update_second(self.mlm)
                self.io()
            time2 = time.time()
            print('T = ' + str(self.TS.t) + ' dt = ' + str(self.TS.dt) + ' walltime = ' + str(time2 - time1))

        return

    def io(self):
        cdef:
            double stats_dt = 0.0

        if self.TS.t > 0.0 and self.TS.rk_step == self.TS.n_rk_steps - 1:
            stats_dt = self.StatsIO.last_output_time + self.StatsIO.frequency - self.TS.t

            if self.StatsIO.last_output_time + self.StatsIO.frequency == self.TS.t:
                self.StatsIO.last_output_time = self.TS.t
                self.StatsIO.open_files()
                self.StatsIO.write_simulation_time(self.TS.t)

                self.mlm.stats_io(self.StatsIO)
                self.Ra.stats_io(self.StatsIO)

                self.StatsIO.close_files()

        return

