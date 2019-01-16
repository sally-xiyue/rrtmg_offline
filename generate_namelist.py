import argparse
import json
import pprint
from sys import exit
import uuid

def main():
    parser = argparse.ArgumentParser(prog='Namelist Generator')
    parser.add_argument('case_name')
    args = parser.parse_args()

    case_name = args.case_name

    if case_name == 'GCMVarying':
        namelist = GCMVarying()
    else:
        print('Not a valid case name')
        exit()

    write_file(namelist)


def GCMVarying():

    namelist = {}

    namelist['input'] = {}
    namelist['input']['path'] = '/Users/xiyue/Clouds/GCM/output/Output.GCMVarying_seaice_1_50x_lat70_lon0_sea_ice.73e98/stats/'
    namelist['input']['case'] = 'GCMVarying_seaice_1_50x_lat70_lon0_sea_ice'
    namelist['input']['file'] = 'Stats.'+namelist['input']['case']+'.nc'
    namelist['input']['t1'] = 0*4
    namelist['input']['t2'] = 25*4
    namelist['input']['time_average'] = False

    namelist['input']['albedo_path'] = '/Users/xiyue/Clouds/GCM/output/Output.GCMVarying_seaice_1_00x_lat70_lon0_sea_ice.f099b/stats/'
    namelist['input']['albedo_case'] = 'GCMVarying_seaice_1_00x_lat70_lon0_sea_ice'
    namelist['input']['albedo_file'] = 'Stats.'+namelist['input']['albedo_case']+'.nc'

    namelist['input']['fix_T'] = True
    namelist['input']['fix_qv'] = True
    namelist['input']['fix_cloud'] = False
    namelist['input']['fix_albedo'] = True

    namelist['output'] = {}
    namelist['output']['file'] = 'rrtm_output_0_liq.pkl'

    namelist['radiation'] = {}
    namelist['radiation']['co2_factor'] = 1.0*4
    namelist['radiation']['use_RRTM'] = True
    namelist['radiation']['frequency'] = 300.0
    namelist['radiation']['n_buffer'] = 2 # adjust according to dz
    namelist['radiation']['stretch_factor'] = 1.0 # adjust according to dz

    namelist['time_stepping'] = {}
    namelist['time_stepping']['t'] = 0.0
    namelist['time_stepping']['dt_initial'] = 300.0
    # namelist['time_stepping']['dt_max'] = 60.0
    namelist['time_stepping']['t_max'] = 3600.0 * 24.0

    namelist['stats_io'] = {}
    namelist['stats_io']['frequency'] = 300.0
    # namelist['stats_io']['output_root'] = './output/'
    namelist['stats_io']['output_root'] = '/Users/xiyue/Clouds/mlm/output/data/'

    namelist['meta'] = {}
    namelist['meta']['simname'] = 'GCMVarying'
    namelist['meta']['casename'] = 'GCMVarying'

    return namelist


def write_file(namelist):

    try:
        type(namelist['meta']['simname'])
    except:
        print('Casename not specified in namelist dictionary!')
        print('FatalError')
        exit()

    namelist['meta']['uuid'] = str(uuid.uuid4())

    fh = open(namelist['meta']['simname'] + '.in', 'w')
    pprint.pprint(namelist)
    json.dump(namelist, fh, sort_keys=True, indent=4)
    fh.close()

    return


if __name__ == '__main__':
    main()
