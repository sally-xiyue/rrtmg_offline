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
        exit()0

    write_file(namelist)


def GCMVarying():

    namelist = {}

    namelist['input'] = {}
    namelist['input']['path'] = '/Users/sallyz/Dropbox/WORK/GCM/output/Output.GCMVarying_seaice_1_00x_lat70_lon90_zero.d3123/stats/'
    namelist['input']['case'] = 'GCMVarying_seaice_1_00x_lat70_lon90_zero'
    namelist['input']['file'] = 'Stats.'+namelist['input']['case']+'.nc'
    namelist['input']['t1'] = 0*4
    namelist['input']['t2'] = 25*4
    namelist['input']['time_average'] = False

    namelist['input']['albedo_path'] = '/Users/sallyz/Dropbox/WORK/GCM/output/Output.GCMVarying_seaice_1_00x_lat70_lon90_zero.d3123/stats/'
    namelist['input']['albedo_case'] = 'GCMVarying_seaice_1_00x_lat70_lon90_zero'
    namelist['input']['albedo_file'] = 'Stats.'+namelist['input']['albedo_case']+'.nc'

    namelist['input']['fix_T'] = True
    namelist['input']['fix_qv'] = True
    namelist['input']['fix_cloud'] = False
    namelist['input']['fix_albedo'] = True
    namelist['input']['no_ice'] = False

    namelist['radiation'] = {}
    namelist['radiation']['co2_factor'] = 1.0#*4.0
    namelist['radiation']['use_RRTM'] = True
    namelist['radiation']['frequency'] = 300.0
    namelist['radiation']['n_buffer'] = 2 # adjust according to dz
    namelist['radiation']['stretch_factor'] = 1.0 # adjust according to dz

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
