import subprocess
import generate_namelist


def main():
    # schemes = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    dTis = [3]
    dSSTs = [0]
    rhs = [80, 90]

    for dTi in dTis:
        for dSST in dSSTs:
            for rh in rhs:
                nml = generate_namelist.Isdac()
                nml['meta']['simname'] = nml['meta'][
                    'simname'] + '_dTi_' + str(dTi) + '_dSST_' + str(dSST) \
                                         + '_rh_' + str(rh)
                nml['initial']['dTi'] = float(dTi)
                nml['initial']['dSST'] = float(dSST)
                nml['initial']['rh'] = rh/100.0

                generate_namelist.write_file(nml)
                run_str = 'bsub python main.py ' + \
                    nml['meta']['simname'] + '.in'
                print(run_str)
                subprocess.call([run_str], shell=True)

    return

if __name__ == "__main__":
    main()