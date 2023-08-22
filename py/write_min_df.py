import numpy as np
import pandas as pd
import parameters as p
import main as rw
import perplexdata as px

# phase_columns = ['gt', 'cpx', 'opx', 'hpcpx', 'ol', 'wad', 'ring', 'pv', 'qtz', 'coes', 'st', 'wus', 'dvm']
phase_columns = ['O', 'Sp', 'Cpx', 'Wad', 'Ring', 'Pv', 'Wus', 'C2/c', 'Opx', 'Aki', 'Ppv', 'Gt', 'CF', 'st', 'ca-pv']
ox_columns = px.oxide_list_default
location = 'apollo'


def write_from_dir(output_parent_path, write_path, pressures, M_p=1, core_eff=88, Tp=1600, phase_columns=phase_columns,
                   ox_columns=ox_columns, output_path_override=None, sep=',', **kwargs):

    if output_path_override is None:
        output_path = output_parent_path + 'hypatia{:d}M_{:d}K_{:d}Fe_hires/'.format(M_p, Tp, core_eff)
    else:
        output_path = output_parent_path + output_path_override

    # read each planet in dir
    dats = rw.read_dir(output_path, subsample=2, verbose=True, prevent_blank_dir=True)

    # initialise list of dfs for each pressure
    df_list = []
    for pressure in pressures:
        df_list.append(pd.DataFrame(columns=['star', 'M_p(M_E)', 'Fe_core/Fe_tot', 'P(GPa)', 'T(K)', 'Tp(K)']
                                            + [ph + '(wt%)' for ph in phase_columns]
                                            + [ox + '(wt%)' for ox in ox_columns],
                            index=range(len(dats))))

        # add constants
        df_list[-1]['M_p(M_E)'] = M_p
        df_list[-1]['Fe_core/Fe_tot'] = core_eff
        df_list[-1]['Tp(K)'] = Tp

    # load data for this planet at all pressures
    for irow, dat in enumerate(dats):

        if dat.df_comp is not None:
            pass
        else:
            print('dat.df_comp is None', dat.name)

        for z, pressure in enumerate(pressures):
            flag = False

            # find pressure
            ii = dat.df_comp['P(bar)'].sub(pressure * 1e4).abs().idxmin()
            ser = dat.df_comp.iloc[ii]

            df_list[z].loc[irow, 'star'] = dat.star
            df_list[z].loc[irow, 'P(GPa)'] = ser['P(bar)'] / 1e4
            df_list[z].loc[irow, 'T(K)'] = ser['T(K)']

            for ph in phase_columns:
                try:
                    df_list[z].loc[irow, ph + '(wt%)'] = ser[ph]
                except KeyError as e:
                    print(e, 'not found in', dat.name, 'mineralogy output at', pressure, 'GPa')
                    flag = True
                    pass
            if flag:
                print(ser)

            for ox in ox_columns:
                try:
                    df_list[z].loc[irow, ox + '(wt%)'] = dat.wt_oxides[ox]
                except KeyError as e:
                    print(e, 'not found in', dat.name, 'wt oxides')
                    pass

    # write csv
    for z, pressure in enumerate(pressures):
        fname = 'hypatia_mineralogies_{:d}GPa_{:d}K'.format(pressure, Tp)
        df_list[z].to_csv(write_path + fname, sep=sep)
        print('wrote to', write_path + fname)


if location == 'apollo':
    output_parent_path = '/raid1/cmg76/perple_x/output/rocky-water/'
    write_from_dir(output_parent_path, '/home/cmg76/Works/rocky-water/csv/', pressures=[1, 4, 30], M_p=1,
                   core_eff=88, Tp=1600)
