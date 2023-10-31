import numpy as np
import pandas as pd
import parameters as p
import main as rw
import perplexdata as px

phase_columns = ['O', 'Opx', 'Cpx', 'Sp', 'Gt', 'Wad', 'Ring', 'Aki', 'Pv', 'Wus', 'C2/c', 'Ppv', 'CF', 'ca-pv', 'q', 'coes', 'st', ]
ox_columns = px.oxide_list_default
location = 'apollo'
pressures = [1]  # [1, 4, 30]
M_p = 1
core_effs = [88]  # [70, 80, 88, 95, 99]
Tp = 1600


def get_df_to_print_mineralogy(output_path, pressures, M_p=1, core_eff=88, Tp=1600, phase_columns=phase_columns,
                   ox_columns=ox_columns, subsample=None, test_mode=False, include_phases=True):

    dats = rw.read_dir(output_path, subsample=subsample, verbose=True, prevent_blank_dir=True)

    # initialise list of dfs for each pressure
    df_list = []
    for pressure in pressures:
        df_list.append(pd.DataFrame(columns=['star', 'P(GPa)', 'T(K)', 'Tp(K)']
                                            + [ph + '(wt%)' for ph in phase_columns] + ['M_p(M_E)', 'Fe_core/Fe_tot']
                                            + [ox + '(wt%)' for ox in ox_columns],
                                    index=range(len(dats))))

        # add constants
        df_list[-1]['M_p(M_E)'] = M_p
        df_list[-1]['Fe_core/Fe_tot'] = core_eff
        df_list[-1]['Tp(K)'] = Tp
        df_list[-1].fillna(0, inplace=True)

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

            if include_phases:
                for ph in phase_columns:
                    try:
                        df_list[z].loc[irow, ph + '(wt%)'] = ser[ph]
                    except KeyError as e:
                        if test_mode:
                            print(e, 'not found in', dat.name, 'mineralogy output at', pressure, 'GPa')
                        flag = True
                        pass
                if flag and test_mode:
                    print(ser)

            for ox in ox_columns:
                try:
                    df_list[z].loc[irow, ox + '(wt%)'] = dat.wt_oxides[ox]
                except KeyError as e:
                    if test_mode:
                        print(e, 'not found in', dat.name, 'wt oxides')
                    pass
    return df_list


def write_from_dir(output_parent_path, write_path, pressures, M_p=1, core_effs=[88], Tp=1600, xfer=3,
                   phase_columns=phase_columns, ox_columns=ox_columns,
                   output_path_override=None, sep=',', test_mode=False,
                   fname_base='hypatia_mineralogies',
                   include_earthsun=True, output_path_earthsun='/home/claire/Works/perple_x/output/references/',
                   **kwargs):

    for core_eff in core_effs:
        if output_path_override is None:
            # output_path = output_parent_path + 'hypatia{:d}M_{:d}K_{:d}Fe_hires/'.format(M_p, Tp, core_eff)
            output_path = output_parent_path + 'hypatia_{:d}coreeff_{:d}ferric_ext/'.format(core_eff, xfer)
        else:
            output_path = output_parent_path + output_path_override

        # read each planet in dir
        if test_mode:
            subsample = 3
        else:
            subsample = None

        df_list = get_df_to_print_mineralogy(output_path, pressures, M_p=M_p, core_eff=core_eff, Tp=Tp,
                                             phase_columns=phase_columns, ox_columns=ox_columns, subsample=subsample,
                                             test_mode=test_mode, **kwargs)

        if include_earthsun:
            df_earthsun_list = get_df_to_print_mineralogy(output_path_earthsun, pressures, M_p=M_p, core_eff=core_eff,
                                                          Tp=Tp, phase_columns=phase_columns, ox_columns=ox_columns,
                                                          subsample=subsample, test_mode=test_mode, **kwargs)
            df_list = [pd.concat([df_earthsun_list[ii], df_list[ii]], ignore_index=True) for ii in range(len(df_list))]

        # write csv
        for z, pressure in enumerate(pressures):
            fname = fname_base + '_{:d}GPa_{:d}K_{:d}core.csv'.format(pressure, Tp, core_eff)
            df_list[z].to_csv(write_path + fname, sep=sep)
            print('wrote to', write_path + fname)


if location == 'apollo':
    output_parent_path = '/raid1/cmg76/perple_x/output/rocky-water/'
    write_from_dir(output_parent_path, '/home/cmg76/Works/rocky-water/csv/',
                   fname_base='mantle_compositions_extended_',
                   include_earthsun=False, pressures=pressures, M_p=M_p,
                   core_effs=core_effs, Tp=Tp)

elif location == 'starlite':
    output_parent_path = '/home/claire/Works/perple_x/output/'
    write_from_dir(output_parent_path, '/home/claire/Works/rocky-water/csv/',
                   include_earthsun=True, output_path_earthsun='/home/claire/Works/perple_x/output/references/',
                   pressures=pressures, M_p=M_p, core_effs=core_effs, Tp=Tp)
