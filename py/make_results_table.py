import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import plot_perplex as plotpx
import saturation as sat
import main as rw
import parameters as p


def across_masses(core_eff, masses=None, Tp=1600, output_path=px.perplex_path_default + 'output/', end='',
                  scale=p.TO ** -1, q=[2.27, 50, 97.93]):
    s_Fe = '{:.2f}'.format(1 - core_eff)
    line_m = ' & ' + s_Fe + ' '
    line_um = ' & ' + s_Fe + ' '

    mass_str = []
    for Mp in masses:
        if isinstance(Mp, float):  # Tp = 1600
            mass_str.append(str(Mp).replace('.', ','))
        elif isinstance(Mp, int):
            mass_str.append(str(Mp))
        else:
            print('MP is', type(Mp))
    dirs = ['hypatia' + ms + 'M_' + str(Tp) + 'K_' + str(int(core_eff * 100)) + 'Fe' + end for ms in mass_str]

    for ms, d in zip(mass_str, dirs):
        dats = rw.read_dir(output_path + d + '/', verbose=True, prevent_blank_dir=False)

        # get water masses
        w_m = []
        w_um = []
        for dat in dats:
            try:
                w_m.append(dat.mass_h2o_total * scale)
                w_um.append(dat.mass_h2o_um * scale)
            except (AttributeError, KeyError):
                pass  # blank data, didn't work because star not measured maybe
            except TypeError:
                # attr is None?
                pass

        # get stats
        try:
            w_m_min, w_m_med, w_m_max = np.percentile(w_m, q, method='nearest')
            w_um_min, w_um_med, w_um_max = np.percentile(w_um, q, method='nearest')
        except IndexError:  # size 0
            w_m_min, w_m_med, w_m_max = 0, 0, 0
            w_um_min, w_um_med, w_um_max = 0, 0, 0

        line_m = line_m + '& {:.2f} ({:.2f}) '.format(w_m_med, w_m_max - w_m_min)
        line_um = line_um + '& {:.2f} ({:.2f}) '.format(w_um_med, w_um_max - w_um_min)

    line_m = line_m + r'\\'
    line_um = line_um + r'\\'
    return line_m, line_um


def across_fe(Mp, core_effs=[0.7, 0.8, 0.88, 0.95, 0.99], Tp=1600, output_path=px.perplex_path_default + 'output/',
              end='',
              scale=p.TO ** -1, q=[2.27, 50, 97.93]):
    # s = '{:.1f}$\,M_\oplus$'.format(Mp)
    s = '{:.1f}'.format(Mp)
    line_m = ' & ' + s + ' '
    line_um = ' & ' + s + ' '

    if isinstance(Mp, float):  # Tp = 1600
        mass_str = str(Mp).replace('.', ',')
    elif isinstance(Mp, int):
        mass_str = str(Mp)
    else:
        raise NotImplementedError('MP is', type(Mp))
    dirs = ['hypatia' + mass_str + 'M_' + str(Tp) + 'K_' + str(int(fe * 100)) + 'Fe' + end for fe in core_effs]

    for fe, d in zip(core_effs, dirs):
        dats = rw.read_dir(output_path + d + '/', prevent_blank_dir=False)

        # get water masses
        w_m = []
        w_um = []
        for dat in dats:
            try:
                w_m.append(dat.mass_h2o_total * scale)
                w_um.append(dat.mass_h2o_um * scale)
            except (AttributeError, KeyError):
                pass  # blank data, didn't work because star not measured maybe
            except TypeError:
                # attr is None?
                pass

        # get stats
        try:
            w_m_min, w_m_med, w_m_max = np.percentile(w_m, q, method='nearest')
            w_um_min, w_um_med, w_um_max = np.percentile(w_um, q, method='nearest')
        except IndexError:  # size 0
            w_m_min, w_m_med, w_m_max = 0, 0, 0
            w_um_min, w_um_med, w_um_max = 0, 0, 0

        line_m = line_m + '& {:.2f} ({:.2f}) '.format(w_m_med, w_m_max - w_m_min)
        line_um = line_um + '& {:.2f} ({:.2f}) '.format(w_um_med, w_um_max - w_um_min)

    line_m = line_m + r'\\'
    line_um = line_um + r'\\'
    return line_m, line_um


def write_table(masses=None,
                         core_effs=None,
                         temps=[1600, 1900]):
    if core_effs is None:
        core_effs = [0.7, 0.8, 0.88, 0.95, 0.99]
    if masses is None:
        masses = [0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5]
    s_m = [r'\multicolumn{' + str(len(core_effs) + 2) + r'}{c}{\textbf{Whole mantle water capacity (OM)}} \\']
    s_um = [r'\multicolumn{' + str(len(core_effs) + 2) + r'}{c}{\textbf{Upper mantle water capacity (OM)}} \\']

    for Tp in temps:
        tab_t_m = []
        tab_t_um = []
        for im, Mp in enumerate(masses):
            lm, lum = across_fe(Mp, core_effs, Tp=Tp, output_path=px.perplex_path_default + 'output/apollo/',
                                end='_hires', scale=p.TO ** -1, q=[2.27, 50, 97.93])
            if im == 0:
                lm = r'\multirow{' + str(len(masses)) + r'}{*}{' + str(Tp) + '}' + lm
                lum = r'\multirow{' + str(len(masses)) + r'}{*}{' + str(Tp) + '}' + lum
                # lm = r'\multirow{' + str(len(masses)) + r'}{*}{$' + str(Tp) + r'\,{\rm K}$}' + lm
                # lum = r'\multirow{' + str(len(masses)) + r'}{*}{$' + str(Tp) + r'\,{\rm K}$}' + lum
            tab_t_m.append(lm)
            tab_t_um.append(lum)

        s_m.extend(tab_t_m)
        s_m.append(r'\hline')

        s_um.extend(tab_t_um)
        s_um.append(r'\hline')

    print('\n\n\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%\n% start copy\n\n\n\n')
    for s in s_m:
        print(s)  # remove last hline

    for s in s_um:
        print(s)  # remove last hline

    print('\n\n\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%\n% end copy\n\n')


write_table()
