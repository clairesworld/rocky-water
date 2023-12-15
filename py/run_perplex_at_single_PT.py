import ask_hypatia as hyp
import main as rw
import perplexdata as px

perplex_path = '/home/g/guimond/Work/perple_x/'
oxide_list = ['MgO', 'SiO2', 'CaO', 'Al2O3', 'FeO']  # 'Na2O

core_Si_wtpt=0  # percentage
core_eff=0.5  # frac

""" run perplex at 4 GPa only """
all_comps = hyp.bulkcomp_from_hypatiatsv(
    core_eff=core_eff, core_Si_wtpt=core_Si_wtpt, oxide_list=oxide_list, n_from_top=-1,
    path_to_tsv='/home/g/guimond/Work/hypatia-compositions/hypatia-04122023.tsv')

# P, T = 4, 1659
P, T = 30, 1907
for attrs in all_comps:
    pl = rw.phases_at_single_pt(pressure=P, temperature=T,
                                output_parent_path=perplex_path + 'output/coreSi' + str(core_Si_wtpt) + '_coreeff' + str(core_eff).replace('.', ',') + '/',
                                test_oxides=attrs['wt_oxides'], test_CMF=attrs['CMF'], pickle=True, verbose=True,
                                perplex_path=perplex_path,
                                suppress_output=True, clean=True, **attrs)


