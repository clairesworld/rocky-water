import ask_hypatia as hyp
import main as rw
import perplexdata as px

perplex_path = '/home/g/guimond/Work/perple_x/'
oxide_list = ['MgO', 'SiO2', 'CaO', 'Al2O3', 'FeO']  # 'Na2O

""" run perplex at 4 GPa only """
all_comps = hyp.bulkcomp_from_hypatiatsv(
    core_eff=0.999, core_Si_wtpt=0, oxide_list=oxide_list, n_from_top=-1,
    path_to_tsv='/home/g/guimond/Work/hypatia-compositions/hypatia-04122023.tsv')

for attrs in all_comps:
    pl = rw.phases_at_single_pt(pressure=4, temperature=1659,
                                output_parent_path=perplex_path + 'output/coreSi0_coreeff0,999/',
                                test_oxides=attrs['wt_oxides'], test_CMF=attrs['CMF'], pickle=True, verbose=True,
                                perplex_path=perplex_path,
                                suppress_output=True, clean=True, **attrs)


