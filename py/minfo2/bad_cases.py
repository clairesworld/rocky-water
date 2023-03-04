

bad_cases = ['1M_88Ceff_HIP13291_999K_3,0fer', '1M_88Ceff_HIP84607_999K_3,0fer', '1M_88Ceff_GaiaDR25242077854771196544_999K_3,0fer', '1M_88Ceff_BD+4602629_999K_3,0fer']
bad_cases = ['1M_88Ceff_HIP13291_999K_9,0fer', '1M_88Ceff_HIP84607_999K_9,0fer', '1M_88Ceff_GaiaDR25242077854771196544_999K_9,0fer', '1M_88Ceff_BD+4602629_999K_9,0fer']


""" 
MELTS

88Ceff 3fer
alphamelts infinite loop: HIP 102409 - also at 5fer
has kyanite: 1M_88Ceff_2MASS19141179+3833548_999K_3,0fer (also quartz so fucked) --> removing (also 7, 1, 5 fer) (also 70ceff)

others:
1M_88Ceff_2MASS19211746+4402089_999K_3,0fer did not complete, no Phase_mass_tbl.txt found in /raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/hypatia_88coreeff_3ferric_ext/1M_88Ceff_2MASS19211746+4402089_999K_3,0fer/10000bar/
1M_88Ceff_2MASS19490218+4650354_999K_3,0fer did not complete, no Phase_mass_tbl.txt found in /raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/hypatia_88coreeff_3ferric_ext/1M_88Ceff_2MASS19490218+4650354_999K_3,0fer/10000bar/
1M_88Ceff_2MASS14474655+0103538_999K_3,0fer did not complete, no Phase_mass_tbl.txt found in /raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/hypatia_88coreeff_3ferric_ext/1M_88Ceff_2MASS14474655+0103538_999K_3,0fer/10000bar/
1M_88Ceff_2MASS19213437+4321500_999K_3,0fer did not complete, no Phase_mass_tbl.txt found in /raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/hypatia_88coreeff_3ferric_ext/1M_88Ceff_2MASS19213437+4321500_999K_3,0fer/10000bar/
1M_88Ceff_2MASS19125618+4031152_999K_3,0fer did not complete, no Phase_mass_tbl.txt found in /raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/hypatia_88coreeff_3ferric_ext/1M_88Ceff_2MASS19125618+4031152_999K_3,0fer/10000bar/
1M_88Ceff_HD89484_999K_3,0fer did not complete, no Phase_mass_tbl.txt found in /raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/hypatia_88coreeff_3ferric_ext/1M_88Ceff_HD89484_999K_3,0fer/10000bar/

70Ceff 3fer
alphamelts infinite loop: HIP 15578, HIP 17054, HIP 69888

HIP 13291 star no Ti

2MASS19064452+4705535_999K_3,0fer/14285bar/ has rhm-oxide but not 10 kbar, 40 kbar
1M_80Ceff_2MASS04215269+5749018_999K_3,0fer/10000 bar has rhm-oxide, removing

Spinel missing (Al is in cpx) - 480 cases 88/3 ?!?!?

Cr
1M_88Ceff_2MASS19141179+3833548_999K_3,0fer kyanite after 10 kbar <- removed
a bunch of these (dozens) had alphamelts not complete
--> might be that:
 Your choice: Superliquidus (1) or subsolidus (0) initial guess ? Initial calculation failed (40000.000000 bars, 2273.150000 K)!


-----------------------------------------------------------------------------------------------
PERPLE_X

cases that crashed/didn't finish on initial runs:

['1M_88Ceff_HIP13291_999K_3,0fer', '1M_88Ceff_HIP84607_999K_3,0fer', '1M_88Ceff_GaiaDR25242077854771196544_999K_3,0fer', '1M_88Ceff_BD+4602629_999K_3,0fer']
1M_88Ceff_HIP84607_999K_3,0fer$ 
['1M_85Ceff_HIP13291_999K_3,0fer', '1M_85Ceff_HIP84607_999K_3,0fer', '1M_85Ceff_GaiaDR25242077854771196544_999K_3,0fer', '1M_85Ceff_BD+4602629_999K_3,0fer']
['1M_80Ceff_HIP13291_999K_3,0fer', '1M_80Ceff_HIP84607_999K_3,0fer', '1M_80Ceff_GaiaDR25242077854771196544_999K_3,0fer', '1M_80Ceff_BD+4602629_999K_3,0fer']
['1M_70Ceff_HIP13291_999K_3,0fer', '1M_70Ceff_HIP84607_999K_3,0fer', '1M_70Ceff_GaiaDR25242077854771196544_999K_3,0fer', '1M_70Ceff_BD+4602629_999K_3,0fer']
['1M_99Ceff_HIP13291_999K_3,0fer', '1M_99Ceff_HIP84607_999K_3,0fer', '1M_99Ceff_GaiaDR25242077854771196544_999K_3,0fer', '1M_99Ceff_BD+4602629_999K_3,0fer']


88/3%: have run but missing vertex files:
HIP52733
HIP13291


Cr
2MASS19190557+4048026 crashed unknown reason


88/3% Fe3+ compositon did not complete to 4 GPa for some reason, maybe badly formatted werami input file
cases_missingfe2o3 = ['1M_88Ceff_HIP69888_999K_3,0fer', '1M_88Ceff_HIP69888_999K_3,0fer', '1M_88Ceff_HIP69888_999K_3,0fer',
                             '1M_88Ceff_2MASS19492623+4947511_999K_3,0fer', '1M_88Ceff_2MASS19492623+4947511_999K_3,0fer',
                             '1M_88Ceff_2MASS19492623+4947511_999K_3,0fer', '1M_88Ceff_HIP31039_999K_3,0fer',
                             '1M_88Ceff_HIP31039_999K_3,0fer', '1M_88Ceff_HIP31039_999K_3,0fer', '1M_88Ceff_HIP77838_999K_3,0fer',
                             '1M_88Ceff_HIP77838_999K_3,0fer', '1M_88Ceff_HIP77838_999K_3,0fer', '1M_88Ceff_2MASS19330772+4817092_999K_3,0fer',
                             '1M_88Ceff_2MASS19330772+4817092_999K_3,0fer', '1M_88Ceff_2MASS19330772+4817092_999K_3,0fer',
                             '1M_88Ceff_2MASS19231995+3811036_999K_3,0fer', '1M_88Ceff_2MASS19231995+3811036_999K_3,0fer',
                             '1M_88Ceff_2MASS19231995+3811036_999K_3,0fer']




"""