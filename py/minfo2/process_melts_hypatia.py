import numpy as np
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(os.path.dirname(PARENT_DIR))

import py.minfo2.meltsfugacitydata as mfug

"""
source /raid1/cmg76/venv/bin/activate
cd ~/Works/rocky-water/
git pull
python3 py/minfo2/process_melts_hypatia.py
"""

# set paths
alphamelts_path_apollo = '/raid1/cmg76/alphamelts/'
opp_apollo = '/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/'
alphamelts_path_starlite = '/home/claire/Works/alphamelts/'
opp_starlite = '/home/claire/Works/min-fo2/alphamelts_output/hypatia_local2/'
opp_galson = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'

""" vvv UNCOMMENT TO RUN vvv """

# set these
T_of_interest = 1573.15  # 1673.15
X_ferric = [float(x) for x in input('Enter X_ferric, separated by spaces (e.g. 0.03): ').split()] #[0.07]  #, 0.03, 0.05, 0.07, 0.09]  # , 0.01, 0.05, 0.07, 0.09]  #[0.01, 0.03, 0.05, 0.07, 0.09]
core_eff = [float(x) for x in input('Enter core_eff, separated by spaces (e.g. 0.88): ').split()] #[0.88]
# core_eff = [0.88]
# X_ferric = [0.03]
location = 'apollo'  # 'starlite'

# run
if location == 'apollo':
    source = opp_apollo
    alphamelts_path = alphamelts_path_apollo
    perplex_path = '/raid1/cmg76/perple_x/'
elif location == 'starlite':
    source = opp_starlite
    alphamelts_path = alphamelts_path_starlite
    perplex_path = '/home/claire/Works/perple_x/'

for ce in core_eff:
    for Xf in X_ferric:
        output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
        output_parent_path = source + output_sub
        # output_parent_path = '/raid1/cmg76/alphamelts/output/rocky-fo2/hypatia_88coreeff_3ferric_ext_Cr/'

        # calculate mantle fo2 only
        mfug.fo2_from_local(output_parent_path, core_efficiency=ce, X_ferric=Xf, alphamelts_path=alphamelts_path,
                            compare_buffer='qfm', perplex_path=perplex_path, T_of_interest=T_of_interest, save=True,
                            verbose=False,
                            # names=['1M_88Ceff_2MASS00182469-1516022_999K_3,0fer']
                            # restart='2MASS19064452+4705535',
                            # names=#['1M_88Ceff_HIP80680_999K_3,0fer', '1M_88Ceff_HIP53719_999K_3,0fer', '1M_88Ceff_HIP7513_999K_3,0fer', '1M_88Ceff_HIP25191_999K_3,0fer', '1M_88Ceff_2MASS19172334+4412307_999K_3,0fer', '1M_88Ceff_2MASS19282502+3946044_999K_3,0fer', '1M_88Ceff_2MASS19500236+4657405_999K_3,0fer', '1M_88Ceff_2MASS19115949+5056395_999K_3,0fer', '1M_88Ceff_2MASS19455215+4235555_999K_3,0fer', '1M_88Ceff_HIP20723_999K_3,0fer', '1M_88Ceff_2MASS19202179+3949011_999K_3,0fer', '1M_88Ceff_2MASS14123753+0403359_999K_3,0fer', '1M_88Ceff_2MASS19245404+4455385_999K_3,0fer', '1M_88Ceff_2MASS10401438-6227201_999K_3,0fer', '1M_88Ceff_2MASS19344235+4436560_999K_3,0fer', '1M_88Ceff_2MASS19511083+4025037_999K_3,0fer', '1M_88Ceff_HIP72339_999K_3,0fer', '1M_88Ceff_HIP49699_999K_3,0fer', '1M_88Ceff_2MASS19012332+4145429_999K_3,0fer', '1M_88Ceff_2MASS19235374+4810413_999K_3,0fer', '1M_88Ceff_2MASS19154176+4603274_999K_3,0fer', '1M_88Ceff_BD+2002184_999K_3,0fer', '1M_88Ceff_2MASS18484802+4425202_999K_3,0fer', '1M_88Ceff_2MASS16163403-2024019_999K_3,0fer', '1M_88Ceff_2MASS19122295+4021574_999K_3,0fer', '1M_88Ceff_HIP74890_999K_3,0fer', '1M_88Ceff_2MASS18553242+4738139_999K_3,0fer', '1M_88Ceff_2MASS19091838+4840243_999K_3,0fer', '1M_88Ceff_2MASS19345587+4154030_999K_3,0fer', '1M_88Ceff_HIP14810_999K_3,0fer', '1M_88Ceff_HIP83949_999K_3,0fer', '1M_88Ceff_HIP68162_999K_3,0fer', '1M_88Ceff_HIP5806_999K_3,0fer', '1M_88Ceff_2MASS19461201+4906083_999K_3,0fer', '1M_88Ceff_2MASS19323844+4852522_999K_3,0fer', '1M_88Ceff_2MASS19580041+4040148_999K_3,0fer', '1M_88Ceff_HIP6993_999K_3,0fer', '1M_88Ceff_HIP108859_999K_3,0fer', '1M_88Ceff_2MASS19413057+3902529_999K_3,0fer', '1M_88Ceff_2MASS11013589-2351382_999K_3,0fer', '1M_88Ceff_2MASS19141179+3833548_999K_3,0fer', '1M_88Ceff_HIP6511_999K_3,0fer', '1M_88Ceff_HIP95740_999K_3,0fer', '1M_88Ceff_2MASS19213437+4321500_999K_3,0fer', '1M_88Ceff_2MASS19025218+4445303_999K_3,0fer', '1M_88Ceff_2MASS19480452+5024323_999K_3,0fer', '1M_88Ceff_2MASS18172957-0322517_999K_3,0fer', '1M_88Ceff_HIP103527_999K_3,0fer', '1M_88Ceff_HIP85294_999K_3,0fer', '1M_88Ceff_HIP111136_999K_3,0fer', '1M_88Ceff_2MASS19063018+3732142_999K_3,0fer', '1M_88Ceff_2MASS19203966+4955259_999K_3,0fer', '1M_88Ceff_HIP99825_999K_3,0fer', '1M_88Ceff_HIP80337_999K_3,0fer', '1M_88Ceff_HIP5529_999K_3,0fer', '1M_88Ceff_HIP80250_999K_3,0fer', '1M_88Ceff_2MASS19125618+4031152_999K_3,0fer', '1M_88Ceff_2MASS19234989+4724226_999K_3,0fer', '1M_88Ceff_HIP80902_999K_3,0fer', '1M_88Ceff_HIP106353_999K_3,0fer', '1M_88Ceff_2MASS19224155+3841276_999K_3,0fer', '1M_88Ceff_HIP801_999K_3,0fer', '1M_88Ceff_HIP67246_999K_3,0fer', '1M_88Ceff_HIP102125_999K_3,0fer', '1M_88Ceff_2MASS19393877+4629292_999K_3,0fer', '1M_88Ceff_2MASS19515301+4743540_999K_3,0fer', '1M_88Ceff_2MASS18504798+4525327_999K_3,0fer', '1M_88Ceff_2MASS19345473+4607449_999K_3,0fer', '1M_88Ceff_HIP84787_999K_3,0fer', '1M_88Ceff_2MASS18584254+4447516_999K_3,0fer', '1M_88Ceff_HIP46076_999K_3,0fer', '1M_88Ceff_HIP22336_999K_3,0fer', '1M_88Ceff_2MASS18525313+4025185_999K_3,0fer', '1M_88Ceff_2MASS19392772+4617090_999K_3,0fer', '1M_88Ceff_2MASS19440088+4416392_999K_3,0fer', '1M_88Ceff_2MASS19474849+4316190_999K_3,0fer', '1M_88Ceff_2MASS19294147+3815587_999K_3,0fer', '1M_88Ceff_2MASS19284107+4054587_999K_3,0fer', '1M_88Ceff_HD17092_999K_3,0fer', '1M_88Ceff_2MASS19471282+4644281_999K_3,0fer', '1M_88Ceff_2MASS19172608+5035545_999K_3,0fer', '1M_88Ceff_2MASS19005780+4640057_999K_3,0fer', '1M_88Ceff_2MASS19530049+4029458_999K_3,0fer', '1M_88Ceff_2MASS19083376+4706547_999K_3,0fer', '1M_88Ceff_2MASS19304566+3750059_999K_3,0fer', '1M_88Ceff_2MASS19052593+4224234_999K_3,0fer', '1M_88Ceff_2MASS19192612+4841378_999K_3,0fer', '1M_88Ceff_2MASS19251253+4741519_999K_3,0fer', '1M_88Ceff_2MASS18515993+4033268_999K_3,0fer', '1M_88Ceff_2MASS19174668+4119039_999K_3,0fer', '1M_88Ceff_2MASS14474655+0103538_999K_3,0fer', '1M_88Ceff_HIP77740_999K_3,0fer', '1M_88Ceff_2MASS19213918+3820375_999K_3,0fer', '1M_88Ceff_2MASS19190557+4048026_999K_3,0fer', '1M_88Ceff_2MASS07324421+3350061_999K_3,0fer', '1M_88Ceff_2MASS16000805-2311213_999K_3,0fer', '1M_88Ceff_2MASS19404359+4956448_999K_3,0fer', '1M_88Ceff_HIP80484_999K_3,0fer', '1M_88Ceff_2MASS18581576+4146425_999K_3,0fer', '1M_88Ceff_2MASS19481641+4650034_999K_3,0fer', '1M_88Ceff_2MASS19004386+4349519_999K_3,0fer', '1M_88Ceff_2MASS19190999+3917070_999K_3,0fer', '1M_88Ceff_2MASS18573463+4614566_999K_3,0fer', '1M_88Ceff_2MASS18583244+4043113_999K_3,0fer', '1M_88Ceff_HIP118319_999K_3,0fer', '1M_88Ceff_2MASS19404641+3932228_999K_3,0fer', '1M_88Ceff_2MASS19393833+4256069_999K_3,0fer', '1M_88Ceff_2MASS19215082+4033448_999K_3,0fer', '1M_88Ceff_2MASS08461929-0801370_999K_3,0fer', '1M_88Ceff_HIP8102_999K_3,0fer', '1M_88Ceff_2MASS19022767+5008087_999K_3,0fer', '1M_88Ceff_HIP57087_999K_3,0fer', '1M_88Ceff_2MASS19371607+4337456_999K_3,0fer', '1M_88Ceff_2MASS19121622+4221193_999K_3,0fer', '1M_88Ceff_HIP83983_999K_3,0fer', '1M_88Ceff_HIP20199_999K_3,0fer', '1M_88Ceff_2MASS04245952+3927382_999K_3,0fer', '1M_88Ceff_2MASS19514517+4652530_999K_3,0fer', '1M_88Ceff_2MASS19441136+4244348_999K_3,0fer', '1M_88Ceff_2MASS19152144+4659122_999K_3,0fer', '1M_88Ceff_2MASS19102533+4931237_999K_3,0fer', '1M_88Ceff_B-1003166_999K_3,0fer', '1M_88Ceff_2MASS19173311+5035493_999K_3,0fer', '1M_88Ceff_2MASS19422610+3844086_999K_3,0fer', '1M_88Ceff_2MASS19441436+4239217_999K_3,0fer', '1M_88Ceff_HIP50786_999K_3,0fer',
                            #        # ['1M_88Ceff_B-0103943_999K_3,0fer']
                            # ['1M_88Ceff_2MASS19330772+4817092_999K_3,0fer', '1M_88Ceff_HIP78169_999K_3,0fer', '1M_88Ceff_2MASS19050572+4837073_999K_3,0fer', '1M_88Ceff_HD240210_999K_3,0fer', '1M_88Ceff_2MASS19444402+4721315_999K_3,0fer', '1M_88Ceff_2MASS19245834+4400313_999K_3,0fer', '1M_88Ceff_2MASS19175659+5056389_999K_3,0fer', '1M_88Ceff_HIP65721_999K_3,0fer', '1M_88Ceff_BD+4602726_999K_3,0fer', '1M_88Ceff_HIP47202_999K_3,0fer', '1M_88Ceff_2MASS19055315+4938564_999K_3,0fer', '1M_88Ceff_2MASS19272912+4405145_999K_3,0fer', '1M_88Ceff_2MASS19490693+4343264_999K_3,0fer', '1M_88Ceff_2MASS22493256-1040320_999K_3,0fer', '1M_88Ceff_2MASS19445520+5016308_999K_3,0fer', '1M_88Ceff_2MASS19211746+4402089_999K_3,0fer', '1M_88Ceff_2MASS19295686+4611463_999K_3,0fer', '1M_88Ceff_HIP21109_999K_3,0fer', '1M_88Ceff_2MASS19231995+3811036_999K_3,0fer', '1M_88Ceff_2MASS22151722-1402593_999K_3,0fer', '1M_88Ceff_2MASS19535546+4136552_999K_3,0fer', '1M_88Ceff_2MASS19013718+4959410_999K_3,0fer', '1M_88Ceff_HIP88348_999K_3,0fer']

                            )
        restart = None  # reset

        # print('\n\n\nfinding common T min\n--------------------\n')
        # mfug.common_Tmin(output_parent_path, include_p=[1e4, 4e4])

""" ^^^ UNCOMMENT TO RUN ^^^ """

""" custom """
# output_parent_path = '/raid1/cmg76/alphamelts/output/rocky-fo2/mgsi_from_earth/'
# mfug.common_Tmin(output_parent_path)
#
# # calculate mantle fo2 only
# ce, X_ferric = 0.88, 0.03
# mfug.fo2_from_local(output_parent_path, core_efficiency=ce, X_ferric=X_ferric, alphamelts_path=alphamelts_path,
#                     compare_buffer='qfm', perplex_path=perplex_path,
#                     T_of_interest=T_of_interest,  # reload_TP=True,
#                     verbose=False)

#
""" for a min T look at all the cases that got there """
# source = opp_galson
# ce = 0.88
# Xf = 0.03
# T = 1373.15
# P = 1
# output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
# output_parent_path = source + output_sub
#
# mfug.create_isothermal_csv(output_parent_path, T, P, Xf, ce)
#
