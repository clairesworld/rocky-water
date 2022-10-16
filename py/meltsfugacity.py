import numpy as np
import perplexdata as px
import bulk_composition as bulk

output_parent_default = '/home/claire/Works/min-fo2/alpha_output/'
melts_path_default = '/home/claire/Works/alphamelts/'
solution_phases_default = ['olivine', 'spinel', 'garnet']  # opx and cpx automatic?

class MeltsFugacity():

    def __init__(self, name='test', star=None, solution_phases=solution_phases_default,
                 output_parent_path=output_parent_default, wt_oxides=px.wt_oxides_DMM, X_ferric=0.03,
                 melts_path=melts_path_default,
                 **kwargs):
        # super().__init__(output_parent_path=output_parent_path, **kwargs)

        self.name = name
        self.star = star

        if 'O2' not in wt_oxides:
            # get O2 mantle fraction
            self.wt_oxides = bulk.o2_from_wt_comp(wt_oxides, X_ferric=X_ferric)
        else:
            self.wt_oxides = wt_oxides
        self.oxide_list = [k for k in self.wt_oxides.keys()]
        self.get_mgsi()
        self.get_mg_number()
        self.solution_phases = solution_phases
        self.output_path = output_parent_path + self.name + '/'
        self.melts_path = melts_path
        self.melts_fname = self.output_path + self.name + '.melts'


    def command_melts(self, options_file=None, **kwargs):
        """ string for alphamelts command file to get fo2 """

        s = '1\n'  #  1. Read MELTS file to set composition of system
        s = s + self.melts_fname + '\n' # MELTS filename:
        s = s + '9\n'  #  9. Turn liquid on / off
        s = s + '0\n'  # Turn liquid on (1) or suppress liquid (0)?
        s = s + '4\n'  #  4. Execute (follow path, mineral isograd or melt contour)
        s = s + '0\n'  #  Superliquidus (1) or subsolidus (0) initial guess ?
        # Phase to include (by name, lower case; 'x' when done):
        for ph in self.solution_phases:
            s = s + ph + '\n'
        s = s + 'x\n'

        return s
