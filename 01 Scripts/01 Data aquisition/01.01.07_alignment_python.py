from prody import *
from pylab import *
import os
ion()

#inputs
directory = 'C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\20231211_KS_dimers'
input_list = os.listdir(directory)
Template = parsePDB('C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\20231211_KS_dimers\\abyB3_Mod.1.pdb')
pdb2 = parsePDB('C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\20231211_KS_dimers\\abyB1_Mod.1.pdb')

Template.numCoordsets()
showProtein(Template)
showProtein(pdb2);
legend()

for i in input_list:
    pdb2 = parsePDB([directory + "\\" + i])
    matchAlign(target = Template, mobile = pdb2, seqid = 1, overlap = 1)
    showProtein(Template)
    showProtein(pdb2)
    writePDB(filename = f'C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\20231211_KS_dimers_aligned\\{i}_aligned', atoms= pdb2)
    legend()
    #savefig(f'C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\Aligned\\figs\\{i}.', )

    figure(figsize = (15,15))