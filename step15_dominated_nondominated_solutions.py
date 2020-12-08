# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:23:34 2019

@author: hongl
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:12:32 2019

@author: hongl
"""

# quickly plot Sx CDF
import numpy as np
from pylab import *
import matplotlib.pyplot as plt 
import os
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

attrib_raw = 'C:/Users/hongl/Documents/2019-06-21Discretization/optimize/source/model/input/step1_raw_data_us_asp180_no_round.csv'
ostrich_output1 = 'C:/Users/hongl/Documents/2019-06-21Discretization/optimize/slp_asp_7param_Sx_trial4/OstOutput0.txt'
ostrich_output2 = 'C:/Users/hongl/Documents/2019-06-21Discretization/optimize/slp_asp_7param_Sx_trial4/OstNonDomSolutions0.txt'

opath =  'C:/Users/hongl/Documents/2019-06-21Discretization/scripts_metrics/step15dominated_nondominatede_solutions'
ofile_obj = 'step15_dominated_nondominatede_solutions.png'

dpi_value = 100
if not os.path.exists(opath):
    os.makedirs(opath)

# dominated solutions
with open(ostrich_output1, 'r') as f:
    for iline, line in enumerate(f):
        line = line.strip()
        if line and line.startswith('Dominated Solutions'):
            line_start = iline+2
        if line and line.startswith('Number of Non-Dominated Solutions'):
            line_end = iline-2
            break
dom_slt = np.zeros((line_end-line_start+1, 10))
count = 0
with open(ostrich_output1, 'r') as f:
    for iline, line in enumerate(f):
        line = line.strip()
        if iline>= line_start and iline<=line_end:
            print(iline)
            slt_i = map(float, line.split())
            dom_slt[count,:]=np.asarray(slt_i)
            count = count+1
        if iline>line_end:
            break

nhru_dom = dom_slt[:,0]
sw_err_dom = dom_slt[:,1]*100
sx_err_dom = dom_slt[:,2]*100

# non-dominated solutions
diag = np.loadtxt(ostrich_output2, skiprows=3)
nhru_nondom = diag[:, 1]
sw_err_nondom = diag[:, 2]*100
sx_err_nondom = diag[:, 3]*100


# Plot
fig = figure(figsize=(7.48*0.7, 7.48*0.7))

ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(sw_err_nondom, sx_err_nondom, nhru_nondom, marker='o', s=8, c= 'r')
ax1.scatter(sw_err_dom, sx_err_dom, nhru_dom, marker='o', s=8, c= 'k')

ax1.set_xlabel('Sw error (%)')
ax1.set_ylabel('Sx error (%)')
ax1.set_zlabel('HRU number')

ax1.set_title('Objective functions')
#ax1.set_xlim(0, 10)
#ax1.set_ylim(0, 10)
fig.tight_layout()
fig.savefig(os.path.join(opath, ofile_obj), dpi=dpi_value)
plt.show()
#plt.close(fig) 
#
