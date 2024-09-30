import numpy as np
import matplotlib.pyplot as pl
from pybemt.solver import Solver

# Run sweep of tip speed ratio with the BEMT method
input_name = 'NREL5MW'
s = Solver(input_name+".ini",input_name)
low = 0
high = 0
interval = 0.1
calc_num = int((high-low)/(interval)+1)

df, section_df = s.run_sweep('twist', calc_num, low, high)
np.savetxt('output/'+input_name+'-df.txt', df)
np.savetxt('output/'+input_name+'-section.txt', section_df[0], delimiter='\t', header="radius\tchord\t\
pitch\tCl\tCd\tCl_Cd\tCn\tCt\tPn\tPt\tdT\tdQ\tF\ta\tap\tRe\tv\tv_theta\tv_rel\talpha\tphi")
s.OutputBlade(df, section_df, input_name, "output/"+input_name+"-Airfoil.txt")

