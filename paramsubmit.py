import os
import numpy as np

runnumbers=np.genfromtxt('0940_3deg_samcheck_postjira121.lis')[:,0]
#runnumbers=np.array([180699,180702,181025])
#runnumbers=np.array([160723,160736])
print(runnumbers)

for runnumber in runnumbers:
    os.system('python submit_short.py \"python processor_params.py '+str(int(runnumber))+' /lfs/l7/hess/users/spencers/realdata/0940params\"')

