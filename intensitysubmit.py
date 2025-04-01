import os
import numpy as np

#runnumbers=np.genfromtxt('hess_all.txt')[:,0]
runnumbers=np.array([157341.])
print(runnumbers)

for runnumber in runnumbers:
    os.system('python submit_short.py \"python processor_intensity.py '+str(int(runnumber))+' /lfs/l7/hess/users/spencers/realdata\"')

