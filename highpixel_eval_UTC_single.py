import os 
import numpy
import sys

nrun = int(sys.argv[1])
os.system('root -l -b -q highpixel_eval_UTC_single.C+\('+str(int(nrun))+'\)')
os.system('python scripts/convert_txt_to_hdf5_and_selection_cut.py '+str(int(nrun)))
