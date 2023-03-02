import sys
import os
import numpy as np
import h5py
from astropy.time import Time
from astropy import units as u

def single_selection_cut(data_high_pixel):
    max_counts = 100
    values, counts = np.unique(data_high_pixel[0], return_counts=True)
    mask_max_counts = counts<max_counts
    st = set(values[mask_max_counts])
    result = [i for i, e in enumerate(data_high_pixel[0]) if e in st]
    new_pix = np.array(data_high_pixel[0])[result]
    new_time = np.array(data_high_pixel[1])[result]
    new_brightness = np.array(data_high_pixel[2])[result]
    high_pixel_cut = np.array([new_pix, new_time, new_brightness])
    return high_pixel_cut

path_high_pixel = "/home/tholang/eval_high_data_utc/high_pixel"
path_high_pixel_cut = "/home/tholang/eval_high_data_utc_cut/high_pixel"
if len(sys.argv)<2:
    print("need number as entry")
    exit()

nrun = int(sys.argv[1])
print("Convert "+str(nrun)+"...")
subpath_to_run = "/run"+str(int(nrun-nrun%200))+"-"+str(int(nrun-nrun%200+199))+"/"+str(int(nrun))+"/"
path_to_run = path_high_pixel + subpath_to_run
path_to_run_cut = path_high_pixel_cut + subpath_to_run

try:
    txt_files = os.listdir(path_to_run)
    print(txt_files)
except:
    print("os.listdir does not work with this run")
    exit()
if len(txt_files) ==0:
    print("No Dataset found for this run")
    exit()

if "high_pix_"+str(nrun)+"_CT_5.h5" in txt_files:
    print("high_pix_"+str(nrun)+"_CT_5.h5 already exists!")
    os.remove(path_to_run + "high_pix_"+str(nrun)+"_CT_5.txt")
    exit()

pix = np.loadtxt(path_to_run+"high_pix_"+str(int(nrun))+"_CT_5.txt",
                 usecols = 0, delimiter = ";", ndmin = 1)
brightness = np.loadtxt(path_to_run+"high_pix_"+str(int(nrun))+"_CT_5.txt",
                        usecols = 1, delimiter = ";", ndmin = 1)
time = np.loadtxt(path_to_run+"high_pix_"+str(int(nrun))+"_CT_5.txt",
                  usecols = 2, delimiter = ";", dtype = np.str, ndmin = 1)
times = [item.strip("UTC: ") for item in time]
t = Time(times, format = "iso", scale = "utc")
t.format ="unix"
t = t.value
os.remove(path_to_run + "high_pix_"+str(nrun)+"_CT_5.txt")


data = [pix, t, brightness]
data = single_selection_cut(data)

if len(data[0])==0:
    data = [[-1], [-1], [-1]]

if not os.path.exists(path_to_run_cut+str(nrun)):
    os.makedirs(path_to_run_cut+str(nrun))
hfile = h5py.File(path_to_run_cut + "high_pix_"+str(nrun)+"_CT_5_cut.h5", "w")
hfile.create_dataset("Pix ID", data = data[0], compression = "gzip", compression_opts=9)
hfile.create_dataset("Time", data = data[1], compression = "gzip", compression_opts=9)
hfile.create_dataset("Brightness", data = data[2], compression = "gzip", compression_opts=9)
hfile.close()
os.remove(path_to_run + "high_pix_"+str(nrun)+"_CT_5.txt")
