import os 
import sys
import numpy as np
import h5py

#Idea: .h5 file with various trail properties,
#i.e.: 
#duration, 
#unix time of beginning of trail 
# => convert to time in night
#average brightness
#maximum brightness (+minimum brightness?)
#zenith, azimuth
#width
#N pixels
#Distance on cam of average first and average last pixel in cam
#Angle of direction on cam (calc by using vectors)
 

#To parallelize the above job: write all these properties out to txt files while searching for the trails themselves.
#This script just takes all those files together and writes them into one big file of all properties
home_path = "/home/tholang"
path_to_files = home_path+"/trail_properties/"
filenames = sorted(os.listdir(path_to_files))
#change to hdf5 output
with open(path_to_files+"run_properties.txt", 'w') as outfile:
    outfile.write("nrun, zenith, azimuth, time-in-night, N-track, "+
                  "t0, start, duration, unique pix, width, velocity, avg-brightness, max-brightness \n")
    for fname in filenames:
        with open(path_to_files+fname) as infile:
            outfile.write(infile.read())

nruns = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 0, skiprows=1, delimiter = ",", dtype = int)
zenith = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 1, skiprows=1, delimiter = ",", dtype = float)
azimuth = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 2, skiprows=1, delimiter = ",", dtype = float)
time_in_night = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 3, skiprows=1, delimiter = ",", dtype = int)
N_track = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 4, skiprows=1, delimiter = ",", dtype = float)
t0 = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 5, skiprows=1, delimiter = ",", dtype = float)
start = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 6, skiprows=1, delimiter = ",", dtype = float)
duration = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 7, skiprows=1, delimiter = ",", dtype = float)
unique_pix = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 8, skiprows=1, delimiter = ",", dtype = int)
width = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 9, skiprows=1, delimiter = ",", dtype = float)
velocity = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 10, skiprows=1, delimiter = ",", dtype = float)
avg_brightness = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 11, skiprows=1, delimiter = ",", dtype = float)
max_brightness = np.loadtxt("/home/tholang/trail_properties/run_properties.txt", usecols = 12, skiprows=1, delimiter = ",", dtype = float)

if os.path.exists(path_to_files+"run_properties.h5"):
  os.remove(path_to_files+"run_properties.h5")
hfile = h5py.File(path_to_files+"run_properties.h5","w")
hfile.create_dataset("nrun", data = nruns, compression = "gzip", compression_opts=9)
hfile.create_dataset("zenith", data = zenith, compression = "gzip", compression_opts=9)
#hfile.create_dataset("azimuth", data = azimuth, compression = "gzip", compression_opts=9)
hfile.create_dataset("time_in_night", data = time_in_night, compression = "gzip", compression_opts=9)
hfile.create_dataset("N_track", data = N_track, compression = "gzip", compression_opts=9)
hfile.create_dataset("t0", data = t0, compression = "gzip", compression_opts=9)
hfile.create_dataset("start", data = start, compression = "gzip", compression_opts=9)
hfile.create_dataset("duration", data = duration, compression = "gzip", compression_opts=9)
hfile.create_dataset("unique_pix", data = unique_pix, compression = "gzip", compression_opts=9)
#hfile.create_dataset("width", data = width, compression = "gzip", compression_opts=9)
hfile.create_dataset("velocity", data = velocity, compression = "gzip", compression_opts=9)
hfile.create_dataset("avg_brightness", data = avg_brightness, compression = "gzip", compression_opts=9)
hfile.create_dataset("max_brightness", data = max_brightness, compression = "gzip", compression_opts=9)
hfile.close()
os.remove(path_to_files+"run_properties.txt")
