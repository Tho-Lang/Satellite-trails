# haptools imports
from haptools.coordinate_lookup import CoordinateLookup
from haptools.runselector import RunList, RunSelector
from haptools.run_quality_criteria import (HESS1_DETECTION_CRITERIA, HESS1_SPECTRAL_CRITERIA,
                                           HESS2_DETECTION_CRITERIA, HESS2_SPECTRAL_CRITERIA)
from haptools import utils
import os
import numpy as np
import h5py
from astropy.time import Time
import astropy.units as u

# Create an instance of the CoordinateLookup class
coordinate_lookup = CoordinateLookup()
#LMC and eta carinea coordinates
RA_LMC = 80.8942
Dec_LMC = -69.7561
RA_eta_car = 161.206
Dec_eta_car = -59.704

rs2 = RunSelector(hess2=True) 
#Remove every entry of the following two from all_runs
LMC_hess2 = rs2.select_runs_by_position(RA_LMC, Dec_LMC, ring_radius_max = 5)
eta_car_hess2 = rs2.select_runs_by_position(RA_eta_car, Dec_eta_car, ring_radius_max = 5)

#all_runs is complete runselector.runlist
all_runs = rs2.select_runs_by_position(83.6333, 22.0144, ring_radius_max=3000) 
rs2.query_required_data(HESS2_SPECTRAL_CRITERIA, runs=all_runs)
all_runs = rs2.select_runs_by_quality(HESS2_SPECTRAL_CRITERIA, all_runs, require_CT5=True, mintels = 5)
all_runs_run = np.array(rs2.data.Monitor_Run_Data.loc[all_runs.index]["Run"])
all_runs_dur = np.array(rs2.data.Monitor_Run_Data.loc[all_runs.index]["Duration"])
all_runs_pat = np.array(rs2.data.Monitor_Run_Data.loc[all_runs.index]["Telescope_Pattern"])
all_runs_runtype = np.array(rs2.data.Monitor_Run_Data.loc[all_runs.index]["RunType"])
mask_runs = np.logical_and(np.logical_and(np.logical_and(all_runs_dur >600,
                                                         all_runs_pat == 62),
                                          all_runs_runtype == "ObservationRun"),
                           np.logical_and(all_runs_run > 155000, all_runs_run <=999999))#177585)) #value used for previous comparison


runlist = all_runs_run[mask_runs]
#Remove every entry of the following two from all_runs
LMC_hess2 = rs2.select_runs_by_position(RA_LMC, Dec_LMC, ring_radius_max = 5)
eta_car_hess2 = rs2.select_runs_by_position(RA_eta_car, Dec_eta_car, ring_radius_max = 5)
runlist_red = np.array([nrun for nrun in runlist if nrun not in LMC_hess2 if nrun not in eta_car_hess2])

#Find Runs where Quality_Flag is "Good" 
good_runs = np.where(np.array(rs2.data.Monitor_Run_Trigger.Array["Quality_Flag"][runlist_red])== "Good")
runlist_red = runlist_red[good_runs]

#mean zenith angles
zen_runlist_red = np.array(rs2.data.Monitor_Run_Trigger.Array["Mean_Zenith"][runlist_red])
start_time = Time(np.array(rs2.data.Monitor_Run_Data.loc[runlist_red]["Run_Start_Time"]))
start_time.format = "unix"
start_time = start_time.value

if os.path.exists("hessall.h5"):
  os.remove("hessall.h5")
hfile = h5py.File("hessall.h5", "w")
hfile.create_dataset("run", data = runlist_red)
hfile.create_dataset("zenith", data = zen_runlist_red.astype(float))
hfile.create_dataset("start_time", data =start_time)
hfile.close()
