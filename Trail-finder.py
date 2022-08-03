#!/usr/bin/env python
# coding: utf-8

# In[34]:


import numpy as np
import pandas as pd
from scipy import ndimage
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tables import *
import os
import operator
import sys
import csv
import time
from astropy.time import Time
import astropy.units as u


# In[35]:


###Initialisation of nruns and dictionaries:
nruns = []
trails_high_avg = {}
avg_pix = {}
avg_tmp = {}
delta_pixel = {}
high_pixel = {}
roll_avg = {}
zenith = {}
azimuth = {}
broken_pixel = {}

path = {}


# In[98]:


#appending  nrun to nruns 
os.chdir("D:\\Masterarbeit ECAP\\First plots\\eval_data\\high_pixel")
first_run = 145000
runs_directory = "run"+str(first_run)+"-"+str(first_run+199)
if (first_run % 200==0):
    pass
else:
    print("First run not divisible by 200")
    sys.exit()
os.chdir(runs_directory)
files = os.listdir(os.curdir)
nruns = list(nruns)
new_nruns = []
for i in range(len(files)):
    try:
        new_nruns.append(int(files[i]))
        nruns.append(int(files[i]))
    except:
        print("No runs in this folder")
nruns = sorted(np.unique(nruns))
print(new_nruns)
print(nruns)
os.chdir("D:\\Masterarbeit ECAP\\First plots\\eval_data")
print(os.getcwd())


# In[92]:


get_ipython().run_cell_magic('time', '', '#delete all empty files\nfor nrun in new_nruns:\n    avg_pix[nrun] = {}\n    avg_tmp[nrun] = {}\n    delta_pixel[nrun] = {}\n    high_pixel[nrun] = {}\n    roll_avg[nrun] = {}\n    broken_pixel[nrun] = {}\n    pix_path = "avg_per_pix\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    tmp_path = "avg_per_tmp\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    delta_path = "delta_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    high_path = "high_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    roll_avg_path = "roll_avg\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    broken_pixel_path = "broken_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    path[nrun] = {}\n    path[nrun]["pix"] = "avg_per_pix\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    path[nrun]["tmp"] = "avg_per_tmp\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    path[nrun]["delta"] = "delta_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    path[nrun]["high"] = "high_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    path[nrun]["roll_avg"] = "roll_avg\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    path[nrun]["broken_pixel"] = "broken_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    \n    for key in path[nrun]:\n        try:\n            file_path = path[nrun][key]\n            #print(file_path)\n            for i in range(len(file_path)):\n                file_name = os.listdir(path[nrun][key])[i]\n                #print(file_name)\n                if os.path.getsize(file_path+file_name) == 0:\n                    print("Empty file", file_path+file_name ,"is being removed")\n                    os.remove(file_path+file_name)\n        except:\n            pass')


# In[103]:


get_ipython().run_cell_magic('time', '', '#Reading data of new_nruns, not overwriting the old ones\nfor nrun in new_nruns:\n    avg_pix[nrun] = {}\n    avg_tmp[nrun] = {}\n    delta_pixel[nrun] = {}\n    high_pixel[nrun] = {}\n    roll_avg[nrun] = {}\n    broken_pixel[nrun] = {}\n    telId=1\n    pix_path = "avg_per_pix\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    tmp_path = "avg_per_tmp\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    delta_path = "delta_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    high_path = "high_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    roll_avg_path = "roll_avg\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n    broken_pixel_path = "broken_pixel\\\\"+runs_directory+"\\\\"+str(nrun)+"\\\\"\n\n    print(pix_path)\n    for telId in range(1,6):\n        try:\n            print("Reading delta_pix_"+str(nrun)+"_CT_"+str(telId)+".txt")    \n            delta_pixel[nrun]["pix"+str(telId)] = np.loadtxt(delta_path + "delta_pix_"+str(nrun)+"_CT_"+str(telId)+".txt",\n                                                             usecols = 0, delimiter = ";")\n            delta_pixel[nrun][telId] = np.loadtxt(delta_path + "delta_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                             usecols = 1, delimiter = ";")\n            delta_pixel[nrun]["tmp"+str(telId)] = np.loadtxt(delta_path + "delta_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                             usecols = 2, delimiter = ";")\n            delta_pixel[nrun]["x-pos"+str(telId)] = np.loadtxt(delta_path + "delta_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                             usecols = 3, delimiter = ";")\n            delta_pixel[nrun]["y-pos"+str(telId)] = np.loadtxt(delta_path + "delta_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                             usecols = 4, delimiter = ";")\n        except:\n             print("No file named delta_pix_"+str(nrun)+"_CT_"+str(telId)+".txt found")\n\n        try:\n            print("Reading broken_pixel"+str(nrun)+"_CT_"+str(telId)+".txt")        \n            broken_pixel[nrun]["pix"+str(telId)] = np.loadtxt(broken_pixel_path + "broken_pixel"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 0, delimiter = ";")\n            broken_pixel[nrun]["tmp"+str(telId)] = np.loadtxt(broken_pixel_path + "broken_pixel"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 1, delimiter = ";")\n            broken_pixel[nrun]["x-pos"+str(telId)] = np.loadtxt(broken_pixel_path + "broken_pixel"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 2, delimiter = ";")\n            broken_pixel[nrun]["y-pos"+str(telId)] = np.loadtxt(broken_pixel_path + "broken_pixel"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 3, delimiter = ";")\n            broken_pixel[nrun][telId] = np.loadtxt(broken_pixel_path + "broken_pixel"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 4, delimiter = ";")\n        except:\n            print("No file named broken_pixel"+str(nrun)+"_CT_"+str(telId)+".txt found")\n            \n        try:\n            print("Reading high_pix_"+str(nrun)+"_CT_"+str(telId)+".txt")        \n            high_pixel[nrun]["pix"+str(telId)] = np.loadtxt(high_path + "high_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 0, delimiter = ";")\n            high_pixel[nrun][telId] = np.loadtxt(high_path + "high_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 1, delimiter = ";")\n            high_pixel[nrun]["tmp"+str(telId)] = np.loadtxt(high_path + "high_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 2, delimiter = ";")\n            high_pixel[nrun]["x-pos"+str(telId)] = np.loadtxt(high_path + "high_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 3, delimiter = ";")\n            high_pixel[nrun]["y-pos"+str(telId)] = np.loadtxt(high_path + "high_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 4, delimiter = ";")\n            \n        except:\n            print("No file named high_pix_"+str(nrun)+"_CT_"+str(telId)+".txt found")\n            \n        \n        try:\n            print("Reading roll_avg_"+str(nrun)+"_CT_"+str(telId)+".txt")      \n            roll_avg[nrun]["pix"+str(telId)] = np.loadtxt(roll_avg_path + "roll_avg_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 0, delimiter = ";")\n            roll_avg[nrun][telId] = np.loadtxt(roll_avg_path + "roll_avg_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 1, delimiter = ";")\n            roll_avg[nrun]["tmp"+str(telId)] = np.loadtxt(roll_avg_path + "roll_avg_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 2, delimiter = ";")\n            roll_avg[nrun]["x-pos"+str(telId)] = np.loadtxt(roll_avg_path + "roll_avg_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 3, delimiter = ";")\n            roll_avg[nrun]["y-pos"+str(telId)] = np.loadtxt(roll_avg_path + "roll_avg_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                            usecols = 4, delimiter = ";")  \n        except:\n            print("No file named roll_avg_"+str(nrun)+"_CT_"+str(telId)+".txt found")\n        \n        \n\n        try:\n            print("Reading avg_per_pix_"+str(nrun)+"_CT_"+str(telId)+".txt")        \n            avg_pix[nrun][telId] = np.loadtxt(pix_path + "avg_per_pix_"+str(nrun)+"_CT_"+str(telId)+".txt", usecols = 1)\n        except:\n            print("No file named avg_per_pix_"+str(nrun)+"_CT_"+str(telId)+".txt found")\n\n        try:\n            print("Reading avg_per_tmp_"+str(nrun)+"_CT_"+str(telId)+".txt")\n            avg_tmp[nrun][telId] = np.loadtxt(tmp_path +"avg_per_tmp_"+str(nrun)+"_CT_"+str(telId)+".txt", \n                                                         usecols = 1, delimiter = ";")\n            avg_tmp[nrun]["tmp"+str(telId)] = np.loadtxt(tmp_path +"avg_per_tmp_"+str(nrun)+"_CT_"+str(telId)+".txt",\n                                                         usecols = 0, delimiter = ";")\n           \n        except:\n            print("No file named avg_per_tmp_"+str(nrun)+"_CT_"+str(telId)+".txt found")\n ')


# In[64]:


#Define camera geometry:
fin_data = open_file("D:/Masterarbeit ECAP/2022-04-06/Thomas/gamma_sat_postselect.h5", mode="r")
geom_hess1_xc_from_root = np.squeeze(fin_data.get_node('/configuration/instrument/telescope/camera/geometry_0').col('pix_x'))
geom_hess1_yc_from_root = np.squeeze(fin_data.get_node('/configuration/instrument/telescope/camera/geometry_0').col('pix_y'))
geom_hess5_xc_from_root = np.squeeze(fin_data.get_node('/configuration/instrument/telescope/camera/geometry_1').col('pix_x'))
geom_hess5_yc_from_root = np.squeeze(fin_data.get_node('/configuration/instrument/telescope/camera/geometry_1').col('pix_y'))
a = sorted(geom_hess5_xc_from_root)
b = sorted(geom_hess5_yc_from_root)
for i in range(len(a)-1):
    if a[i+1]-a[i]>0:
        min_x_dif = round(a[i+1]-a[i], 5)
        print(min_x_dif)
        break
for i in range(len(b)-1):
    if b[i+1]-b[i]>0:
        min_y_dif = round(b[i+1]-b[i], 5)
        print(min_y_dif)
        break


# In[65]:


###activate/deactivate print statements (Work in progress):
temp_stdout = None

# Disable
def blockPrint():
    global temp_stdout
    temp_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    global temp_stdout
    sys.stdout = temp_stdout


# In[ ]:





# In[66]:


def lin_fit(x, m, t):
    return m*x+t
files = os.listdir(os.curdir)
def plot_scatter(x, y, z):
    #fig= plt.figure(figsize=(6,5), )
    plt.scatter(x,y,c=z, cmap = "jet")
    plt.xlim(-1.25,1.25)
    plt.ylim(-1.25,1.25)
    plt.xlabel("Camera x-pos")
    plt.ylabel("Camera y-pos")
    plt.title("Run {}".format(nrun)+", CT {}".format(telId))
    cbar = plt.colorbar()
    cbar.set_label("Time  in run [s]")

def get_next_neighbours(dict_x, dict_y, dict_z, width):
    nnn_x = []
    nnn_y = []
    nnn_z = []
    for k in range(len(dict_x)):
        x_x_cut = []
        y_x_cut = []
        z_x_cut = []
        x_y_cut = []
        y_y_cut = []
        z_y_cut = []
        x = dict_x[k]
        y = dict_y[k]
        z = dict_z[k]
        for i in range(len(dict_x)):
            if x-width*min_x_dif < dict_x[i] < x+width*min_x_dif:
                x_x_cut.append(dict_x[i])
                y_x_cut.append(dict_y[i])
                z_x_cut.append(dict_z[i])

        for i in range(len(x_x_cut)):
            if y-width*min_y_dif < y_x_cut[i] < y+width*min_y_dif:
                x_y_cut.append(x_x_cut[i])
                y_y_cut.append(y_x_cut[i])
                z_y_cut.append(z_x_cut[i])
        nnn_x.append(x_y_cut)
        nnn_y.append(y_y_cut)
        nnn_z.append(z_y_cut)
    return nnn_x, nnn_y, nnn_z
def check_z_neighbours(dict_x, dict_y, dict_z):
    nnn_x = []
    nnn_y = []
    nnn_z = []
    for i in range(len(dict_z)):
        close_times = np.concatenate((np.where(np.array(dict_z) == dict_z[i]-1)[0],
                                    np.where(np.array(dict_z) == dict_z[i])[0],
                                    np.where(np.array(dict_z) == dict_z[i]+1)[0]))
        nnn_x.append(np.array(dict_x)[close_times])
        nnn_y.append(np.array(dict_y)[close_times])
        nnn_z.append(np.array(dict_z)[close_times])
    return nnn_x, nnn_y, nnn_z

def cleaning_cut(dict_x,dict_y,dict_z, dict_brightness, width, length):
    x_reduced = []
    y_reduced = []
    z_reduced = []
    brightness_reduced = []
    nn_x, nn_y, nn_z = get_next_neighbours(dict_x, dict_y, dict_z, width)
    for k in range(len(dict_x)):
        if len(np.unique(nn_y[k]))>length:
            x_reduced.append(dict_x[k])
            y_reduced.append(dict_y[k])
            z_reduced.append(dict_z[k])
            brightness_reduced.append(dict_brightness[k])
    return x_reduced, y_reduced, z_reduced, brightness_reduced
def cleaning_z_cut(dict_x,dict_y,dict_z, width):
    nnn_x, nnn_y, nnn_z = get_next_neighbours(dict_x, dict_y, dict_z, width)
    for i in range(len(dict_z)):
        counter = 0
        for j in range(len(nnn_z[i])):
            if dict_z[i]-2 < nnn_z[i][j] < dict_z[i]+2:
                    counter = counter + 1
        if counter < 7: 
            print("less than 7 values in nn_z for index", i)
def sort_into_tracks(dict_x,dict_y,dict_z,dict_brightness, width):
    nnn_x, nnn_y, nnn_z = get_next_neighbours(dict_x, dict_y, dict_z, width)
    tracks = {}
    max_tracks = 9
    for N in range(1,max_tracks):
        tracks["x-"+str(N)] = []
        tracks["y-"+str(N)] = []
        tracks["z-"+str(N)] = []
        tracks["brightness-"+str(N)] = []
    k=0
    for i in range(len(dict_x)):
        appended_to_track = False
        empty_first_track = True
        if k == 0:
            tracks["x-"+str(1)].append(dict_x[i])
            tracks["y-"+str(1)].append(dict_y[i])
            tracks["z-"+str(1)].append(dict_z[i])
            tracks["brightness-"+str(1)].append(dict_brightness[i])
            appended_to_track  = True 
            k = k+1
            continue
        for N in range(1,max_tracks):
            if len(tracks["x-"+str(N)])==0:
                continue            
            if tracks["x-"+str(N)][-1] in nnn_x[i] and tracks["y-"+str(N)][-1] in nnn_y[i]:
                if dict_z[i] -2< tracks["z-"+str(N)][-1]:
                    tracks["x-"+str(N)].append(dict_x[i])
                    tracks["y-"+str(N)].append(dict_y[i])
                    tracks["z-"+str(N)].append(dict_z[i])
                    tracks["brightness-"+str(N)].append(dict_brightness[i])
                    appended_to_track  = True 
                    break
                else:
                    if len(tracks["x-"+str(N)])<3:
                        tracks["x-"+str(N)] = []
                        tracks["y-"+str(N)] = []
                        tracks["z-"+str(N)] = []
                        tracks["brightness-"+str(N)] = []
                        print("track", N, " was deleted at time", dict_z[i] )
                        
        if appended_to_track == False:
            for N in range(1,max_tracks):
                if len(tracks["x-"+str(N)])==0:
                    tracks["x-"+str(N)].append(dict_x[i])
                    tracks["y-"+str(N)].append(dict_y[i])
                    tracks["z-"+str(N)].append(dict_z[i])
                    tracks["brightness-"+str(N)].append(dict_brightness[i])
                    break
    possible_meteorites = {}
    for N in range(1,max_tracks):
        if len(tracks["x-"+str(N)]) != 0 and len(np.unique(tracks["z-"+str(N)]))<2:
            print("Track", N, "from", np.unique(tracks["z-"+str(N)]),
                  ",too fast", len(tracks["z-"+str(N)]), " entires, hence no track")
            possible_meteorites["x-"+str(N)] = tracks["x-"+str(N)] 
            possible_meteorites["y-"+str(N)] = tracks["y-"+str(N)] 
            possible_meteorites["z-"+str(N)] = tracks["z-"+str(N)] 
            possible_meteorites["brightness-"+str(N)] = tracks["brightness-"+str(N)]
            tracks["x-"+str(N)] = []
            tracks["y-"+str(N)] = []
            tracks["z-"+str(N)] = []
            tracks["brightness-"+str(N)] = []
        try:
            if possible_meteorites["z-"+str(N)]:
                pass # add tracks with same time
        except:
            pass
        
            
    return tracks, possible_meteorites



def sort_into_tracks_old(dict_x,dict_y,dict_z, width):
    nnn_x, nnn_y, nnn_z = get_next_neighbours(dict_x, dict_y, dict_z, width)
    tracks = {}
    for N in range(1,6):
        tracks["x-"+str(N)] = []
        tracks["y-"+str(N)] = []
        tracks["z-"+str(N)] = []
    k=0
    for i in range(len(dict_x)):
        appended_to_track = False
        empty_first_track = True
        if k == 0:
            tracks["x-"+str(1)].append(dict_x[i])
            tracks["y-"+str(1)].append(dict_y[i])
            tracks["z-"+str(1)].append(dict_z[i])
            appended_to_track  = True 
            k = k+1
            continue
        for N in range(1,6):
            if len(tracks["x-"+str(N)])==0:
                continue            
            if tracks["x-"+str(N)][-1] in nnn_x[i] and tracks["y-"+str(N)][-1] in nnn_y[i]:
                if dict_z[i] -2< tracks["z-"+str(N)][-1]:
                    tracks["x-"+str(N)].append(dict_x[i])
                    tracks["y-"+str(N)].append(dict_y[i])
                    tracks["z-"+str(N)].append(dict_z[i])
                    appended_to_track  = True 
                    break
                else:
                    if len(tracks["x-"+str(N)])<3:
                        tracks["x-"+str(N)] = []
                        tracks["y-"+str(N)] = []
                        tracks["z-"+str(N)] = []
                        print("track", N, " was deleted at time", dict_z[i] )
                        
        if appended_to_track == False:
            for N in range(1,6):
                if len(tracks["x-"+str(N)])==0:
                    tracks["x-"+str(N)].append(dict_x[i])
                    tracks["y-"+str(N)].append(dict_y[i])
                    tracks["z-"+str(N)].append(dict_z[i])
                    break
    for N in range(1,6):
        if len(tracks["x-"+str(N)]) != 0 and len(np.unique(tracks["z-"+str(N)]))<2:
            print("Track", N, "from", np.unique(tracks["z-"+str(N)]),
                  ",too fast", len(tracks["z-"+str(N)]), " entires, hence no track")
            tracks["x-"+str(N)] = []
            tracks["y-"+str(N)] = []
            tracks["z-"+str(N)] = []
            
    return tracks

def get_track_width(dict_x, dict_y):
    popt, pcov = curve_fit(lin_fit, dict_x, dict_y)
    d_from_lin_fit = []
    for i in range(len(dict_x)):
        x_p = dict_x[i]
        y_p = dict_y[i]
        t_2 = x_p/popt[0]+y_p
        x_s = popt[0]/(popt[0]**2+1)*(t_2 -popt[1])
        y_s = popt[0]*x_s+popt[1]
        d_from_lin_fit.append(np.sqrt(((x_p-x_s)/min_x_dif)**2+((y_p-y_s)/min_y_dif)**2))
    print("Average distance from linfit: ",  round(np.average(d_from_lin_fit), 4))
    print("Average width = ", round(np.average(d_from_lin_fit)*2, 4), "Pixels")
    


# In[215]:


def draw_box_new(m, t):
    mask_x_top = np.around(np.array(geom_hess5_xc_from_root), 6) == np.max(np.around(geom_hess5_xc_from_root, 6))
    mask_x_bot = np.around(np.array(geom_hess5_xc_from_root), 6) == np.min(np.around(geom_hess5_xc_from_root, 6))
    mask_y_top = np.around(np.array(geom_hess5_yc_from_root), 6) == np.max(np.around(geom_hess5_yc_from_root, 6))
    mask_y_bot = np.around(np.array(geom_hess5_yc_from_root), 6) == np.min(np.around(geom_hess5_yc_from_root, 6))
    
    y_secondtolast = sorted(np.unique(np.around(geom_hess5_yc_from_root, 6)))[-2]
    y_secondtofirst = sorted(np.unique(np.around(geom_hess5_yc_from_root, 6)))[1]
    mask1 = np.around(np.array(geom_hess5_yc_from_root), 6) == y_secondtolast
    mask2 = np.around(np.array(geom_hess5_yc_from_root), 6) == y_secondtofirst
    
    box_x = [geom_hess5_xc_from_root[mask_x_top][0],
             np.max(np.around(geom_hess5_xc_from_root[mask1], 6)),
            np.max(geom_hess5_xc_from_root[mask_y_top]),
            np.min(geom_hess5_xc_from_root[mask_y_top]),
            geom_hess5_xc_from_root[mask_x_bot][0],
            geom_hess5_xc_from_root[mask_x_bot][-1],
            np.min(geom_hess5_xc_from_root[mask_y_bot]),
            np.max(geom_hess5_xc_from_root[mask_y_bot]),
             np.max(np.around(geom_hess5_xc_from_root[mask2], 6)),
            geom_hess5_xc_from_root[mask_x_top][0],
            geom_hess5_xc_from_root[mask_x_top][0]]
    box_y = [geom_hess5_yc_from_root[mask_x_top][0],
             y_secondtolast,
             geom_hess5_yc_from_root[mask_y_top][0],
             geom_hess5_yc_from_root[mask_y_top][-1],
             np.max(geom_hess5_yc_from_root[mask_x_bot]),
             np.min(geom_hess5_yc_from_root[mask_x_bot]),
             geom_hess5_yc_from_root[mask_y_bot][0],
             geom_hess5_yc_from_root[mask_y_bot][-1],
             y_secondtofirst,
             geom_hess5_yc_from_root[mask_x_top][0],
             geom_hess5_yc_from_root[mask_x_top][0]]
    intersections = []
    for i in range(len(box_x)-2):
        if (box_x[i+1]-box_x[i]) !=0:
            if (box_x[i+1]-box_x[i])> 0:
                x_vals = np.arange(box_x[i], box_x[i+1], 0.001)
            else:
                x_vals = np.arange(box_x[i+1], box_x[i], 0.001)
            m_box = (box_y[i+1]-box_y[i])/(box_x[i+1]-box_x[i])
            t_box = box_y[i]- m_box*box_x[i]
            y_vals = m_box*x_vals + t_box
            plt.plot(x_vals, y_vals, c = "black")
            idx = np.argwhere(np.diff(np.sign(m*x_vals+t - y_vals))).flatten()
#             if len(idx) !=0:
#                 intersections.append([x_vals[idx][0], y_vals[idx][0]])
#             plt.scatter(x_vals[idx], y_vals[idx], facecolors='none', edgecolors="black")
        else:
            x_vals = [box_x[i]]
            if box_y[i]<box_y[i+1]:
                y_vals = np.arange(box_y[i], box_y[i+1], 0.001)
            else:
                y_vals = np.arange(box_y[i+1], box_y[i], 0.001)
            for j in range(len(y_vals)-1):
                x_vals.append(x_vals[0])
            plt.vlines(box_x[i], box_y[i], box_y[i+1], colors ="black")
            idx = np.argwhere(np.diff(np.sign(m*np.array(x_vals)+t - y_vals))).flatten()
#             try:
#                 if len(idx) !=0:
#                     idx = idx[0]
#                     intersections.append([x_vals[idx][0], y_vals[idx][0]])
#             except:
#                 print(idx)
#     intersections = np.around(intersections, 4)
#     print("intersections are", intersections[0], ", ", intersections[1])
    return intersections   


# In[68]:


def get_nn_pix(size):
    nn_pix = []
    for pix in range(len(geom_hess5_xc_from_root)):
        mask_x = np.logical_and(np.around(geom_hess5_xc_from_root,4)<
                                min_x_dif*size*0.6+np.around(geom_hess5_xc_from_root,4)[pix],
                                np.around(geom_hess5_xc_from_root,4)>
                                -min_x_dif*size*0.6+np.around(geom_hess5_xc_from_root,4)[pix])

        mask_xy = np.logical_and(np.around(geom_hess5_yc_from_root,4)[mask_x]<
                                 min_y_dif*(size-0.1)+np.around(geom_hess5_yc_from_root,4)[pix],
                                 np.around(geom_hess5_yc_from_root,4)[mask_x]>
                                 -min_y_dif*size+np.around(geom_hess5_yc_from_root,4)[pix])


        st_x = set(np.around(geom_hess5_xc_from_root,4)[mask_x][mask_xy])
        st_y = set(np.around(geom_hess5_yc_from_root,4)[mask_x][mask_xy])
        indices_x = np.array([j for j, e in enumerate(np.around(geom_hess5_xc_from_root,4)) if e in st_x])
        indices_y = np.array([j for j, e in enumerate(np.around(geom_hess5_yc_from_root,4)) if e in st_y])

        indices_xy = np.array([j for j, e in enumerate(indices_x) if e in indices_y ])
        nn_pix.append(indices_x[indices_xy])
    return nn_pix
nn_pix = get_nn_pix(10)


# In[105]:


print(nn_pix[1090])


# In[69]:


def draw_box(m, t):
    y_top = np.where(geom_hess5_yc_from_root == np.max(geom_hess5_yc_from_root))[0]
    y_bot = np.where(geom_hess5_yc_from_root == np.min(geom_hess5_yc_from_root))[0]
    x_top = np.where(geom_hess5_xc_from_root == np.max(geom_hess5_xc_from_root))[0]
    x_bot = np.where(geom_hess5_xc_from_root == np.min(geom_hess5_xc_from_root))[0]

    #draw a box 
    box_x = []
    box_y = []
    box_x.append(sorted(geom_hess5_xc_from_root[y_top])[0])
    box_x.append(sorted(geom_hess5_xc_from_root[y_top])[-1])
    box_x.append(sorted(geom_hess5_xc_from_root[x_top])[0])
    box_x.append(sorted(geom_hess5_xc_from_root[y_bot])[-1])
    box_x.append(sorted(geom_hess5_xc_from_root[y_bot])[0])
    box_x.append(sorted(geom_hess5_xc_from_root[x_bot])[0])
    box_x.append(sorted(geom_hess5_xc_from_root[x_bot])[1])
    box_x.append(sorted(geom_hess5_xc_from_root[y_top])[0])

    box_y.append(sorted(geom_hess5_yc_from_root[y_top])[0])
    box_y.append(sorted(geom_hess5_yc_from_root[y_top])[-1])
    box_y.append(sorted(geom_hess5_yc_from_root[x_top])[0])
    box_y.append(sorted(geom_hess5_yc_from_root[y_bot])[-1])
    box_y.append(sorted(geom_hess5_yc_from_root[y_bot])[0])
    box_y.append(sorted(geom_hess5_yc_from_root[x_bot])[0])
    box_y.append(sorted(geom_hess5_yc_from_root[x_bot])[1])
    box_y.append(sorted(geom_hess5_yc_from_root[y_top])[0])
    # 7 lines to draw the box:    
    plt.plot(box_x, box_y, c = "blue")
    #plt.plot((box_x[-1], box_x[0]),(box_y[-1], box_y[0]), c = "blue")
    #find intersections with lin_fit:
    box_lines_param = {}
    for i in range(0,7):
        if (box_x[i+1]-box_x[i]) !=0:
            m_box = (box_y[i+1]-box_y[i])/(box_x[i+1]-box_x[i])
            t_box = box_y[i]- m_box*box_x[i]
            box_lines_param[i] =[m_box, t_box] 
        else:
            box_lines_param[i] = [box_x[i]]
#     m_box = (box_y[0]-box_y[-1])/(box_x[0]-box_x[-1])
#     t_box = box_y[-1]- m_box*box_x[-1]
#     box_lines_param[6] =[m_box, t_box] 
    intersections = []
    for i in range(len(box_lines_param)):
        if len(box_lines_param[i])==1:
            pass
        else: 
            x_sp = (box_lines_param[i][1] - t)/(m - box_lines_param[i][0])
        y_sp = m * x_sp + t
        if -1.25< x_sp < 1.25:
            if -1.25< y_sp < 1.25:
                intersections.append([x_sp, y_sp])
    intersections = np.array(intersections)
    plt.scatter(intersections[:,0],intersections[:,1], marker ="o", facecolors='none', edgecolors='black', s = 100)
    print("Coordinates of intersections:", np.around(intersections[0], 4), ", ", np.around(intersections[1], 4))
    return(intersections)
    


# In[70]:


def print_important_params(track_x, track_y, track_z):
    print("Track", N, end =": ")    
    start_time = track_z[0]
    stop_time = track_z[-1]
    print("Start:", int(start_time), "s , Stop: ", int(stop_time), "s , Duration: ", int(stop_time-start_time), "s" )
    popt_param, pcov_param = curve_fit(lin_fit, track_x, track_y)
    print("linear fit : y =",round(popt_param[0],4),"* x +",round(popt_param[1],4)) 
    x_y_len = np.sqrt((track_x[-1]-track_x[0])**2+
                      (lin_fit(np.array(track_x), *popt_param)[-1] - 
                       lin_fit(np.array(track_x), *popt_param)[0])**2 
                             )
    print("Length on cam:", round(x_y_len, 4), "m(?)")
    print("Speed on cam: ", round(x_y_len/dur,4), "m/s")
    x_y = []
    for i in range(len(track_x)):
        x_y.append((round(track_x[i],4), 
                    round(track_y[i],4)))
    print("Unique pixels:",len(list(set(x_y))))
    print("Pixels per unit of length:", round(len(list(set(x_y)))/x_y_len, 4) )   


# In[71]:


def conversion_horizontal_to_equatorial(zen, az, UTC):
    #RA = alpha, DEC = delta, Local hour angle = H = LST - RA = > RA = LST - H, latitude = phi, altitude = alt
    alt = 1800
    lon = -16.5
    lat = 23.2715
    LST = 
    
    DEC = np.arcsin(np.sin(alt)*np.sin(lat) + np.cos(alt)*np.cos(lat)*np.cos(az) )
    LHA = np.arcsin(-np.sin(az)*np.cos(alt)/np.cos(DEC) )
    RA = LST - LHA
    sin(δ) = sin(a)sin(φ) + cos(a) cos(φ) cos(A)
    sin(H) = - sin(A) cos(a) / cos(δ)
    cos(H) = { sin(a) - sin(δ) sin(φ)} / cos(δ) cos(φ)
    α = t – H


# In[85]:


def get_neighbour_pixel_entries(dict_x,dict_y,dict_z, dict_pix, dict_brightness):
    #check für jeden Eintrag, ob in den next_neighbours einträge innerhalb der nächsten 2s liegen
    #nn_pix ist array mit array von neighbours und next neigbours als einträgen (2d-array), nn_pix[N] gibt nn von pixId N an
    min_len = 0
    x = []
    y = []
    z = []
    pix = []
    brightness = []
    for i in range(len(dict_pix)):
        mask_time = np.logical_and( dict_z>dict_z[i]-3,dict_z<dict_z[i]+3)
        st = set(nn_pix[int(dict_pix[i])])
        mask_space = [i for i, e in enumerate(dict_pix[mask_time]) if e in st]
        
        if len(np.unique(dict_pix[mask_time][mask_space]))>min_len:
            x.append(dict_x[i])
            y.append(dict_y[i])
            z.append(dict_z[i])
            pix.append(dict_pix[i])
            brightness.append(dict_brightness[i])
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    pix = np.array(pix)
    brightness = np.array(brightness)
    
    time_values, time_counts = np.unique(z, return_counts=True)
    if len(z)>0:  
        mask_min_counts = time_counts>0
        set_min_counts =set(time_values[mask_min_counts])
        mask = [i for i, e in enumerate(z) if e in set_min_counts] 
        # mask is for all values in the previously appended arrays to get only the values which also fulfill mask_min_counts 
        plot_scatter(x[mask],y[mask],z[mask])
        draw_box_new(0,0)
        plt.show()
#         plt.scatter(time_values[mask_min_counts], time_counts[mask_min_counts])
#         plt.show()
        
    else:
        print("No suitable values found for", nrun)
        return
    if (len(x[mask])==0):
        print("Too few values in run", nrun)
        return
    
    
    max_tracks = 100 #stupid, but if neccessary just add a 0
    for N in range(1,max_tracks):
        tracks[N] = {}
        tracks[N]["x"] = []
        tracks[N]["y"] = []
        tracks[N]["z"] = []
        tracks[N]["pix"] = []
        tracks[N]["brightness"] = []
    for i in range(len(x[mask])):
        appended_to_track = False
        if i == 0:
            tracks[1]["x"].append(x[mask][i])
            tracks[1]["y"].append(y[mask][i])
            tracks[1]["z"].append(z[mask][i])
            tracks[1]["pix"].append(pix[mask][i])
            tracks[1]["brightness"].append(brightness[mask][i])
            appended_to_track  = True 
            continue
        for N in range(1,max_tracks):
            if N not in tracks:
                continue       
            if len(tracks[N]["pix"])==0:
                continue
            if tracks[N]["pix"][-1] in nn_pix[int(pix[mask][i])]:
                if z[mask][i] -2< tracks[N]["z"][-1]:
                    tracks[N]["x"].append(x[mask][i])
                    tracks[N]["y"].append(y[mask][i])
                    tracks[N]["z"].append(z[mask][i])
                    tracks[N]["pix"].append(pix[mask][i])
                    tracks[N]["brightness"].append(brightness[mask][i])
                    appended_to_track  = True 
                    break
            else:
                if len(tracks[N]["x"])<2:
                    tracks.pop(N, None)
          
                        
        if appended_to_track == False:
            for N in range(1,max_tracks):
                if N not in tracks:
                    continue
                if len(tracks[N]["x"])==0:
                    tracks[N]["x"].append(x[mask][i])
                    tracks[N]["y"].append(y[mask][i])
                    tracks[N]["z"].append(z[mask][i])
                    tracks[N]["pix"].append(pix[mask][i])
                    tracks[N]["brightness"].append(brightness[mask][i])
                    break
               
                        
    possible_meteorites = {}
    N_track = 1
    for N in range(1,max_tracks):
        if N not in tracks:
            continue
        if len(tracks[N]["x"]) == 0:
            tracks.pop(N, None)
            continue
        possible_meteorites[N_track] = {}
        #if len(tracks[N]["x"]) != 0 and len(np.unique(tracks[N]["z"]))<3:
        if len(np.unique(tracks[N]["pix"]))<3:
            print("Track", N, "from", np.unique(tracks[N]["z"]),
                  ",too few pixels", len(tracks[N]["pix"]))
            possible_meteorites[N_track]["x"] = tracks[N]["x"] 
            possible_meteorites[N_track]["y"] = tracks[N]["y"] 
            possible_meteorites[N_track]["z"] = tracks[N]["z"] 
            possible_meteorites[N_track]["pix"] = tracks[N]["pix"] 
            possible_meteorites[N_track]["brightness"] = tracks[N]["brightness"]
            tracks.pop(N, None)
        elif len(np.unique(tracks[N]["z"]))<3:
            print("Track", N, "from", np.unique(tracks[N]["z"]),
                  ",too few different times", set(tracks[N]["z"]))
            possible_meteorites[N_track]["x"] = tracks[N]["x"] 
            possible_meteorites[N_track]["y"] = tracks[N]["y"] 
            possible_meteorites[N_track]["z"] = tracks[N]["z"] 
            possible_meteorites[N_track]["pix"] = tracks[N]["pix"] 
            possible_meteorites[N_track]["brightness"] = tracks[N]["brightness"]
            tracks.pop(N, None)
            
        try:
            if possible_meteorites[N_track]["z"]:
                pass # add tracks with same time
        except:
            if len(np.unique(tracks[N]["pix"]))>2:
                plot_scatter(tracks[N]["x"], tracks[N]["y"], tracks[N]["z"])
                draw_box_new(0,0)
                plt.show()
        N_track+=1
        for N in reversed(range(1,max_tracks)):
            try:
                #print(N, np.unique(possible_meteorites[N]["z"]))
                if len(set(np.unique(possible_meteorites[N]["z"]) + np.min(possible_meteorites[N]["z"])-1) & 
                       set(np.unique(possible_meteorites[N-1]["z"]) + np.max(possible_meteorites[N-1]["z"])+1 ))>0:
                    print("yes")
                    possible_meteorites[N-1]["x"] = possible_meteorites[N-1]["x"]+possible_meteorites[N]["x"]
                    possible_meteorites[N-1]["y"] = possible_meteorites[N-1]["y"]+possible_meteorites[N]["y"]
                    possible_meteorites[N-1]["z"] = possible_meteorites[N-1]["z"]+possible_meteorites[N]["z"]
                    possible_meteorites[N-1]["pix"] = possible_meteorites[N-1]["pix"]+possible_meteorites[N]["pix"]
                    possible_meteorites[N-1]["brightness"] = possible_meteorites[N-1]["brightness"]+possible_meteorites[N]["brightness"]
                    possible_meteorites.pop(N, None)

                else:
                    if len(np.unique(possible_meteorites[N]["pix"]))<5:
                        possible_meteorites.pop(N, None)
            except:
                pass
        
            
    return tracks, possible_meteorites
    
   


# In[104]:


get_ipython().run_cell_magic('time', '', '# Remove all pixels that have more than max_counts occurences, since those are most of the time no satellite trails.\n# losses along a  trail have to be accepted. Max_counts is deliberately chosen generously\n# if unique pixels too many make new threshold\n#################Improved version#################\n#################Use this one#################\n\ntelId = 5\nmax_counts = 50\nfor nrun in new_nruns:\n    trails_high_avg[nrun] = {}\n    trails_high_avg[nrun][telId] = {}\n    try:\n        values, counts = np.unique(high_pixel[nrun]["pix"+str(telId)], return_counts=True)\n        #print(len(values))\n        mask_max_counts = counts<max_counts\n        #print(len(values[mask_max_counts]))\n        i=0\n        while len(values)>400:\n            mask = high_pixel[nrun][telId]>900+i*50\n            values, counts = np.unique(high_pixel[nrun]["pix"+str(telId)][mask], return_counts=True)\n            #print(len(values))\n            i+=1\n\n        mask_max_counts = counts<max_counts\n        print("new treshold at", 900+i*50)        \n        st = set(values[mask_max_counts])\n        result = [i for i, e in enumerate(high_pixel[nrun]["pix"+str(telId)]) if e in st]\n        \n        trails_high_avg[nrun][telId]["x"] = high_pixel[nrun]["x-pos"+str(telId)][result]\n        trails_high_avg[nrun][telId]["y"] = high_pixel[nrun]["y-pos"+str(telId)][result]\n        trails_high_avg[nrun][telId]["z"] = high_pixel[nrun]["tmp"+str(telId)][result]-high_pixel[nrun]["tmp"+str(telId)][0]\n        trails_high_avg[nrun][telId]["pix ID"] = high_pixel[nrun]["pix"+str(telId)][result]\n        trails_high_avg[nrun][telId]["brightness"] = high_pixel[nrun][telId][result]\n        plot_scatter(trails_high_avg[nrun][telId]["x"],\n                     trails_high_avg[nrun][telId]["y"],\n                     trails_high_avg[nrun][telId]["z"])\n        plt.show()\n\n    except:\n        print("No CT 5 for run ", nrun )\n        pass')


# In[112]:


get_ipython().run_cell_magic('time', '', 'telId = 5\ntrackruns = {}\npossible_meteorites = {}\nfor nrun in nruns:\n    trackruns[nrun] = {}\n    possible_meteorites[nrun] = {}\n    try:\n        \n        trackruns[nrun], possible_meteorites[nrun] = get_neighbour_pixel_entries(np.array(trails_high_avg[nrun][telId]["x"]),\n                                np.array(trails_high_avg[nrun][telId]["y"]),\n                                np.array(trails_high_avg[nrun][telId]["z"]),\n                                np.array(trails_high_avg[nrun][telId]["pix ID"]),\n                                np.array(trails_high_avg[nrun][telId]["brightness"]) )\n        \n        \n        \n        \n    except:\n        print(nrun, "ERROR")\n        pass\n    ')


# In[406]:


get_ipython().run_cell_magic('time', '', 'nrun = 174949\ntelId=5\n\nfig1, ax1 = plt.subplots(figsize=(8,6))\ndraw_box_new(0,0)\nfig2, ax2 = plt.subplots(figsize=(8,6))\ndraw_box_new(0,0)\nthresh = 2650\n\nprint(len(high_pixel[nrun]["x-pos"+str(telId)]))\nmask = high_pixel[nrun][telId]<6100\nb =  np.where(np.array(high_pixel[nrun][telId])>thresh, thresh, np.array(high_pixel[nrun][telId]))\ncax2 = ax2.scatter(high_pixel[nrun]["x-pos"+str(telId)], \n                 high_pixel[nrun]["y-pos"+str(telId)], \n                 c = b, cmap = "jet")\ncax1 = ax1.scatter(high_pixel[nrun]["x-pos"+str(telId)], \n                 high_pixel[nrun]["y-pos"+str(telId)], \n                 c = high_pixel[nrun]["tmp"+str(telId)] -high_pixel[nrun]["tmp"+str(telId)][0], cmap = "jet")\n\n\nax1.set_title("Timing in s")\nax1.set_xlabel("Camera x-pos")\nax1.set_ylabel("Camera y-pos")\ncbar1 = fig1.colorbar(cax1)\nax2.set_title("Brightness")\nax2.set_xlabel("Camera x-pos")\nax2.set_ylabel("Camera y-pos")\ncbar2 = fig2.colorbar(cax2)\nfig1.savefig("D:\\\\Masterarbeit ECAP\\\\Run-174949-uncleaned-timing.jpg")\nfig2.savefig("D:\\\\Masterarbeit ECAP\\\\Run-174949-uncleaned-brightness.jpg")\nplt.show()')


# In[464]:


get_ipython().run_cell_magic('time', '', '\nnrun = 174949\ntelId = 5\n# x = np.array(trails_high_avg[nrun][telId]["x"])\n# y = np.array(trails_high_avg[nrun][telId]["y"])\n# z = np.array(trails_high_avg[nrun][telId]["z"])\n# pix = np.array(trails_high_avg[nrun][telId]["pix ID"])\n# brightness = np.array(trails_high_avg[nrun][telId]["brightness"])\n\nx = np.array(high_pixel[nrun]["x-pos"+str(telId)])\ny = np.array(high_pixel[nrun]["y-pos"+str(telId)])\nz = np.array(high_pixel[nrun]["tmp"+str(telId)]- high_pixel[nrun]["tmp"+str(telId)][0])\npix = np.array(high_pixel[nrun]["pix"+str(telId)])\nbrightness = np.array(high_pixel[nrun][telId])\nPix = [ (pix[i],i) for i in np.arange(len(pix))]\nPix.sort()\n#Probiere nochmal mit List comprehension?\nsorted_pix,permutation = zip(*Pix)\nsorted_pix = np.array(sorted_pix)\npermutation = np.array(permutation)\n#print(len(set(high_pixel[nrun]["pix"+str(telId)])))\n# for pix in set(high_pixel[nrun]["pix"+str(telId)]):\n#     print(int(pix), end = ", ")\nmean_pix_pix = []\nmean_pix_x = []\nmean_pix_y = []\nmean_pix_z = []\nmean_pix_brightness = []\ncounter = 0\nprint(len(set(pix)))\nfor i in range(len(x)):\n    if i == 0:\n        #mean_x = [x[permutation][i]]\n        #mean_y = [y[permutation][i]]\n        mean_z = [z[permutation][i]]\n        mean_pix = [pix[permutation][i]]\n        mean_brightness = [brightness[permutation][i]]\n        continue\n    if pix[permutation][i] == pix[permutation][i-1]:\n        \n        #mean_x.append(x[permutation][i])\n        #mean_y.append(y[permutation][i])\n        mean_z.append(z[permutation][i])\n        mean_pix.append(pix[permutation][i])\n        mean_brightness.append(brightness[permutation][i])\n    else:\n        counter +=1\n        print(counter)\n        #mean_pix_x.append(np.average(mean_x))\n        #mean_pix_y.append(np.average(mean_y))\n        mean_pix_z.append(np.average(mean_z))\n        mean_pix_pix.append(np.average(mean_pix))\n        mean_pix_brightness.append(np.average(mean_brightness))\n        \n        #mean_x = [x[permutation][i]]\n        #mean_y = [y[permutation][i]]\n        mean_z = [z[permutation][i]]\n        mean_pix = [pix[permutation][i]]\n        mean_brightness = [brightness[permutation][i]]\nfig, ax = plt.subplots(figsize =(8,6))    \ndraw_box_new(0,0)\n\ncax = ax.scatter(np.array(mean_pix_x),\n                 np.array(mean_pix_y),\n                 c = np.array(mean_pix_brightness), cmap = "jet" )\nfig.colorbar(cax)\nplt.show()')


# In[484]:


print(mean_pix_pix)
for i in range(len(mean_pix_pix)):
    mean_pix_pix[i] = int(mean_pix_pix[i])
thresh = 4100
fig, ax = plt.subplots(figsize =(8,6))
b = np.where(np.array(mean_pix_brightness) > thresh, thresh, np.array(mean_pix_brightness))
    
cax =  ax.scatter(np.array(geom_hess5_xc_from_root)[np.array(mean_pix_pix)], 
           np.array(geom_hess5_yc_from_root)[np.array(mean_pix_pix)], 
           c = b, cmap = "jet")
draw_box_new(0,0) 
ax.set_xlim(-1.25,1.25)
ax.set_ylim(-1.25,1.25)
ax.set_title("Average brightness in MHz")
ax.set_xlabel("Camera x-pos")
ax.set_ylabel("Camera y-pos")
cbar = fig.colorbar(cax)
fig.savefig("D:\\Masterarbeit ECAP\\First Plots\\run174949-average-brightness.jpg")
plt.show()


# In[196]:


get_ipython().run_cell_magic('time', '', 'nrun = 174949\ntrackruns[nrun] = {}\npossible_meteorites[nrun] = {}\ntrackruns[nrun], possible_meteorites[nrun] = get_neighbour_pixel_entries(np.array(trails_high_avg[nrun][telId]["x"]),\n                                np.array(trails_high_avg[nrun][telId]["y"]),\n                                np.array(trails_high_avg[nrun][telId]["z"]),\n                                np.array(trails_high_avg[nrun][telId]["pix ID"]),\n                                np.array(trails_high_avg[nrun][telId]["brightness"]) )\nfor N in list(possible_meteorites[nrun].keys()):\n    print(0)\n    print(possible_meteorites[nrun][list(possible_meteorites[nrun].keys())][0])\n    #plot_scatter()')


# In[337]:


nrun = 174949
print(possible_meteorites[nrun].keys())
for N in reversed(list(possible_meteorites[nrun].keys())):
    #print(possible_meteorites[nrun][np.array(list(possible_meteorites[nrun].keys()))[0]].keys())
    try:
        try:
            if len(set(possible_meteorites[nrun][N]["z"]) & set(possible_meteorites[nrun][N-1]["z"]))>0:
                possible_meteorites[nrun][N-1]["x"] =  possible_meteorites[nrun][N-1]["x"]+ possible_meteorites[nrun][N]["x"]
                possible_meteorites[nrun][N-1]["y"] =  possible_meteorites[nrun][N-1]["y"]+ possible_meteorites[nrun][N]["y"]
                possible_meteorites[nrun][N-1]["z"] =  possible_meteorites[nrun][N-1]["z"]+ possible_meteorites[nrun][N]["z"]
                possible_meteorites[nrun][N-1]["pix"] =  possible_meteorites[nrun][N-1]["pix"]+ possible_meteorites[nrun][N]["pix"]
                possible_meteorites[nrun][N-1]["brightness"] =  possible_meteorites[nrun][N-1]["brightness"] + possible_meteorites[nrun][N]["brightness"]
                possible_meteorites[nrun].pop(N, None)
                continue
        except:
            pass
        print(set(possible_meteorites[nrun][N]["z"]))
        popt, pcov = curve_fit(lin_fit,possible_meteorites[nrun][N]["x"],possible_meteorites[nrun][N]["y"])
        draw_box_new(popt[0],popt[1])
        plot_scatter(possible_meteorites[nrun][N]["x"],
                    possible_meteorites[nrun][N]["y"],
                    possible_meteorites[nrun][N]["z"])
        plt.plot(np.arange(-1.25,1.25, 0.01), lin_fit(np.arange(-1.25,1.25, 0.01), *popt))
        try:
            print("Saving...")
            plt.savefig("D:\\Masterarbeit ECAP\\Run-174949-meteorite.jpg")
        except:
            print("saving did not work, try closing the file")
        plt.show()
        plot_scatter(possible_meteorites[nrun][N]["x"],
                    possible_meteorites[nrun][N]["y"],
                    possible_meteorites[nrun][N]["brightness"])
        draw_box_new(popt[0],popt[1])
        plt.show()
    except:
        try:
            print("error", set(possible_meteorites[nrun][N]["z"]))
            pass
        except:
            possible_meteorites[nrun].pop(N, None)
            print("no time found")
            
print(trackruns[nrun].keys())
counter = 0
for N in (list(trackruns[nrun].keys())):
    if N not in trackruns[nrun].keys():
        continue
    if "x" not in trackruns[nrun][N].keys():
        trackruns[nrun].pop(N, None)
        continue
    if len(set(trackruns[nrun][N]["pix"]))>4:
        popt, pcov = curve_fit(lin_fit,trackruns[nrun][N]["x"],trackruns[nrun][N]["y"])
        print(popt)
        draw_box_new(popt[0], popt[1])
        plot_scatter(trackruns[nrun][N]["x"],
                    trackruns[nrun][N]["y"],
                    trackruns[nrun][N]["z"])

        plt.plot(np.arange(-1.25,1.25, 0.01), lin_fit(np.arange(-1.25,1.25, 0.01), *popt))
        counter+=1
        plt.savefig("D:\\Masterarbeit ECAP\\Run-174949-track "+str(counter)+".jpg")
        plt.show()
        mask = np.array(trackruns[nrun][N]["brightness"])<np.average(np.array(trackruns[nrun][N]["brightness"]))+0.1*np.std(np.array(trackruns[nrun][N]["brightness"]))
        mask = np.array(trackruns[nrun][N]["brightness"])>np.average(np.array(trackruns[nrun][N]["brightness"]))-0.1*np.std(np.array(trackruns[nrun][N]["brightness"]))
        mask = np.array(trackruns[nrun][N]["brightness"]) > 0
        plot_scatter(np.array(trackruns[nrun][N]["x"])[mask],
                    np.array(trackruns[nrun][N]["y"])[mask],
                    np.array(trackruns[nrun][N]["brightness"])[mask])
        draw_box_new(popt[0],popt[1])
        plt.show()

    else:
        trackruns[nrun].pop(N, None)
   


# In[502]:


# print(possible_meteorites[nrun].keys())
fig, ax = plt.subplots(figsize=(8,6))
draw_box_new(0,0)
fig2, ax2 = plt.subplots(figsize=(8,6))
draw_box_new(0,0)
fig3, ax3 = plt.subplots(figsize=(8,6))
draw_box_new(0,0)
thresh = 2500
for N in range(np.max(list(possible_meteorites[nrun].keys()))+1):
    if N not in possible_meteorites[nrun]:
        continue
    if len(set(possible_meteorites[nrun][N]["pix"]))< 6:
        possible_meteorites[nrun].pop(N, None)
        continue
    print(set(possible_meteorites[nrun][N]["pix"]))
    ax.scatter(possible_meteorites[nrun][N]["x"],
                possible_meteorites[nrun][N]["y"],
                c = possible_meteorites[nrun][N]["z"]-possible_meteorites[nrun][N]["z"][0],  cmap = "jet")
    pix = list(possible_meteorites[nrun][N]["pix"])
    z = list(possible_meteorites[nrun][N]["z"]-possible_meteorites[nrun][N]["z"][0])
    brightness = list(possible_meteorites[nrun][N]["brightness"])
    print(set(possible_meteorites[nrun][N]["z"]- possible_meteorites[nrun][N]["z"][0]))
    b = np.where(np.array(possible_meteorites[nrun][N]["brightness"]) > thresh, thresh, np.array(possible_meteorites[nrun][N]["brightness"]))
    ax2.scatter(possible_meteorites[nrun][N]["x"],
                possible_meteorites[nrun][N]["y"],
                c = b, cmap = "jet")
    
    popt, pcov = curve_fit(lin_fit, possible_meteorites[nrun][N]["x"],possible_meteorites[nrun][N]["y"])
    ax.plot(np.arange(-1.,1.01), lin_fit(np.arange(-1.,1.01), *popt), "--", c = "black")
    ax2.plot(np.arange(-1.,1.01), lin_fit(np.arange(-1.,1.01), *popt), "--", c = "black")
    ax3.plot(np.arange(-1.,1.01), lin_fit(np.arange(-1.,1.01), *popt), "--", c = "black")
print(trackruns[nrun].keys())
for N in trackruns[nrun].keys():
    cax = ax.scatter(trackruns[nrun][N]["x"],
                trackruns[nrun][N]["y"],
                c = trackruns[nrun][N]["z"]- trackruns[nrun][N]["z"][0], cmap = "jet")
    pix = pix + list(trackruns[nrun][N]["pix"])
    z_pre = np.array(trackruns[nrun][N]["z"])-trackruns[nrun][N]["z"][0]
    z = z + list(trackruns[nrun][N]["z"]-trackruns[nrun][N]["z"][0])
    brightness = brightness + list(trackruns[nrun][N]["brightness"])
    print(set(trackruns[nrun][N]["z"]- trackruns[nrun][N]["z"][0]))
    popt, pcov = curve_fit(lin_fit, trackruns[nrun][N]["x"],trackruns[nrun][N]["y"])
    
    b = np.where(np.array(trackruns[nrun][N]["brightness"]) > thresh, thresh, np.array(trackruns[nrun][N]["brightness"]))
    
    
    cax2 = ax2.scatter(trackruns[nrun][N]["x"],
                trackruns[nrun][N]["y"],
                c = b, cmap = "jet")
    popt, pcov = curve_fit(lin_fit, trackruns[nrun][N]["x"],trackruns[nrun][N]["y"])
    ax.plot(np.arange(-1.,1.01), lin_fit(np.arange(-1.,1.01), *popt),  "--", c = "black")
    ax2.plot(np.arange(-1.,1.01), lin_fit(np.arange(-1.,1.01), *popt),  "--", c = "black")
    ax3.plot(np.arange(-1.,1.01), lin_fit(np.arange(-1.,1.01), *popt),  "--", c = "black")
    
ax.set_xlim(-1.25,1.25)
ax.set_ylim(-1.25,1.25)
ax2.set_xlim(-1.25,1.25)
ax2.set_ylim(-1.25,1.25)
draw_box_new(popt[0], popt[1])
cbar = fig.colorbar(cax)
#cbar.ax.set_yticklabels(['first in track', 'last in track '])
ax.set_title("Timing during each track")
ax.set_xlabel("Camera x-pos")
ax.set_ylabel("Camera y-pos")
cbar2 = fig2.colorbar(cax2) #ticks = [1000,2000,3000,4000, 6000])
cbar2.set_label("MHz")
#cbar2.ax2.set_yticklabels(['first in track', 'last in track '])
ax2.set_title("Brightness")
ax2.set_xlabel("Camera x-pos")
ax2.set_ylabel("Camera y-pos")
fig.savefig("D:\\Masterarbeit ECAP\\Run-174949-all-tracks-timing.jpg")
fig2.savefig("D:\\Masterarbeit ECAP\\Run-174949-all-tracks-brightness.jpg")
for i in range(len(pix)):
    pix[i] = int(pix[i])

ax3.set_xlim(-1.25,1.25)
ax3.set_ylim(-1.25,1.25)
cax3 = ax3.scatter(np.array(geom_hess5_xc_from_root)[np.array(pix)], 
            np.array(geom_hess5_yc_from_root)[np.array(pix)],
            c = z, cmap ="rainbow")
cbar3 = fig3.colorbar(cax3)
ax3.set_title("Timing during each track")
ax3.set_xlabel("Camera x-pos")
ax3.set_ylabel("Camera y-pos")
cbar3.set_label("Time in track [s]")

fig3.savefig("D:\\Masterarbeit ECAP\\Run-174949-all-tracks-timing-better.jpg")
plt.show()


# In[192]:


get_ipython().run_cell_magic('time', '', 'telId = 5\nnrun = 174949\nplot_scatter(-high_pixel[nrun]["y-pos"+str(telId)], \n             -high_pixel[nrun]["x-pos"+str(telId)], \n             high_pixel[nrun]["tmp"+str(telId)])\nplt.show()')


# In[976]:


for nrun in nruns:
    print(nrun)
    try:
        for N in range(1,np.max(list(possible_meteorites[nrun].keys()))+1):
            try:
                if len(np.unique(possible_meteorites[nrun][N]["pix"]))>0:
                    print(N, np.unique(possible_meteorites[nrun][N]["z"]))
                    popt, pcov = curve_fit(lin_fit, possible_meteorites[nrun][N]["x"], possible_meteorites[nrun][N]["y"])
                    
                    intersections = draw_box_new(popt[0],popt[1])
                    print("values for m, t:", np.around(popt, 6))
                    print("Covariances",np.around(np.diag(pcov),6))
                    
                    plt.plot(np.arange(-1.25,1.25, 0.01), lin_fit(np.arange(-1.25,1.25, 0.01), *popt))
                    plot_scatter(possible_meteorites[nrun][N]["x"], 
                                 possible_meteorites[nrun][N]["y"], 
                                 possible_meteorites[nrun][N]["z"])
                    plt.show()
            except:
                print("ERROR")
                pass
    except:
        pass


# In[29]:


get_ipython().run_cell_magic('time', '', 'for nrun in nruns:\n    telId = 1\n    for telId in range(1,5):\n        fig= plt.figure(figsize=(6.5,5), )\n        try:\n            print("CT", telId, "of run", nrun, end =": ")\n            plt.scatter(geom_hess1_xc_from_root, geom_hess1_yc_from_root, c=avg_pix[nrun][telId], s=70, cmap=\'turbo\')\n            plt.colorbar()\n            plt.show()\n            print\n        except:\n            print("No Data")\n            plt.close()\n    fig= plt.figure(figsize=(10,8), )\n    try:\n        \n        telId = 5\n        print("CT", telId, "of run", nrun, end =": ")\n        plt.scatter(geom_hess5_xc_from_root, geom_hess5_yc_from_root, c=avg_pix[nrun][telId], s=30, cmap=\'rainbow\')\n        plt.colorbar()\n        plt.show()\n    except: \n        print("No Data")\n        plt.close()\n        ')


# In[36]:


possible_peak_times = {}
for nrun in nruns:
    print("Run number ",nrun)
    telId = 1
    #test = avg_tmp[nrun]["tmp"+str(telId)]
    try:
        print(avg_tmp[nrun]["tmp"+str(telId)][0])
    except:
        pass
    fig = plt.figure(figsize=(12,6), )
   
    possible_peak_times[nrun]={}
    
    for telId in range(1,6):
    #    print("1st entry: " ,avg_tmp[nrun]["tmp"+str(telId)][0])
        try:
            possible_peak_times[nrun][telId] = {}
            tmp0 = avg_tmp[nrun]["tmp"+str(telId)][0]
            for i in range(len(avg_tmp[nrun]["tmp"+str(telId)])):            
                avg_tmp[nrun]["tmp"+str(telId)][i] = avg_tmp[nrun]["tmp"+str(telId)][i]- tmp0

                if telId==5:
                    if avg_tmp[nrun][telId][i]>np.average(avg_tmp[nrun][telId][i-10:i+10])+1:
                        possible_peak_times[nrun][telId][round(avg_tmp[nrun]["tmp"+str(telId)][i],1)] =  avg_tmp[nrun][telId][i]
                        #print("CT", telId, ": ", round(avg_tmp[nrun]["tmp"+str(telId)][i],1), ", " , avg_tmp[nrun][telId][i])
                else:
                    if avg_tmp[nrun][telId][i]>np.average(avg_tmp[nrun][telId][i-3:i+3])+4:
                        possible_peak_times[nrun][telId][round(avg_tmp[nrun]["tmp"+str(telId)][i],1)] =  avg_tmp[nrun][telId][i]
                        #print("CT", telId, ": ", round(avg_tmp[nrun]["tmp"+str(telId)][i],1), ", " , avg_tmp[nrun][telId][i])

            print("number of unique times vs time diff: ", len(np.unique(avg_tmp[nrun]["tmp"+str(telId)])), "vs", 
                  round(np.max(avg_tmp[nrun]["tmp"+str(telId)])-np.min(avg_tmp[nrun]["tmp"+str(telId)]),1))
            plt.scatter(avg_tmp[nrun]["tmp"+str(telId)], avg_tmp[nrun][telId], label = "CT"+str(telId))
            if telId == 4:
                #plt.xlim(100,200)
                plt.xlabel("Time from \"GetTime()\" [s]")
                plt.ylabel("NSB [MHz]")
                plt.legend()
                plt.show()

                fig = plt.figure(figsize=(12,6), )

        except:
            print("No CT"+str(telId)+" Data")
    plt.legend()
    plt.xlabel("Time from \"GetTime()\" [s]")
    #plt.xlim(117, 121)
    #plt.ylim(450, 470)
    plt.ylabel("NSB [MHz]")
    plt.show()
print(possible_peak_times)


# In[123]:


get_ipython().run_cell_magic('time', '', 'telId = 5\nfor nrun in nruns:\n    plot_scatter(high_pixel[nrun]["x-pos"+str(telId)], \n                 high_pixel[nrun]["y-pos"+str(telId)], \n                 high_pixel[nrun]["tmp"+str(telId)]-high_pixel[nrun]["tmp"+str(telId)][0] )\n    plt.show()\n    print("test")\n\n    #first cleaning:\n    x_cleaned_1,y_cleaned_1,z_cleaned_1 = cleaning_cut(high_pixel[nrun]["x-pos"+str(telId)], \n                                                       high_pixel[nrun]["y-pos"+str(telId)], \n                                                       high_pixel[nrun]["tmp"+str(telId)]-high_pixel[nrun]["tmp"+str(telId)][0], 3, 3)\n\n    #second cleaning\n    result = cleaning_cut(x_cleaned_1, y_cleaned_1, z_cleaned_1, 4, 2)\n\n\n    tracks = sort_into_tracks(result[0],result[1], result[2], 6)\n    for N in range(1,6):\n        if len(tracks["x-"+str(N)]) !=0:\n            number_of_tracks +=1\n            print("Track", N, end =": ")\n            popt, pcov = curve_fit(lin_fit, tracks["x-"+str(N)], tracks["y-"+str(N)])\n            print("linear fit : y =",round(popt[0],4),"* x +",round(popt[1],4)) \n            dur = tracks["z-"+str(N)][-1]- tracks["z-"+str(N)][0]\n            print("Duration: ", int(dur), "s" )\n            x_y_len = np.sqrt((tracks["x-"+str(N)][-1]-tracks["x-"+str(N)][0])**2+\n                              (lin_fit(np.array(tracks["x-"+str(N)]), *popt)[-1] - \n                               lin_fit(np.array(tracks["x-"+str(N)]), *popt)[0])**2 \n                             )\n            print("Length:", round(x_y_len, 4), "m(?)")\n            print("Speed on Cam: ", round(x_y_len/dur,4), "m/s")\n            x_y = []\n            for i in range(len(tracks["x-"+str(N)])):\n                x_y.append((round(tracks["x-"+str(N)][i],4), \n                            round(tracks["y-"+str(N)][i],4)))\n            print("Unique pixels:",len(list(set(x_y))))\n            print("Pixels per unit of length:", round(len(list(set(x_y)))/x_y_len, 4) )\n            plt.plot(np.array(tracks["x-"+str(N)]), lin_fit(np.array(tracks["x-"+str(N)]), *popt))\n            plt.plot()\n            plot_scatter(tracks["x-"+str(N)],tracks["y-"+str(N)], tracks["z-"+str(N)])\n            plt.show()\n    print("number of tracks: ", number_of_tracks)')


# In[25]:


get_ipython().run_cell_magic('time', '', 'telId = 5\nfor nrun in nruns:\n    \n    brightness_hist0, brightness_counts = np.unique(np.around(high_pixel[nrun][telId],-1), return_counts=True)    \n    plt.plot(brightness_hist0, brightness_counts)\n    plt.title(nrun)\n    plt.show()')


# In[66]:


get_ipython().run_cell_magic('time', '', 'nrun = 158214\ntelId = 5\n# for i in range(0,13):\n#     mask = (high_pixel[nrun]["tmp"+str(telId)]-high_pixel[nrun]["tmp"+str(telId)][0])>i*30\n#     mask = np.logical_and(mask,(high_pixel[nrun]["tmp"+str(telId)]-high_pixel[nrun]["tmp"+str(telId)][0])<(i+1)*30)\n#     plt.scatter(high_pixel[nrun]["x-pos"+str(telId)][mask], high_pixel[nrun]["y-pos"+str(telId)][mask], \n#                 c = high_pixel[nrun]["tmp"+str(telId)][mask]-high_pixel[nrun]["tmp"+str(telId)][0], cmap = "jet")\n#     plt.colorbar()\n#     plt.show()\nprint(nruns)\navg_from_high_pixel = []\nnruns_avg_from_high_pixel = []\nfor nrun in nruns:\n    print(np.average(high_pixel[nrun][telId]))\n    avg_from_high_pixel.append(np.average(high_pixel[nrun][telId]))\n    nruns_avg_from_high_pixel.append(nrun)\nplt.plot(nruns_avg_from_high_pixel, avg_from_high_pixel )\nplt.show()')


# In[308]:


get_ipython().run_cell_magic('time', '', 'telId = 5\nnrun = 158358\nmask = np.logical_and(high_pixel[nrun][telId]<10000,high_pixel[nrun][telId]>2000 )\nvalues, counts = np.unique(high_pixel[nrun]["pix"+str(telId)][mask], return_counts=True)\nmask2 = counts<30\nplt.scatter(values[mask2], counts[mask2])\nplt.show()\nfor i in range(1,17):\n    mask_time = np.logical_and((high_pixel[nrun]["tmp"+str(telId)][mask]-high_pixel[nrun]["tmp"+str(telId)][0])<100*(i+1),\n                              (high_pixel[nrun]["tmp"+str(telId)][mask]-high_pixel[nrun]["tmp"+str(telId)][0])>100*i)\n    \n    plot_scatter(high_pixel[nrun]["x-pos"+str(telId)][mask][mask_time],\n                 high_pixel[nrun]["y-pos"+str(telId)][mask][mask_time],\n                 high_pixel[nrun]["tmp"+str(telId)][mask][mask_time]-high_pixel[nrun]["tmp"+str(telId)][0])\n\n    plt.show()')


# In[945]:


get_ipython().run_cell_magic('time', '', 'nrun = nruns[0]\ntelId = 5\ntime_in_run = high_pixel[nrun]["tmp"+str(telId)]-high_pixel[nrun]["tmp"+str(telId)][0]\ntest_x = []\ntest_y = []\ntest_z = []\ntest_z_1 = []\ntest_z_2 = []\nvalues, counts = np.unique(high_pixel[nrun]["pix"+str(telId)], return_counts=True)\nval_num = []\nfor i in range(len(values)):\n    if counts[i]<30:\n        #print(values[i], counts[i])\n        val_num.append(values[i])\nprint(len(val_num))\nfor i in range(len(high_pixel[nrun]["x-pos"+str(telId)])):\n    if high_pixel[nrun][telId][i]<3000:\n        test_x.append(high_pixel[nrun]["x-pos"+str(telId)][i])\n        test_y.append(high_pixel[nrun]["y-pos"+str(telId)][i])\n        test_z_1.append(high_pixel[nrun][telId][i])\n        test_z_2.append(time_in_run[i])\nmask = np.array(test_z_2) < 345.\nmask = np.logical_and(mask, np.array(test_z_2)> 320.)\nprint("Unique times:", np.unique(np.array(test_z_2)[mask]))\nprint("minimum",np.min(high_pixel[nrun][telId]))\nprint(len(high_pixel[nrun]["x-pos"+str(telId)]))\nprint(len(test_z_1))\nplt.scatter(np.array(test_x)[mask],np.array(test_y)[mask],c=np.array(test_z_1)[mask])\nplt.colorbar()\nplt.show()\nplot_scatter(np.array(test_x)[mask],np.array(test_y)[mask],np.array(test_z_2)[mask])\nplt.show()\nbrightness_hist0, brightness_counts = np.unique(np.around(high_pixel[nrun][telId],0), return_counts=True)\nplt.plot(brightness_hist0, brightness_counts)\nplt.xlim(800,5000)\nplt.show()\n# x_cleaned_1,y_cleaned_1,z_cleaned_1 = cleaning_cut(high_pixel[nrun]["x-pos"+str(telId)], \n#                                                        high_pixel[nrun]["y-pos"+str(telId)], \n#                                                        time_in_run, 3, 3)\n# plot_scatter(x_cleaned_1, y_cleaned_1, z_cleaned_1)\n# plt.show()')


# In[21]:


get_ipython().run_cell_magic('time', '', '#Sort by pixels\ntmp_from_high = {}\nfor nrun in nruns:\n    tmp_from_high[nrun] = {}\n    for telId in range(1,6):\n        tmp_from_high[nrun][telId] = {}\n        tmp_from_high[nrun][telId]["all pixels"] = {}\n        try:\n            print("Run", nrun,"CT"+str(telId)+":", high_pixel[nrun]["pix"+str(telId)])\n            counter = 0\n            for pix in high_pixel[nrun]["pix"+str(telId)]:\n                if "tmp"+str(pix) in tmp_from_high[nrun][telId]:\n                    pass\n                else:\n                    tmp_from_high[nrun][telId]["tmp"+str(pix)] = []\n                    tmp_from_high[nrun][telId]["nsb"+str(pix)] = []\n                    tmp_from_high[nrun][telId]["x-pos"+str(pix)] = []\n                    tmp_from_high[nrun][telId]["y-pos"+str(pix)] = []\n                    tmp_from_high[nrun][telId]["all pixels"][pix] = pix \n                tmp_from_high[nrun][telId]["tmp"+str(pix)].append(high_pixel[nrun]["tmp"+str(telId)][counter])\n                tmp_from_high[nrun][telId]["nsb"+str(pix)].append(high_pixel[nrun][telId][counter])\n                tmp_from_high[nrun][telId]["x-pos"+str(pix)].append(high_pixel[nrun]["x-pos"+str(telId)][counter])\n                tmp_from_high[nrun][telId]["y-pos"+str(pix)].append(high_pixel[nrun]["y-pos"+str(telId)][counter])\n                #print(pix,tmp_from_roll[nrun][telId]["tmp"+str(pix)], tmp_from_roll[nrun][telId]["nsb"+str(pix)])\n                \n                counter = counter+1\n        except:\n            print("Run", nrun,",CT", telId,"does not exist")\n\n\n\nfor nrun in nruns:\n    print("")\n    print("Run number", nrun)\n    for telId in range(1,6):\n        print("")\n        print("Run", nrun,", CT", telId, end =": ")\n        counter = 0\n        for pix in sorted(tmp_from_high[nrun][telId]["all pixels"]):\n            print(int(tmp_from_high[nrun][telId]["all pixels"][pix]), end = ", ")\n            x = tmp_from_high[nrun][telId]["tmp"+str(pix)]\n            y = tmp_from_high[nrun][telId]["nsb"+str(pix)]\n            plt.plot(x,y, label = int(pix), lw = 0.5)\n            counter+=1\n            if counter >8:                \n                plt.legend()\n                plt.show()\n                counter = 0\n                print("Run", nrun,", CT", telId, end=": ")\n        if counter >0:\n            plt.legend()\n            plt.show()\n            ')


# In[22]:


get_ipython().run_cell_magic('time', '', 'nrun = 158212\ntelId = 4\nnew_x = []\nnew_y = []\nnew_time = []\ntime = high_pixel[nrun]["tmp"+str(telId)]-high_pixel[nrun]["tmp"+str(telId)][0]\nlowlim = 1430\nuplim = 1445\nfor i in range(len(high_pixel[nrun]["y-pos"+str(telId)])):\n    if lowlim<time[i] <uplim:\n        new_x.append(high_pixel[nrun]["x-pos"+str(telId)][i])\n        new_y.append(high_pixel[nrun]["y-pos"+str(telId)][i])\n        new_time.append(time[i])\nplt.scatter(new_x, new_y, c = new_time, cmap = "jet")\nplt.colorbar(label = "time in s")\nplt.xlabel("x-coordinate")\nplt.ylabel("y-coordinate")\nplt.savefig(str(nrun)+" CT 5.pdf")\nplt.show()\npix_for_now = []\nfor pix in sorted(tmp_from_high[nrun][telId]["all pixels"]):\n    #print(int(tmp_from_high[nrun][telId]["all pixels"][pix]), end = ", ")\n    tmp = tmp_from_high[nrun][telId]["tmp"+str(pix)]-high_pixel[nrun]["tmp"+str(telId)][0]\n    nsb = tmp_from_high[nrun][telId]["nsb"+str(pix)]\n    print("pixel", pix, end=": ")\n    lower = np.where(lowlim<=tmp)\n    upper = np.where(uplim>=tmp)\n    try:\n        print(list(tmp[lower[0]])[0], end ="") \n    except:\n        pass\n    try:\n        print(list(tmp[upper[-1]])[-1])\n    except:\n        print("")\n    if np.std(nsb)>500:\n        print(pix, ", ",np.std(nsb))\n        plt.plot(tmp,nsb, label = int(pix), lw = 0.5)\n        pix_for_now.append(pix)\nplt.xlim(lowlim,uplim)\nplt.ylim(0,7000)\nplt.show()\n#print(pix_for_now)\nnewnew_x = []\nnewnew_y = []\nprint(len(pix_for_now))\nfor i in range(0,1):\n    print(i)\n#     pixel = pix_for_now[i]\n#     newnew_x.append(geom_hess5_xc_from_root[pixel])\n#     newnew_y.append(geom_hess5_yc_from_root[pixel])\n#     print(pixel)\nplt.plot(newnew_x, newnew_y)\nplt.show()')


# In[46]:




for nrun in nruns:
    print("")
    print("Run number", nrun)
    for telId in range(1,6):
        print("")
        print("CT", telId, end =": ")
        counter = 0
        for pix in sorted(tmp_from_roll[nrun][telId]["all pixels"]):
            print(int(tmp_from_roll[nrun][telId]["all pixels"][pix]), end = ", ")
            x = tmp_from_roll[nrun][telId]["tmp"+str(pix)]
            y = tmp_from_roll[nrun][telId]["nsb"+str(pix)]
            plt.scatter(x,y, label = int(pix))
            counter+=1
            if counter >8:                
                plt.legend()
                plt.show()
                counter = 0
                print("CT", telId, end=": ")
        if counter >0:
            plt.legend()
            plt.show()


# In[27]:


get_ipython().run_cell_magic('time', '', 'telId = 5\nnrun = nruns[1]\nprint(len(high_pixel[nrun]["x-pos"+str(telId)]))\nplot_scatter(high_pixel[nrun]["x-pos"+str(telId)],\n            high_pixel[nrun]["y-pos"+str(telId)],\n            high_pixel[nrun]["tmp"+str(telId)])\ndraw_box_new(0,0)\nplt.show()')


# In[122]:


get_ipython().run_cell_magic('time', '', '# Remove all pixels that have more than max_counts occurences, since those are most of the time no satellite trails.\n# losses along a satellite trail have to be accepted. Max_counts is deliberately chosen generously\n#################Old version#################\n#################Do Not Use this one#########\ntrails_high_avg = {}\ntelId = 5\nmax_counts = 50\n\nfor nrun in nruns:\n    trails_high_avg[nrun] = {}\n    for telId in range(5,6):\n        trails_high_avg[nrun][telId] = {}\n        try:\n            values, counts = np.unique(high_pixel[nrun]["pix"+str(telId)], return_counts=True)\n            print("Length of values", len(values))\n            trails_high_avg[nrun][telId]["unique pix ID"] = []\n            trails_high_avg[nrun][telId]["pix ID counts"] = []\n            new_val = []\n            new_counts = []\n            for i in range(len(values)):\n                if counts[i]<max_counts:#Number needs to be calculated for slowest satellite\n                    new_val.append(values[i])\n                    new_counts.append(counts[i])\n\n                    trails_high_avg[nrun][telId]["unique pix ID"].append(values[i])\n                    trails_high_avg[nrun][telId]["pix ID counts"].append(counts[i])\n\n            trails_high_avg[nrun][telId]["pix ID"] = []\n            trails_high_avg[nrun][telId]["x"] = []\n            trails_high_avg[nrun][telId]["y"] = []\n            trails_high_avg[nrun][telId]["z"] = []\n            trails_high_avg[nrun][telId]["brightness"] = []\n\n            print("Run", nrun, ", CT", telId,":")\n            for i in range(len(high_pixel[nrun]["pix"+str(telId)])):\n                if high_pixel[nrun]["pix"+str(telId)][i] in trails_high_avg[nrun][telId]["unique pix ID"]:\n\n                    if high_pixel[nrun]["tmp"+str(telId)][i]-high_pixel[nrun]["tmp"+str(telId)][0] <5:\n                        pass\n\n                    else:\n                        trails_high_avg[nrun][telId]["pix ID"].append(int(high_pixel[nrun]["pix"+str(telId)][i]))\n                        trails_high_avg[nrun][telId]["x"].append(high_pixel[nrun]["x-pos"+str(telId)][i])\n                        trails_high_avg[nrun][telId]["y"].append(high_pixel[nrun]["y-pos"+str(telId)][i])\n                        trails_high_avg[nrun][telId]["z"].append(high_pixel[nrun]["tmp"+str(telId)][i]-high_pixel[nrun]["tmp"+str(telId)][0])\n                        trails_high_avg[nrun][telId]["brightness"].append(high_pixel[nrun][telId][i])\n\n            plt.scatter(trails_high_avg[nrun][telId]["x"],\n                        trails_high_avg[nrun][telId]["y"],\n                        c=trails_high_avg[nrun][telId]["z"], \n                        cmap = "jet")\n            plt.colorbar()\n            plt.xlim(-1.2,1.2)\n            plt.ylim(-1.2,1.2)\n            plt.show()\n        except:\n            pass')


# Plan: GetNeighboursList.C schreiben, um PixId von neighbours and nextneigbours zu erhalten. 
# => Beschleunigung von cleaning

# width of trail (done),
# 
# angular distance, duration (done to second timescale), speed (m/s, pixels/s(done) )
# 
# direction of run (zenith) (still need code snippet to insert in C++ )
# 
# CT1-4 cuts, finding trails

# In[110]:


get_ipython().run_cell_magic('time', '', '#Version for only CT5\n#Remove all pixels that have too few neighbouring pixels by cleaning twice\ntelId = 5\nx_values = np.arange(-1.2,1.2, 0.01)\nfor nrun in nruns:\n    try:\n        print("CT", telId, ", Run", nrun, ":")\n        plot_scatter(trails_high_avg[nrun][telId]["x"], trails_high_avg[nrun][telId]["y"], \n                    trails_high_avg[nrun][telId]["z"])    \n        plt.show()\n        #first  cleaning:\n        x_cleaned_1,y_cleaned_1,z_cleaned_1, brightness_cleaned_1 = cleaning_cut(trails_high_avg[nrun][telId]["x"],\n                                                                                 trails_high_avg[nrun][telId]["y"],\n                                                                                 trails_high_avg[nrun][telId]["z"],\n                                                                                 trails_high_avg[nrun][telId]["brightness"],4, 3)\n\n        #second cleaning\n        result = cleaning_cut(x_cleaned_1, y_cleaned_1, z_cleaned_1, brightness_cleaned_1, 4, 3)\n        trails_high_avg[nrun][telId]["x_cleaned"] = result[0]\n        trails_high_avg[nrun][telId]["y_cleaned"] = result[1]\n        trails_high_avg[nrun][telId]["z_cleaned"] = result[2]\n        trails_high_avg[nrun][telId]["brightness_cleaned"] = result[3]\n\n        tracks, meterorites = sort_into_tracks(trails_high_avg[nrun][telId]["x_cleaned"],\n                                               trails_high_avg[nrun][telId]["y_cleaned"],\n                                               trails_high_avg[nrun][telId]["z_cleaned"], \n                                               trails_high_avg[nrun][telId]["brightness_cleaned"], 6)\n        trails_high_avg[nrun][telId]["tracks"] = tracks\n        print("Trying to separate:")\n        number_of_tracks = 0\n        for N in range(1,6):\n            if len(tracks["x-"+str(N)]) !=0:\n                number_of_tracks +=1\n                trails_high_avg[nrun][telId]["tracks"]["x-"+str(number_of_tracks)] = trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)]\n                trails_high_avg[nrun][telId]["tracks"]["y-"+str(number_of_tracks)] = trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)]\n                trails_high_avg[nrun][telId]["tracks"]["z-"+str(number_of_tracks)] = trails_high_avg[nrun][telId]["tracks"]["z-"+str(N)]\n                trails_high_avg[nrun][telId]["tracks"]["brightness-"+str(number_of_tracks)] = np.around(trails_high_avg[nrun][telId]["tracks"]["brightness-"+str(N)],1)\n                if N != number_of_tracks:\n                    trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)] = []\n                    trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)] = []\n                    trails_high_avg[nrun][telId]["tracks"]["z-"+str(N)] = []\n                    trails_high_avg[nrun][telId]["tracks"]["brightness-"+str(N)] = []\n                print("Track", number_of_tracks, end =": ")\n                popt, pcov = curve_fit(lin_fit, tracks["x-"+str(number_of_tracks)], tracks["y-"+str(number_of_tracks)])\n                print("linear fit : y =",round(popt[0],4),"* x +",round(popt[1],4)) \n                dur = tracks["z-"+str(number_of_tracks)][-1]- tracks["z-"+str(number_of_tracks)][0]\n                print("Duration: ", int(dur), "s" )\n                x_y_len = np.sqrt((tracks["x-"+str(number_of_tracks)][-1]-tracks["x-"+str(number_of_tracks)][0])**2+\n                                  (lin_fit(np.array(tracks["x-"+str(number_of_tracks)]), *popt)[-1] - \n                                   lin_fit(np.array(tracks["x-"+str(number_of_tracks)]), *popt)[0])**2 \n                                 )\n                print("Length:", round(x_y_len, 4), "m(?)")\n                print("Average brightness:", round(np.average(tracks["brightness-"+str(number_of_tracks)]),1) , "+-", \n                      round(np.std(tracks["brightness-"+str(number_of_tracks)]), 1))\n                print("Speed on Cam: ", round(x_y_len/dur,4), "m/s")\n                x_y = []\n                for i in range(len(tracks["x-"+str(number_of_tracks)])):\n                    x_y.append((round(tracks["x-"+str(number_of_tracks)][i],4), \n                                round(tracks["y-"+str(number_of_tracks)][i],4)))\n                print("Unique pixels:",len(list(set(x_y))))\n                print("Pixels per unit of length:", round(len(list(set(x_y)))/x_y_len, 4) )             \n                get_track_width(trails_high_avg[nrun][telId]["tracks"]["x-"+str(number_of_tracks)],\n                               trails_high_avg[nrun][telId]["tracks"]["y-"+str(number_of_tracks)])\n                intersection = draw_box_new(popt[0], popt[1])\n                trails_high_avg[nrun][telId]["tracks"]["IS track "+str(number_of_tracks)] = []\n                trails_high_avg[nrun][telId]["tracks"]["IS track "+str(number_of_tracks)].append(intersection)\n                plt.plot(np.array(tracks["x-"+str(number_of_tracks)]), lin_fit(np.array(tracks["x-"+str(number_of_tracks)]), *popt))\n                plot_scatter(tracks["x-"+str(number_of_tracks)],tracks["y-"+str(number_of_tracks)], tracks["z-"+str(number_of_tracks)])\n    #           \'  \n    #             print_important_params(tracks["x-"+str(N)],\n    #                                    tracks["y-"+str(N)],\n    #                                    tracks["z-"+str(N)])\'\n\n                plt.savefig("run_"+str(nrun)+"_CT_"+str(telId)+"_track_"+str(number_of_tracks)+".pdf")\n                plt.show()\n                brightness_hist0, brightness_counts = np.unique(np.around(tracks["brightness-"+str(number_of_tracks)],-3), return_counts=True)\n                plt.plot(brightness_hist0, brightness_counts)\n                plt.show()\n\n\n        print("Number of satellite tracks found in run", nrun, "CT", telId, ": ", number_of_tracks)\n        print("Number of tracks identified as metiorites:", len(meterorites.keys()))\n    except:\n        print("No CT 5 data available")')


# In[113]:


get_ipython().run_cell_magic('time', '', '#version without cleaning first\ntelId = 5 \nx_values = np.arange(-1.2,1.2, 0.01)\nfor nrun in nruns:\n    try:\n        print("CT", telId, ", Run", nrun, ":")\n        plot_scatter(trails_high_avg[nrun][telId]["x"], trails_high_avg[nrun][telId]["y"], \n                    trails_high_avg[nrun][telId]["z"])    \n        plt.show()\n        tracks, meterorites = sort_into_tracks(trails_high_avg[nrun][telId]["x"],\n                                               trails_high_avg[nrun][telId]["y"],\n                                               trails_high_avg[nrun][telId]["z"], \n                                               trails_high_avg[nrun][telId]["brightness"], 6)\n        trails_high_avg[nrun][telId]["tracks"] = tracks\n        print("Trying to separate:")\n        number_of_tracks = 0\n        for N in range(1,9):\n            if len(tracks["x-"+str(N)]) !=0:\n                number_of_tracks +=1\n                trails_high_avg[nrun][telId]["tracks"]["x-"+str(number_of_tracks)] = trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)]\n                trails_high_avg[nrun][telId]["tracks"]["y-"+str(number_of_tracks)] = trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)]\n                trails_high_avg[nrun][telId]["tracks"]["z-"+str(number_of_tracks)] = trails_high_avg[nrun][telId]["tracks"]["z-"+str(N)]\n                trails_high_avg[nrun][telId]["tracks"]["brightness-"+str(number_of_tracks)] = np.around(trails_high_avg[nrun][telId]["tracks"]["brightness-"+str(N)],1)\n                if N != number_of_tracks:\n                    trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)] = []\n                    trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)] = []\n                    trails_high_avg[nrun][telId]["tracks"]["z-"+str(N)] = []\n                    trails_high_avg[nrun][telId]["tracks"]["brightness-"+str(N)] = []\n                print("Track", number_of_tracks, end =": ")\n                popt, pcov = curve_fit(lin_fit, tracks["x-"+str(number_of_tracks)], tracks["y-"+str(number_of_tracks)])\n                print("linear fit : y =",round(popt[0],4),"* x +",round(popt[1],4)) \n                dur = tracks["z-"+str(number_of_tracks)][-1]- tracks["z-"+str(number_of_tracks)][0]\n                print("Duration: ", int(dur), "s" )\n                x_y_len = np.sqrt((tracks["x-"+str(number_of_tracks)][-1]-tracks["x-"+str(number_of_tracks)][0])**2+\n                                  (lin_fit(np.array(tracks["x-"+str(number_of_tracks)]), *popt)[-1] - \n                                   lin_fit(np.array(tracks["x-"+str(number_of_tracks)]), *popt)[0])**2 \n                                 )\n                print("Length:", round(x_y_len, 4), "m(?)")\n                print("Average brightness:", round(np.average(tracks["brightness-"+str(number_of_tracks)]),1) , "+-", \n                      round(np.std(tracks["brightness-"+str(number_of_tracks)]), 1))\n                print("Speed on Cam: ", round(x_y_len/dur,4), "m/s")\n                x_y = []\n                for i in range(len(tracks["x-"+str(number_of_tracks)])):\n                    x_y.append((round(tracks["x-"+str(number_of_tracks)][i],4), \n                                round(tracks["y-"+str(number_of_tracks)][i],4)))\n                print("Unique pixels:",len(list(set(x_y))))\n                print("Pixels per unit of length:", round(len(list(set(x_y)))/x_y_len, 4) )             \n                get_track_width(trails_high_avg[nrun][telId]["tracks"]["x-"+str(number_of_tracks)],\n                               trails_high_avg[nrun][telId]["tracks"]["y-"+str(number_of_tracks)])\n                intersection = draw_box_new(popt[0], popt[1])\n                trails_high_avg[nrun][telId]["tracks"]["IS track "+str(number_of_tracks)] = []\n                trails_high_avg[nrun][telId]["tracks"]["IS track "+str(number_of_tracks)].append(intersection)\n                plt.plot(np.array(tracks["x-"+str(number_of_tracks)]), lin_fit(np.array(tracks["x-"+str(number_of_tracks)]), *popt))\n                plot_scatter(tracks["x-"+str(number_of_tracks)],tracks["y-"+str(number_of_tracks)], tracks["z-"+str(number_of_tracks)])\n                plt.show()\n        print("Number of satellite tracks found in run", nrun, "CT", telId, ": ", number_of_tracks)\n        print("Number of tracks identified as metiorites:", len(meterorites.keys()))\n    except:\n        print("No CT 5 data available")')


# In[84]:


get_ipython().run_cell_magic('time', '', '#Version for all CTs\n#Remove all pixels that have too few neighbouring pixels by cleaning twice\ntelId = 5\nx_values = np.arange(-1.2,1.2, 0.01)\nfor nrun in nruns:\n    for telId in range(5,6):\n        try:\n            print("CT", telId, ", Run", nrun, ":")\n            plot_scatter(trails_high_avg[nrun][telId]["x"], trails_high_avg[nrun][telId]["y"], \n                        trails_high_avg[nrun][telId]["z"])    \n            plt.show()\n            #first  cleaning:\n            x_cleaned_1,y_cleaned_1,z_cleaned_1, brightness_cleaned_1 = cleaning_cut(trails_high_avg[nrun][telId]["x"],\n                                                                                     trails_high_avg[nrun][telId]["y"],\n                                                                                     trails_high_avg[nrun][telId]["z"],\n                                                                                     trails_high_avg[nrun][telId]["brightness"],3, 3)\n\n            #second cleaning\n            result = cleaning_cut(x_cleaned_1, y_cleaned_1, z_cleaned_1, brightness_cleaned_1, 4, 2)\n            trails_high_avg[nrun][telId]["x_cleaned"] = result[0]\n            trails_high_avg[nrun][telId]["y_cleaned"] = result[1]\n            trails_high_avg[nrun][telId]["z_cleaned"] = result[2]\n            trails_high_avg[nrun][telId]["brightness_cleaned"] = result[3]\n\n            tracks, meterorites = sort_into_tracks(trails_high_avg[nrun][telId]["x_cleaned"],\n                                                   trails_high_avg[nrun][telId]["y_cleaned"],\n                                                   trails_high_avg[nrun][telId]["z_cleaned"], \n                                                   trails_high_avg[nrun][telId]["brightness_cleaned"], 6)\n            trails_high_avg[nrun][telId]["tracks"] = tracks\n            print("Trying to separate:")\n            number_of_tracks = 0\n            for N in range(1,6):\n                if len(tracks["x-"+str(N)]) !=0:\n                    number_of_tracks +=1\n                    print("Track", N, end =": ")\n                    popt, pcov = curve_fit(lin_fit, tracks["x-"+str(N)], tracks["y-"+str(N)])\n                    print("linear fit : y =",round(popt[0],4),"* x +",round(popt[1],4)) \n                    dur = tracks["z-"+str(N)][-1]- tracks["z-"+str(N)][0]\n                    print("Duration: ", int(dur), "s" )\n                    x_y_len = np.sqrt((tracks["x-"+str(N)][-1]-tracks["x-"+str(N)][0])**2+\n                                      (lin_fit(np.array(tracks["x-"+str(N)]), *popt)[-1] - \n                                       lin_fit(np.array(tracks["x-"+str(N)]), *popt)[0])**2 \n                                     )\n                    print("Length:", round(x_y_len, 4), "m(?)")\n                    print("Average brightness:", round(np.average(tracks["brightness-"+str(N)]),1) , "+-", \n                          round(np.std(tracks["brightness-"+str(N)]), 1))\n                    print("Speed on Cam: ", round(x_y_len/dur,4), "m/s")\n                    x_y = []\n                    for i in range(len(tracks["x-"+str(N)])):\n                        x_y.append((round(tracks["x-"+str(N)][i],4), \n                                    round(tracks["y-"+str(N)][i],4)))\n                    print("Unique pixels:",len(list(set(x_y))))\n                    print("Pixels per unit of length:", round(len(list(set(x_y)))/x_y_len, 4) )             \n                    get_track_width(trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)],\n                                   trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)])\n                    intersection = draw_box_new(popt[0], popt[1])\n                    trails_high_avg[nrun][telId]["tracks"]["IS track "+str(number_of_tracks)] = []\n                    trails_high_avg[nrun][telId]["tracks"]["IS track "+str(number_of_tracks)].append(intersection)\n                    plt.plot(np.array(tracks["x-"+str(N)]), lin_fit(np.array(tracks["x-"+str(N)]), *popt))\n                    plot_scatter(tracks["x-"+str(N)],tracks["y-"+str(N)], tracks["z-"+str(N)])\n        #           \'  \n        #             print_important_params(tracks["x-"+str(N)],\n        #                                    tracks["y-"+str(N)],\n        #                                    tracks["z-"+str(N)])\'\n\n                    plt.savefig("run_"+str(nrun)+"_CT_"+str(telId)+"_track_"+str(number_of_tracks)+".pdf")\n                    plt.show()\n\n            print("Number of satellite tracks found in run", nrun, "CT", telId, ": ", number_of_tracks)\n            print("Number of tracks identified as metiorites:", len(meterorites.keys()))\n        except:\n            pass')


# In[92]:


nrun = 158385
telId = 5
print(high_pixel.keys())
mask = np.array(trails_high_avg[nrun][telId]["z"]) < 800.
mask = np.logical_and(mask, np.array(trails_high_avg[nrun][telId]["z"])>790.)
print(mask)
plot_scatter(np.array(trails_high_avg[nrun][telId]["x"])[mask],
             np.array(trails_high_avg[nrun][telId]["y"])[mask], 
             np.array(trails_high_avg[nrun][telId]["z"])[mask])
print(np.array(trails_high_avg[nrun][telId]["z"])[mask])
plt.show()


# In[147]:


get_ipython().run_cell_magic('time', '', 'nrun = 158229\ntelId = 5\n#print(high_pixel[nrun].keys())\nprint(len(high_pixel[nrun][telId]))\nmask = high_pixel[nrun][telId] >200. \nprint(len(np.where(mask==True)[0]))\nmask = np.logical_and(mask, high_pixel[nrun][telId] <5000)\nprint(len(np.where(mask==True)[0]))\n#mask = np.logical_and(mask, high_pixel[nrun]["tmp"+str(telId)]-high_pixel[nrun]["tmp"+str(telId)][0]>350)\nprint(len(np.where(mask==True)[0]))\nplot_scatter(high_pixel[nrun]["x-pos"+str(telId)][mask], high_pixel[nrun]["y-pos"+str(telId)][mask],\n             high_pixel[nrun]["tmp"+str(telId)][mask]-high_pixel[nrun]["tmp"+str(telId)][0])\nif telId != 5:\n    plt.xlim(-0.75,0.75)\n    plt.ylim(-0.75,0.75)\nplt.show()')


# get charateristic parameters for each run which are of interest:
# nrun, telid, direction (zenith, azimuth or sky coordinates RA, DEC), trail direction (maybe in appropriate coordinate system, 
# duration, start-/end-time in UTC, mean brightness of trail, 
# 
# applying mask of CT5 on CT1-4 
# 
# mon 11:00

# In[59]:


nrun = 158229
telId = 5
for N in range(1,6):
    if len(trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)])!=0:
        popt, pcov = curve_fit(lin_fit, trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)], 
                              trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)])
        print("linear fit : y =",round(popt[0],4),"* x +",round(popt[1],4)) 
        plot_scatter(trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)],
                     trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)],
                     trails_high_avg[nrun][telId]["tracks"]["z-"+str(N)])
        
        plt.plot(np.array(trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)]), lin_fit(np.array(trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)]), *popt))
        
        d_from_lin_fit = []
        for i in range(len(trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)])):
            x_p = trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)][i]
            y_p = trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)][i]
            t_2 = x_p/popt[0]+y_p
            x_s = popt[0]/(popt[0]**2+1)*(t_2 -popt[1])
            y_s = popt[0]*x_s+popt[1]
            d_from_lin_fit.append(np.sqrt(((x_p-x_s)/min_x_dif)**2+((y_p-y_s)/min_y_dif)**2))
        print(np.average(d_from_lin_fit))
        get_track_width(trails_high_avg[nrun][telId]["tracks"]["x-"+str(N)], 
                       trails_high_avg[nrun][telId]["tracks"]["y-"+str(N)])
        draw_box(popt[0], popt[1])
        plt.show()


# In[ ]:


#linear fits, direction of travel and speed
#popt, pcov = curve_fit(lin_fit, track_x_1, track_y_1)
#vector:
# v =[track_x_1[-1]-track_x_1[0],track_y_1[-1]-track_y_1[-0]]


# In[294]:


get_ipython().run_cell_magic('time', '', 'nrun = 158229\nprint("CT", telId, ", Run", nrun, ":")\nplt.scatter(trails_high_avg[nrun][telId]["x_cleaned"], trails_high_avg[nrun][telId]["y_cleaned"],\n           c= trails_high_avg[nrun][telId]["z_cleaned"], cmap = "jet")\nplt.xlim(-1.2,1.2)\nplt.ylim(-1.2,1.2)\nplt.colorbar()\nplt.show()\ntrack_x_1 = []\ntrack_x_2 = []\ntrack_y_1 = []\ntrack_y_2 = []\ntrack_z_1 = []\ntrack_z_2 = []\nfor i in range(len(trails_high_avg[nrun][telId]["x_cleaned"])):      \n    x = trails_high_avg[nrun][telId]["x_cleaned"][i]\n    y = trails_high_avg[nrun][telId]["y_cleaned"][i]\n    z = trails_high_avg[nrun][telId]["z_cleaned"][i]\n    #neighbours and next neighbours\n    nn_x_xcut = [] \n    nn_y_xcut = []\n    nn_z_xcut = []\n    nn_x = [] \n    nn_y = []\n    nn_z = []\n    for j in range(len(trails_high_avg[nrun][telId]["x_cleaned"])):\n        if (x-2*min_x_dif<\n            trails_high_avg[nrun][telId]["x_cleaned"][j] < \n            x + 2*min_x_dif):\n            nn_x_xcut.append(trails_high_avg[nrun][telId]["x_cleaned"][j])\n            nn_y_xcut.append(trails_high_avg[nrun][telId]["y_cleaned"][j])\n            nn_z_xcut.append(trails_high_avg[nrun][telId]["z_cleaned"][j])\n    for j in range(len(nn_x_xcut)):\n        if (y-2*min_y_dif<\n            nn_y_xcut[j] < \n            y + 2*min_y_dif):\n            nn_x.append(nn_x_xcut[j])\n            nn_y.append(nn_y_xcut[j])\n            nn_z.append(nn_z_xcut[j])\n    \n    if np.max(nn_x)-np.min(nn_x)> 4*min_x_dif:\n        print(np.max(nn_x)-np.min(nn_x), "max allowed: ", 4*min_x_dif)\n    if np.max(nn_y)-np.min(nn_y)> 4*min_y_dif:\n        print(np.max(nn_y)-np.min(nn_y), "max allowed: ", 4*min_y_dif)\n#     if i == 563:\n#         plt.scatter(nn_x,nn_y, c = nn_z)\n#         plt.xlim(-1.2,1.2)\n#         plt.ylim(-1.2,1.2)\n#         plt.colorbar()\n#         plt.show()\n    # Insert into tracks\n    if i == 0:\n        track_x_1.append(x)\n        track_y_1.append(y)\n        track_z_1.append(z)\n    else:\n        if track_x_1[-1] in nn_x and track_y_1[-1] in nn_y:\n            track_x_1.append(x)\n            track_y_1.append(y)\n            track_z_1.append(z)\n        elif len(track_x_2) == 0:\n            track_x_2.append(x)\n            track_y_2.append(y)\n            track_z_2.append(z)\n        elif track_x_2[-1] in nn_x and track_y_2[-1] in nn_y:\n            track_x_2.append(x)\n            track_y_2.append(y)\n            track_z_2.append(z)\n            \n            \n            \n        \n        \nif len(track_x_1)!=0:      \n    print("Trying to separate:")        \n    plt.scatter(track_x_1, track_y_1, c = track_z_1, cmap = "jet")\n    plt.xlim(-1.2,1.2)\n    plt.ylim(-1.2,1.2)\n    plt.colorbar()\n    plt.arrow(track_x_1[0],track_y_1[0], track_x_1[-1]-track_x_1[0],track_y_1[-1]-track_y_1[0])\n    plt.show()\n    print("Time difference:", track_z_1[-1]-track_z_1[0], "s")\n    print("Maximum distance on cam:", round(np.sqrt((track_x_1[-1]-track_x_1[0])**2 +(track_y_1[-1]-track_y_1[0])**2 ),3))\n    if len(np.unique(track_z_1)) <3:\n        print("Probably metiorite")\n        print("")\nprint(len(track_x_1))\nplt.plot(track_x_1[540:580], track_y_1[540:580])\nplt.xlim(-1.2,1.2)\nplt.ylim(-1.2,1.2)\nplt.show()\nprint(len(track_x_2))')


# In[76]:


a = [1,1,1,4,5]
b = [2,1,2,1,4]
a_b = []
for i in range(len(a)):
    a_b.append((a[i], b[i]))
print(a_b)
print(len(list(set(a_b))))


# In[502]:


x_fov = 5
y_fov = 5
x_fov_pp = x_fov/len(np.unique(geom_hess5_xc_from_root)) #per pixel
y_fov_pp = y_fov/len(np.unique(geom_hess5_yc_from_root))
x_fov_px = x_fov/(np.max(geom_hess5_xc_from_root)-np.min(geom_hess5_xc_from_root))
y_fov_py = y_fov/(np.max(geom_hess5_yc_from_root)-np.min(geom_hess5_yc_from_root))
print(x_fov_pp, y_fov_pp)
print(x_fov_px, y_fov_py)


# In[64]:


get_ipython().run_cell_magic('time', '', '#reading zenith files:\nos.chdir("D:\\\\Masterarbeit ECAP\\\\First plots")\naz_zen_nrun = np.loadtxt("az_zen_file.txt", usecols=0, delimiter = ";")\naz_zen_telId = np.loadtxt("az_zen_file.txt", usecols=1, delimiter = ";")\naz_zen_az = np.loadtxt("az_zen_file.txt", usecols=2, delimiter = ";")\naz_zen_zen = np.loadtxt("az_zen_file.txt", usecols=3, delimiter = ";")\naz_zen_utc = np.genfromtxt("az_zen_file.txt", dtype=str, usecols=4, delimiter = ";")\nos.chdir("D:\\\\Masterarbeit ECAP\\\\First plots\\eval_data")')


# In[129]:


get_ipython().run_cell_magic('time', '', '#write to a file the interesting values such as \n#nrun, telId, az, zen, utc,\n#1.track IS times and points and mean brigtness, 2. track IS times and points and mean brigtness, ... \nnrun = nruns[1]\nprint(nrun)\ns = az_zen_utc[np.where(az_zen_nrun ==nrun)[0]][0]\nprint(s[6:])\nheader = ["nrun", "telId", "az", "zen", "utc first event in nrun",\n          "1. track X1 ", "1. track Y1", "1. track X2", "1. track Y2", \n          "1. track time", "1. track duration", "1. track mean brightness",\n          "2. track X1 ", "2. track Y1", "2. track X2", "2. track Y2", \n          "2. track time", "2. track duration", "2. track mean brightness",\n          "3. track X1 ", "3. track Y1", "3. track X2", "3. track Y2", \n          "3. track time", "3. track duration", "3. track mean brightness"]\nfile_params = open("parameters.csv", "w")\nparams_writer = csv.writer(file_params)\nparams_writer.writerow(header) \nfor i in range(len(az_zen_nrun)):\n    row = [int(az_zen_nrun[i]), int(az_zen_telId[i]), int(az_zen_az[i]), int(az_zen_zen[i]), az_zen_utc[i]]\n    for nrun in nruns:\n        for telId in range(1,6):\n\n            try:\n                \n                if (nrun == int(az_zen_nrun[i]) and int(az_zen_telId[i]) == telId):\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][0][0])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][0][1])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][1][0])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][1][1])\n                    print(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][1][1])\n                    row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-1"][0]))\n                    row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-1"][-1] - \n                               trails_high_avg[nrun][telId]["tracks"]["z-1"][0]))\n                    row.append(round(np.average(trails_high_avg[nrun][telId]["tracks"]["brightness-1"]),1))\n                    try:\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][0][0])\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][0][1])\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][1][0])\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][1][1])\n                        row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-2"][0]))\n                        row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-2"][-1] - \n                                       trails_high_avg[nrun][telId]["tracks"]["z-2"][0]))\n                        row.append(round(np.average(trails_high_avg[nrun][telId]["tracks"]["brightness-2"]),1))\n                    except:\n                        pass\n                    try:\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][0][0])\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][0][1])\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][1][0])\n                        row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][1][1])\n                        row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-3"][0]))\n                        row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-3"][-1] - \n                                       trails_high_avg[nrun][telId]["tracks"]["z-3"][0]))\n                        row.append(round(np.average(trails_high_avg[nrun][telId]["tracks"]["brightness-3"]),1))\n                    except:\n                        pass\n            except:\n                pass\n            \n    params_writer.writerow(row)\nfile_params.close()')


# In[153]:


get_ipython().run_cell_magic('time', '', '#write to a file the interesting values such as \n#nrun, telId, az, zen, utc,\n#1.track IS times and points and mean brigtness, 2. track IS times and points and mean brigtness, ... \nnrun = nruns[1]\nprint(nrun)\ns = az_zen_utc[np.where(az_zen_nrun ==nrun)[0]][0]\nprint(s[6:])\nheader_CT5 = ["nrun", "telId", "az", "zen", "utc first event in nrun",\n          "1. track X1 ", "1. track Y1", "1. track X2", "1. track Y2", \n          "1. track time", "1. track duration", "1. track mean brightness",\n          "2. track X1 ", "2. track Y1", "2. track X2", "2. track Y2", \n          "2. track time", "2. track duration", "2. track mean brightness",\n          "3. track X1 ", "3. track Y1", "3. track X2", "3. track Y2", \n          "3. track time", "3. track duration", "3. track mean brightness"]\nfile_params_CT5 = open("parameters_CT5.csv", "w")\nparams_CT5_writer = csv.writer(file_params_CT5)\nparams_CT5_writer.writerow(header_CT5) \nrow = []\nfor i in range(len(az_zen_nrun)):\n    lastrow = row\n    row = [int(az_zen_nrun[i]),int(az_zen_telId[i]), int(az_zen_az[i]), int(az_zen_zen[i]), az_zen_utc[i]]\n    for nrun in nruns:\n        try:\n\n            if (nrun == int(az_zen_nrun[i]) and int(az_zen_telId[i]) == 5):\n                print(nrun , int(az_zen_nrun[i]) , int(az_zen_telId[i]) , 5)\n                row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][0][0])\n                row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][0][1])\n                row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][1][0])\n                row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][1][1])\n                print(trails_high_avg[nrun][telId]["tracks"]["IS track 1"][0][1][1])\n                row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-1"][0]))\n                row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-1"][-1] - \n                           trails_high_avg[nrun][telId]["tracks"]["z-1"][0]))\n                row.append(round(np.average(trails_high_avg[nrun][telId]["tracks"]["brightness-1"]),1))\n                try:\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][0][0])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][0][1])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][1][0])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 2"][0][1][1])\n                    row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-2"][0]))\n                    row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-2"][-1] - \n                                   trails_high_avg[nrun][telId]["tracks"]["z-2"][0]))\n                    row.append(round(np.average(trails_high_avg[nrun][telId]["tracks"]["brightness-2"]),1))\n                except:\n                    pass\n                try:\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][0][0])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][0][1])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][1][0])\n                    row.append(trails_high_avg[nrun][telId]["tracks"]["IS track 3"][0][1][1])\n                    row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-3"][0]))\n                    row.append(int(trails_high_avg[nrun][telId]["tracks"]["z-3"][-1] - \n                                   trails_high_avg[nrun][telId]["tracks"]["z-3"][0]))\n                    row.append(round(np.average(trails_high_avg[nrun][telId]["tracks"]["brightness-3"]),1))\n                except:\n                    pass\n                \n        except:\n            pass\n    try:\n        if row[1] ==5:\n            print(row)\n            params_CT5_writer.writerow(row)\n    except:\n        pass\nfile_params_CT5.close()')

