#load camera geometry from txt file
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
import h5py
import time
from astropy.time import Time
import astropy.units as u
from numba import jit,njit,int32,float32
from numba.typed import List
global flash_geom_id
global flash_geom_x
global flash_geom_y
global min_x_dif
global min_y_dif
global nn_pix
flash_geom_id = np.loadtxt("D:/Masterarbeit ECAP/Trail finder/flashcam_geometry.txt", usecols = 0, delimiter=";")
flash_geom_id = list(flash_geom_id)
for i in range(len(flash_geom_id)):
    flash_geom_id[int(i)] = int(i)
flash_geom_id = np.array(flash_geom_id)
flash_geom_x = np.loadtxt("D:/Masterarbeit ECAP/Trail finder/flashcam_geometry.txt", usecols = 1, delimiter=";")
flash_geom_y = np.loadtxt("D:/Masterarbeit ECAP/Trail finder/flashcam_geometry.txt", usecols = 2, delimiter=";")
a = sorted(flash_geom_x)
b = sorted(flash_geom_y)
for i in range(len(a)-1):
    if a[i+1]-a[i]>0:
        min_x_dif = round(a[i+1]-a[i], 5)
        #print(min_x_dif)
        break
for i in range(len(b)-1):
    if b[i+1]-b[i]>0:
        min_y_dif = round(b[i+1]-b[i], 5)
        #print(min_y_dif)
        break
def lin_fit(x, m, t):
    return m*x+t
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
def plot_scatter(x, y, z):
    #fig= plt.figure(figsize=(6,5), )
    plt.scatter(x,y,c=z, cmap = "jet")
    plt.xlim(-1.25,1.25)
    plt.ylim(-1.25,1.25)
    plt.xlabel("Camera x-pos")
    plt.ylabel("Camera y-pos")
    plt.title("Run {}".format(nrun)+", CT {}".format(telId))
    #cbar = plt.colorbar()
    #cbar.set_label("Time  in run [s]")
def plot_from_pix(pix, z):
    x = []
    y = []
    for i in range(len(pix)):
        x.append(flash_geom_x[int(pix[i])])
        y.append(flash_geom_y[int(pix[i])])
    plt.scatter(x,y, c = z, cmap = "jet")
    plt.xlim(-1.25,1.25)    
    plt.ylim(-1.25,1.25)  
    plt.xlabel("Camera x-pos")
    plt.ylabel("Camera y-pos")
    #plt.title("Run {}".format(nrun)+", CT {}".format(5))
    cbar = plt.colorbar()
    cbar.set_label("Time  in run [s]")
    return
def draw_box():
    mask_x_top = np.around(np.array(flash_geom_x), 6) == np.max(np.around(flash_geom_x, 6))
    mask_x_bot = np.around(np.array(flash_geom_x), 6) == np.min(np.around(flash_geom_x, 6))
    mask_y_top = np.around(np.array(flash_geom_y), 6) == np.max(np.around(flash_geom_y, 6))
    mask_y_bot = np.around(np.array(flash_geom_y), 6) == np.min(np.around(flash_geom_y, 6))
    
    y_secondtolast = sorted(np.unique(np.around(flash_geom_y, 6)))[-2]
    y_secondtofirst = sorted(np.unique(np.around(flash_geom_y, 6)))[1]
    mask1 = np.around(np.array(flash_geom_y), 6) == y_secondtolast
    mask2 = np.around(np.array(flash_geom_y), 6) == y_secondtofirst
    
    box_x = [flash_geom_x[mask_x_top][0],
             np.max(np.around(flash_geom_x[mask1], 6)),
            np.max(flash_geom_x[mask_y_top]),
            np.min(flash_geom_x[mask_y_top]),
            flash_geom_x[mask_x_bot][0],
            flash_geom_x[mask_x_bot][-1],
            np.min(flash_geom_x[mask_y_bot]),
            np.max(flash_geom_x[mask_y_bot]),
             np.max(np.around(flash_geom_x[mask2], 6)),
            flash_geom_x[mask_x_top][0],
            flash_geom_x[mask_x_top][0]]
    box_y = [flash_geom_y[mask_x_top][0],
             y_secondtolast,
             flash_geom_y[mask_y_top][0],
             flash_geom_y[mask_y_top][-1],
             np.max(flash_geom_y[mask_x_bot]),
             np.min(flash_geom_y[mask_x_bot]),
             flash_geom_y[mask_y_bot][0],
             flash_geom_y[mask_y_bot][-1],
             y_secondtofirst,
             flash_geom_y[mask_x_top][0],
             flash_geom_y[mask_x_top][0]]
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
            #idx = np.argwhere(np.diff(np.sign(m*x_vals+t - y_vals))).flatten()
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
            #idx = np.argwhere(np.diff(np.sign(m*np.array(x_vals)+t - y_vals))).flatten()
#             try:
#                 if len(idx) !=0:
#                     idx = idx[0]
#                     intersections.append([x_vals[idx][0], y_vals[idx][0]])
#             except:
#                 print(idx)
#     intersections = np.around(intersections, 4)
#     print("intersections are", intersections[0], ", ", intersections[1])
    return #intersections
def get_x_y_from_pix(pix):
    x = []
    y = []
    for i in range(len(pix)):
        x.append(flash_geom_x[int(pix[i])])
        y.append(flash_geom_y[int(pix[i])])
    return x, y



def read_data_sets(firstrun, lastrun):
    if (firstrun%200!=0  or lastrun%200!=0):
        print("Numbers divisible by 200 needed!")
        return
    if firstrun> lastrun:
        print("Second number needs to be larger than first one!")
        return
    run_set = 0
    new_nruns = []
    new_high_pixel = {}
    os.chdir("D:\\Masterarbeit ECAP\\First plots\\eval_data\\high_pixel")
    while (run_set<(lastrun-firstrun)/200):
        runs_directory = "run"+str(firstrun+run_set*200)+"-"+str(firstrun+run_set*200+199)
        try:
            os.chdir("D:\\Masterarbeit ECAP\\First plots\\eval_data\\high_pixel\\"+runs_directory)
        except:
            print("folder","\""+runs_directory+"\"", "does not exist")
            run_set+=1
            continue
        run_numbers = os.listdir(os.curdir)
        for i in range(len(run_numbers)):
            try:
                new_nruns.append(int(run_numbers[i]))
                try:
                    print(os.getcwd())
#                     print(str(run_numbers[i])+"\\high_pix_"+str(run_numbers[i]))#+"_CT_"+str(telId)+".h5")
                    hfile = h5py.File(str(run_numbers[i])+"\\high_pix_"+str(run_numbers[i])+"_CT_"+str(5)+".h5", "r") 
                    print("Reading high_pix_"+str(run_numbers[i])+"_CT_"+str(5)+".h5")
                    new_high_pixel[int(run_numbers[i])] = {}
                    new_high_pixel[int(run_numbers[i])]["pix"] = hfile["Pix ID"]
                    new_high_pixel[int(run_numbers[i])]["brightness"] = hfile["Brightness"]
                    new_high_pixel[int(run_numbers[i])]["time"] = hfile["Time"]
                    hfile.close()
                    print("Reading successful")
                except:
                    print("No file named high_pix_"+str(run_numbers[i])+"_CT_"+str(5)+".h5 found in folder", run_numbers[i])
            except:
                print("No runs in folder/not a folder:", run_numbers[i])
        
        run_set+=1
    return new_high_pixel

print("Function names are:")
print("flash_geom_x , flash_geom_y")
print("lin_fit(x, m, t)")
print("intersection(lst1, lst2)")
print("plot_scatter(x, y, z)")
print("plot_from_pix(pix, z), uncomment in functions_import.py for title with nrun")
print("draw_box()")
print("get_x_y_from_pix(pix), output x,y")
print("plot_lin_fit(track)")
print("read_data_sets(firstrun, lastrun)")
print("apply_selection_cut(dict_high_pixel)")

def get_nn_pix(size):


    nn_pix = []
    for pix in range(len(flash_geom_x)):
        mask_x = np.logical_and(np.around(flash_geom_x,4)<
                                min_x_dif*size*0.6+np.around(flash_geom_x,4)[pix],
                                np.around(flash_geom_x,4)>
                                -min_x_dif*size*0.6+np.around(flash_geom_x,4)[pix])

        mask_xy = np.logical_and(np.around(flash_geom_y,4)[mask_x]<
                                 min_y_dif*(size-0.1)+np.around(flash_geom_y,4)[pix],
                                 np.around(flash_geom_y,4)[mask_x]>
                                 -min_y_dif*size+np.around(flash_geom_y,4)[pix])


        st_x = set(np.around(flash_geom_x,4)[mask_x][mask_xy])
        st_y = set(np.around(flash_geom_y,4)[mask_x][mask_xy])
        indices_x = np.array([j for j, e in enumerate(np.around(flash_geom_x,4)) if e in st_x])
        indices_y = np.array([j for j, e in enumerate(np.around(flash_geom_y,4)) if e in st_y])

        indices_xy = np.array([j for j, e in enumerate(indices_x) if e in indices_y ])
        nn_pix.append(list(indices_x[indices_xy]))
    return nn_pix


nn_pix = get_nn_pix(10)
print("get_nn_pix(size) mit nn_pix, bzw. typed_nn_pix for @jit")

def plot_lin_fit(track):
    if track[0][0][0] == -1:
        return 
    x_box = np.arange(-1.25, 1.25, 0.1)
    for i in range(len(track)):
        pix_x, pix_y = get_x_y_from_pix(track[i][0])
        popt, pcov = curve_fit(lin_fit, pix_x, pix_y)
        print(popt)
        plt.plot(x_box, lin_fit(x_box, popt[0], popt[1]))
        plt.scatter(pix_x, pix_y)
        plt.xlim(-1.25, 1.25)
        plt.ylim(-1.25, 1.25)
    plt.show
    return

def apply_selection_cut(dict_high_pixel):
    max_counts = 70
    high_pixel_cut = {}
    for nrun in dict_high_pixel.keys():
        high_pixel_cut[nrun] = {}
        values, counts = np.unique(dict_high_pixel[nrun]["pix"], return_counts=True)
        mask_max_counts = counts<max_counts       
        st = set(values[mask_max_counts])
        result = [i for i, e in enumerate(dict_high_pixel[nrun]["pix"]) if e in st]
        high_pixel_cut[nrun]["pix"] = np.array(dict_high_pixel[nrun]["pix"])[result].astype(int)
        high_pixel_cut[nrun]["brightness"] = np.array(dict_high_pixel[nrun]["brightness"])[result]
        high_pixel_cut[nrun]["time"] = np.array(dict_high_pixel[nrun]["time"])[result]-dict_high_pixel[nrun]["time"][0]
    return high_pixel_cut

