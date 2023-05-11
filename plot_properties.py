import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import h5py
from datetime import date
from datetime import datetime
from astropy.time import Time
import astropy.units as u
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize':'x-large'}
pylab.rcParams.update(params)
mpl.rcParams["figure.dpi"]= 300



path_to_prop = "/home/tholang/trail_properties/"
path_to_plots = "/home/tholang/plots/"
hfile = h5py.File(path_to_prop+"run_properties.h5", "r")
props = list(hfile.keys())
prop_dict = {}
mask_time = np.array(hfile["t0"])>1e09 #Removes all values of t0=0 
for key in props:
    prop_dict[key] = np.array(hfile[key])
    prop_dict[key] = prop_dict[key][mask_time]

index_first_run = []
for i in range(len(prop_dict["nrun"])):
    if prop_dict["nrun"][i] not in prop_dict["nrun"][:i]:
        index_first_run.append(i)

index_first_run = np.array(index_first_run)
mask_N = prop_dict["N_track"]>-1
t = Time(prop_dict["t0"], format = "unix")
t.format = "decimalyear"
t = t.value
TIN = prop_dict["time_in_night"]/3600
for i in range(len(TIN)):
    if TIN[i]>12:
        TIN[i] = TIN[i]-24

table = np.genfromtxt("Starlink_launches_from_source.CSV", 
                      delimiter = ";", skip_header = 2)
table = table.T
header =np.genfromtxt("Starlink_launches_from_source.CSV", delimiter = ";", 
                      skip_header = 1, max_rows = 1, dtype = "str")
COSPAR_ID = np.loadtxt("Starlink_launches_from_source.CSV", delimiter = ";",
                       dtype = "str", usecols = 0)
COSPAR_ID = [COSPAR_ID[i][np.char.find(COSPAR_ID[i], "("):] 
             for i in range(len(COSPAR_ID))]
COSPAR_ID = [COSPAR_ID[i][np.char.find(COSPAR_ID[i], ", 20")+2:-1] 
             for i in range(len(COSPAR_ID))]
launch_dates = np.loadtxt("Starlink_launch_dates.CSV", delimiter = ";",dtype = "str", 
                          usecols = 0, skiprows = 1)
launch_dates_all_sats = np.loadtxt("Launch dates.txt", dtype = "str",
                                   usecols = 0, skiprows = 1)
launch_dates_all_sats = [np.array(launch_dates_all_sats[i].split(".")).astype(int) 
                         for i in range(len(launch_dates_all_sats))]
launch_dates_all_sats = [datetime(launch_dates_all_sats[i][2], 
                                  launch_dates_all_sats[i][1], 
                                  launch_dates_all_sats[i][0]) 
                         for i in range(len(launch_dates_all_sats))]
dates_all_sats = Time(launch_dates_all_sats, format= "datetime")

format_date = "%d %B %Y, %H:%M"
format_date_sec = "%d %B %Y, %H:%M:%S"
dates = []
for i in range(len(launch_dates)):
    try:
        date = datetime.strptime(launch_dates[i] ,format_date)
        dates.append(date)
        continue
    except:
        pass
    try:
        date = datetime.strptime(launch_dates[i] ,format_date_sec)
    except:
        date = datetime(1900,1,1,0)
    dates.append(date)

dates = Time(dates, format = "datetime")
dates = dates.value

cumul_sum_launched = [np.nansum(table[1][:i+1]) for i in range(len(table[1]))]
cumul_sum_orbit = [np.nansum(table[5][:i+1]) for i in range(len(table[5]))]




launch_SL = open("Launch_dates_starlink.txt", "r").read().split("\n")
launch_OW = open("Launch_dates_OneWeb.txt", "r").read().split("\n")
launch_ES = open("Launch_dates_E-Space.txt", "r").read().split("\n")
format_launch = "%d.%m.%Y"
launch_all = launch_SL + launch_OW + launch_ES

launch_all = Time(np.sort([datetime.strptime(i, format_launch) for i in launch_all]),
                  format = "datetime")
mask_dates_all_sats = launch_all >np.sort(dates_all_sats.value)[-1]



TIN_bins = np.arange(int(np.min(TIN[mask_N]))-1, int(np.max(TIN[mask_N])+2))
zenith_bins = np.arange(0,61,5)
dur_bins = np.logspace(np.log10(np.min(prop_dict["duration"][mask_N])), np.log10(np.max(prop_dict["duration"])), 20 )
bright_bins = np.logspace(np.log10(900), np.log10(np.max(prop_dict["avg_brightness"])), 20 )
pix_bins = np.logspace(np.log10(np.min(prop_dict["unique_pix"][mask_N])),
                       np.log10(np.max(prop_dict["unique_pix"])), 20 )


mask_TIN = np.logical_or(TIN<-5.5, TIN>2.5)
mask_N_TIN = np.logical_or(TIN[mask_N]<-5.5, TIN[mask_N]>2.5)

figure_size = (6,3)


fig1, ax1 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
plt.rcParams['date.converter'] = 'concise'
ax1.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=12, maxticks=20))
N_runs,b = np.histogram(np.unique(t), bins = np.arange(2019,2023.5,0.08333333))# Runs per month
N_trails,d = np.histogram(t[mask_N], bins = np.arange(2019,2023.5,0.08333333))# Trails per month
mask_zeros = N_trails >0
N_trails = N_trails[mask_zeros]
N_runs = N_runs[mask_zeros]
b = b[:-1][mask_zeros]
b = Time(b, format = "decimalyear")
b.format = "datetime"
b = b.value
err = N_trails/N_runs/np.sqrt(N_trails)
ax1.scatter(b[:-1], N_trails[:-1]/N_runs[:-1], s = 5)
ax1.errorbar(b[:-1], N_trails[:-1]/N_runs[:-1], yerr = err[:-1], fmt=".k", label = "Trails per run")
ax1.set_xlabel("Time [year]")
ax1.set_ylabel("Number of Trails [counts/run]")
#ax1.legend(loc = "upper left")
ax1.set_xlim([datetime(2019, 9,1), datetime(2023, 2, 28)])
#ax1.set_ylim(0, 1.2)
ax1_2 = ax1.twinx()
#ax1_2.plot(dates[45:], np.array(cumul_sum_orbit[45:])+len(dates_all_sats)-cumul_sum_orbit[45], label = "Starlink in orbit", c = "red")
ax1_2.plot(np.sort(dates_all_sats.value), np.arange(1, len(dates_all_sats)+1), 
           label = "Satellites in orbit", c = "red")
ax1_2.plot(launch_all[mask_dates_all_sats].value,
           np.arange(len(dates_all_sats.value),
                     len(launch_all[mask_dates_all_sats])+len(dates_all_sats)),
           label = "Starlink, OneWeb and E-Space", c = "red", ls = "dotted")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax1_2.get_legend_handles_labels()
ax1_2.legend(lines+lines2, labels +labels2, loc = 0)
ax1_2.set_ylabel("Number of Satellites")
ax1_2.set_ylim(0,15000)
#ax1.set_title("All runs")
fig1.savefig(path_to_plots+"Avg_trails_all_runs.jpg")
plt.close()

counter_of_plots = 0
counter_of_plots+=1
print(counter_of_plots,"Plots created")

fig1_1, ax1_1 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
plt.rcParams['date.converter'] = 'concise'
ax1_1.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=12, maxticks=20))
N_runs,b = np.histogram(np.unique(t[mask_N]), bins = np.arange(2019,2023.5,0.08333333))# Runs per month
N_trails,d = np.histogram(t[mask_N], bins = np.arange(2019,2023.5,0.08333333))# Trails per month
mask_zeros = N_trails >0
N_trails = N_trails[mask_zeros]
N_runs = N_runs[mask_zeros]
b = b[:-1][mask_zeros]
b = Time(b, format = "decimalyear")
b.format = "datetime"
b = b.value
err = N_trails/N_runs/np.sqrt(N_trails)
ax1_1.scatter(b, N_trails/N_runs, s = 5)
ax1_1.errorbar(b, N_trails/N_runs, yerr = err, fmt=".k")
ax1_1.set_xlabel("Time [year]")
ax1_1.set_ylabel("Number of Trails [counts/run]") 
#ax1_1.set_title("Only runs with trails")  
fig1_1.savefig(path_to_plots+"Avg_trails_runs_with_trails.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


'''
fig2, ax2 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
c,d,e = np.histogram2d(t, prop_dict["N_track"], 
                       bins = [np.arange(2019,2023.5, 0.08333333), np.arange(0,np.max(prop_dict["N_track"])+1)])
X,Y = np.meshgrid(d,e)
for i in range(len(c)):
    c[i]=c[i]/a_new[i]
c=c.T
plt.pcolormesh(X,Y,c, cmap = "gist_heat_r")
plt.xlabel("Year")
plt.ylabel("N-tracks")
plt.colorbar(label = "Runs with N tracks per nruns in month")
plt.title("")
plt.xlim((2019.75,2022.75))
fig2.savefig(path_to_plots+"n_tracks_hist2d.jpg")
plt.close()
'''

fig3, ax3 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
dur_bins = np.logspace(np.log10(np.min(prop_dict["duration"][mask_N])), np.log10(np.max(prop_dict["duration"])), 20 )
a,b = np.histogram(prop_dict["duration"][mask_N], bins =dur_bins)
bin_avg = (b[:-1]+b[1:])/2
#ax3.scatter(bin_avg,a)
for i in range(len(bin_avg)):
    ax3.hlines(y = a[i], xmin = b[i], xmax = b[i+1], color = "black")

ax3.errorbar(bin_avg,a, yerr = np.sqrt(a), fmt=".k", color = "black" )
ax3.set_xlabel("Duration [s]")
ax3.set_ylabel("Number of Trails [counts]")
#ax3.set_title("Histogram durations")
xticks = np.array([5,10,20,50,100,200])
ax3.set_xscale("log")
ax3.set_xticks(xticks)
xticks_label = xticks.astype("str")
ax3.set_xticklabels(xticks_label)
fig3.savefig(path_to_plots+"hist_durations.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig4, ax4 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
v_max = np.max(prop_dict["velocity"])
v_step = v_max/len(dur_bins)
v_bins = np.arange(0,v_max+0.05, 0.05)
a,b = np.histogram(prop_dict["velocity"], bins = v_bins)
bin_avg = (b[:-1]+b[1:])/2
#ax4.scatter(b[:-1], a)
ax4.errorbar(bin_avg, a, yerr = np.sqrt(a), xerr = round((bin_avg[1]-bin_avg[0])/2,3), fmt = ".k")
ax4.set_xlabel("Velocity [$^\circ$/s]")
ax4.set_ylabel("Number of Trails [counts]")
#ax4.set_title("Histogram velocities")
fig4.savefig(path_to_plots+"hist_velocity.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")

fig4_2, ax4_2 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
v_max = np.max(prop_dict["velocity"])
v_step = v_max/(len(dur_bins)-1)
v_bins = np.arange(0,v_max+v_step, v_step)
a,b,c = np.histogram2d(prop_dict["duration"][mask_N], prop_dict["velocity"][mask_N],
                       bins = [dur_bins, v_bins])
X,Y = np.meshgrid(b,c)
a = a.T
ax4_2.pcolormesh(X,Y,a)#, cmap = "gist_heat_r")
ax4_2.set_xlabel("Duration [s]")
ax4_2.set_ylabel("Velocity [$^\circ$/s]")
#ax4_2.set_title("Duration-velocity distribution")
xticks = np.array([5,10,20,50,100,200])
ax4_2.set_xscale("log")
ax4_2.set_xticks(xticks)
xticks_label = xticks.astype("str")
ax4_2.set_xticklabels(xticks_label)
fig4_2.colorbar(ax4_2.pcolormesh(X,Y,a), label = "Number of Trails [counts]")
fig4_2.savefig(path_to_plots+"hist_2D_dur_v.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig4_3, ax4_3 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
v_max = np.max(prop_dict["velocity"])
v_step = v_max/(len(dur_bins)-1)
v_bins = np.arange(0,v_max+v_step, v_step)
pix_bins = np.logspace(np.log10(np.min(prop_dict["unique_pix"][mask_N])), 
                       np.log10(np.max(prop_dict["unique_pix"])), 20 )
a,b,c = np.histogram2d(prop_dict["unique_pix"][mask_N], prop_dict["velocity"][mask_N],
                       bins = [pix_bins, v_bins])
X,Y = np.meshgrid(b,c)
a = a.T
ax4_3.pcolormesh(X,Y,a)#, cmap = "gist_heat_r")
ax4_3.set_xlabel("Unique Pixels [counts]")
ax4_3.set_ylabel("Velocity [m/s]") 
xticks = np.array([5,10,20,50,100,200])
ax4_3.set_xscale("log")
ax4_3.set_xticks(xticks)
xticks_label = xticks.astype("str")
ax4_3.set_xticklabels(xticks_label)
fig4_3.colorbar(ax4_3.pcolormesh(X,Y,a), label = "Number of Trails [counts]")
fig4_3.savefig(path_to_plots+"hist_2D_pix_v.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")

fig4_4, ax4_4 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
pix_bins = np.logspace(np.log10(np.min(prop_dict["unique_pix"][mask_N])), 
                       np.log10(np.max(prop_dict["unique_pix"])), 20 )
a,b,c = np.histogram2d(prop_dict["unique_pix"][mask_N], prop_dict["duration"][mask_N],
                       bins = [pix_bins, dur_bins])
X,Y = np.meshgrid(b,c)
a = a.T
ax4_4.pcolormesh(X,Y,a)#, cmap = "gist_heat_r")
ax4_4.set_xlabel("Unique pixels [counts]")
ax4_4.set_ylabel("Duration [s]")
xticks = np.array([5,10,20,50,100,200])
ax4_4.set_xscale("log")
ax4_4.set_xticks(xticks)
ax4_4.set_yscale("log")
ax4_4.set_yticks(xticks)
xticks_label = xticks.astype("str")
ax4_4.set_yticklabels(xticks_label)
ax4_4.set_xticklabels(xticks_label)
fig4_4.colorbar(ax4_4.pcolormesh(X,Y,a), label = "Number of Trails [counts]")
fig4_4.savefig(path_to_plots+"hist_2D_pix_dur.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig4_5, ax4_5 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
v_max = np.max(prop_dict["velocity"])
v_step = v_max/(len(dur_bins)-1)
v_bins = np.arange(0,v_max+v_step, v_step)
bright_bins_lin = np.arange(900, np.max(prop_dict["avg_brightness"]), 200) 
a,b,c = np.histogram2d(prop_dict["duration"][mask_N], prop_dict["avg_brightness"][mask_N],
                       bins = [dur_bins, bright_bins])
X,Y = np.meshgrid(b,c)
a = a.T
ax4_5.pcolormesh(X,Y,a)#, cmap = "gist_heat_r")
ax4_5.set_xlabel("Duration [s]")
ax4_5.set_ylabel("Average Brightness [MHz]")
#ax4_2.set_title("Duration-velocity distribution")
xticks = np.array([5,10,20,50,100,200])
ax4_5.set_xscale("log")
ax4_5.loglog()
ax4_5.set_xticks(xticks)
xticks_label = xticks.astype("str")
ax4_5.set_xticklabels(xticks_label)
yticks = np.array([1000,2000, 5000])
ax4_5.set_yticks(yticks)
yticks_label = yticks.astype("str")
ax4_5.set_yticklabels(yticks_label)
ax4_5.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax4_5.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())                                                                                            
ax4_5.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax4_5.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
fig4_5.colorbar(ax4_5.pcolormesh(X,Y,a), label = "Number of Trails [counts]")
fig4_5.savefig(path_to_plots+"hist_2D_dur_bright.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


'''
#Width is not implemented yet
fig5 = plt.figure(figsize=figure_size)
bins_w = np.arange(0,np.max(prop_dict["width"]), 0.01)
plt.hist(prop_dict["width"], bins =bins_w)
plt.xlabel("w [pix]")
plt.ylabel("N")
plt.xlim((0,np.max(bins_w)))
plt.title("Histogram Widths")
fig5.savefig(path_to_plots+"hist_velocity.jpg")
plt.close()
'''

fig6, ax6 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
pix_bins = np.logspace(np.log10(np.min(prop_dict["unique_pix"][mask_N])), 
                       np.log10(np.max(prop_dict["unique_pix"])), 20 )
a,b = np.histogram(prop_dict["unique_pix"], bins =pix_bins)
bin_avg = (b[:-1]+b[1:])/2
ax6.errorbar(bin_avg, a, yerr = np.sqrt(a),fmt = ".k")
for i in range(len(bin_avg)):
    ax6.hlines(y = a[i], xmin = b[i], xmax = b[i+1], color = "black")

ax6.set_xlabel("Unique Pixels [counts]")
ax6.set_ylabel("Number of Trails [counts]")
ax6.set_xscale("log")
#ax6.set_title("Histogram pixels per trail")   
xticks = np.array([5,10,20,50,100,200])
ax6.set_xscale("log")
ax6.set_xticks(xticks)
xticks_label = xticks.astype("str")
ax6.set_xticklabels(xticks_label)
fig6.savefig(path_to_plots+"hist_unique_pix.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


'''
fig7, ax7 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
TIN = prop_dict["time_in_night"]/3600
for i in range(len(TIN)):
    if TIN[i]>12:
        TIN[i] = TIN[i]-24

TIN_bins = np.arange(int(np.min(TIN[mask_N]))-1, int(np.max(TIN[mask_N])+2))
a,b = np.histogram(prop_dict["zenith"][index_first_run], bins = np.arange(0,np.max(prop_dict["zenith"]),5))
#a accounts for every run, whether a trail was detected or not
c,d,e = np.histogram2d(prop_dict["zenith"][mask_N], TIN[mask_N],
                       bins = [zenith_bins, TIN_bins])
#c takes only Trails [counts] into account 
a_2, b_2 = np.histogram(TIN[index_first_run], bins = TIN_bins)
X,Y = np.meshgrid(d,e)
for i in range(len(c)):
    c[i] = c[i]/a[i] #Number of Trails [counts] per zenith angle range
    #c[i] = c[i]/a_2[i] #Number of Trails [counts] per time in night

c=c.T
im7 = ax7.pcolormesh(X,Y,c)#, cmap = "gist_heat_r")
ax7.set_xlabel("Zenith [deg]")
ax7.set_ylabel("Time_in_night [h]")
ax7.set_yticks(np.arange(-6,6,2))
ax7.set_yticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
fig7.colorbar(im7, label = "Average number of Trails [counts]")
#ax7.set_title("2D-histogram of time in night vs. average number of Trails [counts]")
fig7.savefig(path_to_plots+"time_in_night_zenith_hist2d_relative.jpg")
plt.close()


fig7_2, ax7_2 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
TIN = prop_dict["time_in_night"]/3600
for i in range(len(TIN)):
    if TIN[i]>12:
        TIN[i] = TIN[i]-24

TIN_bins = np.arange(int(np.min(TIN[mask_N]))-1, int(np.max(TIN[mask_N])+2))
a,b = np.histogram(prop_dict["zenith"][index_first_run], bins = np.arange(0,np.max(prop_dict["zenith"]),5))
#a accounts for every run, whether a trail was detected or not
c,d,e = np.histogram2d(prop_dict["zenith"][mask_N], TIN[mask_N],
                           bins = [zenith_bins, TIN_bins])
#c takes only Trails [counts] into account
c_err = np.sqrt(c) #Poissonian error
X,Y = np.meshgrid(d,e)
for i in range(len(c_err)):
    c_err[i] = c_err[i]/a[i] #Number of Trails [counts] per zenith angle range

c_err=c_err.T
im7_2 = ax7_2.pcolormesh(X,Y,c_err)#, cmap = "gist_heat_r")
ax7_2.set_xlabel("Zenith [deg]")
ax7_2.set_ylabel("Time_in_night [h]")
ax7_2.set_yticks(np.arange(-6,6,2))
ax7_2.set_yticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
fig7_2.colorbar(im7_2, label = "Error of number of Trails [counts] per zenith range")
ax7_2.set_title("2D-histogram of errors of time in night vs. average number of Trails [counts]")
fig7_2.savefig(path_to_plots+"time_in_night_zenith_hist2d_relative_err.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")

'''


#Idea: Find entry indices corresponding to each normalized bin
fig8_03, ax8_03 = plt.subplots(1,1,figsize = figure_size, constrained_layout = True)
a,b,c = np.histogram2d(prop_dict["zenith"][mask_N], TIN[mask_N],
                           bins = [zenith_bins, TIN_bins])
X,Y = np.meshgrid(b,c)
im8_03 = ax8_03.pcolormesh(X,Y, a.T)
ax8_03.set_xlabel("Zenith angle [$^\circ$]")
ax8_03.set_ylabel("Time in night [UTC]")
ax8_03.set_yticks(np.arange(-6,6,2))
ax8_03.set_yticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
fig8_03.colorbar(im8_03, label = "Number of Trails [counts]")
fig8_03.savefig(path_to_plots+"zenith_TIN_trails.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig8_03_2, ax8_03_2 = plt.subplots(1,1,figsize = figure_size, constrained_layout = True)
a,b,c = np.histogram2d(prop_dict["zenith"][index_first_run], TIN[index_first_run],
                           bins = [zenith_bins, TIN_bins])
X,Y = np.meshgrid(b,c)
im8_03_2 = ax8_03_2.pcolormesh(X,Y, a.T)
ax8_03_2.set_xlabel("Zenith angle [$^\circ$]")
ax8_03_2.set_ylabel("Time in night [UTC]")
ax8_03_2.set_yticks(np.arange(-6,6,2))
ax8_03_2.set_yticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
fig8_03_2.colorbar(im8_03_2, label = "Number of runs")
fig8_03_2.savefig(path_to_plots+"zenith_TIN_runs.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")




fig8_04, ax8_04 = plt.subplots(1,1,figsize = figure_size, constrained_layout = True)
a,b,c = np.histogram2d(prop_dict["zenith"][mask_N], TIN[mask_N],
                           bins = [zenith_bins, TIN_bins])
x,y,z = np.histogram2d(prop_dict["zenith"][index_first_run], TIN[index_first_run],
                           bins = [zenith_bins, TIN_bins])
X,Y = np.meshgrid(b,c)
im8_04 = ax8_04.pcolormesh(X,Y, (a/x).T, norm = mpl.colors.LogNorm(), cmap = "viridis")
ax8_04.set_xlabel("Zenith angle [$^\circ$]")
ax8_04.set_ylabel("Time in night [UTC]")
ax8_04.set_yticks(np.arange(-6,6,2))
ax8_04.set_yticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
cbar8_04 = fig8_04.colorbar(im8_04, label = "Number of Trails [counts/run]", 
                            ticks = [0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5])
cbar8_04.ax.set_yticklabels(np.array([0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]).astype(str))
fig8_04.savefig(path_to_plots+"zenith_TIN_normed.jpg")
plt.close()


counter_of_plots+=1
print(counter_of_plots,"Plots created")




fig8, ax8 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
ax8.set_xticks(np.arange(-6,6,2))
ax8.set_xticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
ax8.hist(TIN[mask_N], bins = TIN_bins)
ax8.set_xlabel("Trail time [UTC]")
ax8.set_ylabel("Number of trails")
#ax8.set_title("Number of tracks")
fig8.savefig(path_to_plots+"TIN_hist_trails.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig8_2, ax8_2 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
ax8_2.set_xticks(np.arange(-6,6,2))
ax8_2.set_xticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
ax8_2.hist(TIN[index_first_run], bins = TIN_bins)
ax8_2.set_xlabel("Time of beginning of run")
ax8_2.set_ylabel("Number of runs")
#ax8_2.set_title("Number of runs per hour of night") 
fig8_2.savefig(path_to_plots+"TIN_hist_runs.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig8_3, ax8_3 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
a,b = np.histogram(TIN[index_first_run], bins = TIN_bins)
c,d = np.histogram(TIN[mask_N], bins = TIN_bins)
bin_avg = (b[:-1]+b[1:])/2
ax8_3.errorbar(bin_avg, c/a, xerr = 0.5 , yerr = np.sqrt(c)/a, fmt = ".k")
ax8_3.set_xticks(np.arange(-6,6,2))
ax8_3.set_xticklabels(["6$\,$pm", "8$\,$pm", "10$\,$pm", "12$\,$am", "2$\,$am", "4$\,$am"])
ax8_3.set_xlabel("Time of beginning of run")
ax8_3.set_ylabel("Average number of trails")
#ax8_3.set_title("Number of runs per hour of night")
fig8_3.savefig(path_to_plots+"TIN_hist_average.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig10, ax10 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
ax10.hist(prop_dict["zenith"][mask_N], bins = zenith_bins)
ax10.set_xlabel("Zenith [$^\circ$]")
ax10.set_ylabel("Number of trails")
#ax10.set_title("Number of trails")
fig10.savefig(path_to_plots+"zenith_hist_trails.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig10_2, ax10_2 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
ax10_2.hist(prop_dict["zenith"][index_first_run], bins = zenith_bins)
ax10_2.set_xlabel("Zenith [$^\circ$]")
ax10_2.set_ylabel("Number of runs")
#ax10_2.set_title("Number of runs")
fig10_2.savefig(path_to_plots+"zenith_hist_runs.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig11, ax11 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
a,b = np.histogram(prop_dict["zenith"][mask_N], bins = zenith_bins)
c,d = np.histogram(prop_dict["zenith"][index_first_run], bins = zenith_bins)
bin_avg = (b[:-1]+b[1:])/2
ax11.errorbar(bin_avg, a/c, xerr = 2.5, yerr = np.sqrt(a)/c, fmt=".k")
ax11.set_xlabel("Zenith [$^\circ$]")
ax11.set_ylabel("Average number of trails")
#ax11.set_title("Average number of trails per run")
fig11.savefig(path_to_plots+"zenith_hist_average.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")



fig12, ax12 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
hist_zenith_dur = []
for i in range(len(zenith_bins)-1):
    mask_zenith_bin = np.logical_and(prop_dict["zenith"][mask_N]>=zenith_bins[i], prop_dict["zenith"][mask_N]<zenith_bins[i+1])
    a,b = np.histogram(prop_dict["duration"][mask_N][mask_zenith_bin], bins =dur_bins)
    hist_zenith_dur.append(a)

X,Y = np.meshgrid(zenith_bins, dur_bins)
im12 = ax12.pcolormesh(X,Y,np.array(hist_zenith_dur).T)
ax12.set_xlabel("Zenith [$^\circ$]")
ax12.set_ylabel("Duration [s]")
ax12.set_yscale("log")
yticks = np.array([5,10,20,50,100])
ax12.set_yticks(yticks)
yticks_label = yticks.astype("str")
ax12.set_yticklabels(yticks_label)
fig12.colorbar(im12, label = "Number of Trails [coucts]")
fig12.savefig(path_to_plots+"hist_2D_zenith_duration.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig13, ax13 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
hist_zenith_dur = []
c,d = np.histogram(prop_dict["zenith"][index_first_run], bins = zenith_bins)
for i in range(len(zenith_bins)-1):
    mask_zenith_bin = np.logical_and(prop_dict["zenith"][mask_N]>=zenith_bins[i], prop_dict["zenith"][mask_N]<zenith_bins[i+1])
    a,b = np.histogram(prop_dict["duration"][mask_N][mask_zenith_bin], bins =dur_bins)
    hist_zenith_dur.append(a/c[i])

X,Y = np.meshgrid(zenith_bins, dur_bins)
im13 = ax13.pcolormesh(X,Y,np.array(hist_zenith_dur).T)
ax13.set_xlabel("Zenith [$^\circ$]")
ax13.set_ylabel("Duration [s]")
ax13.set_yscale("log")
yticks = np.array([5,10,20,50,100])
ax13.set_yticks(yticks)
yticks_label = yticks.astype("str")
ax13.set_yticklabels(yticks_label)
fig13.colorbar(im13, label = "Rate of Trails [counts/run]")
fig13.savefig(path_to_plots+"hist_2D_zenith_duration_rate.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")

fig14_0, ax14_0 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
bright_bins = np.logspace(np.log10(900), np.log10(np.max(prop_dict["avg_brightness"])), 20 )
ax14_0.set_xlabel("Brightness [MHz]")
ax14_0.set_ylabel("Number of Trails [counts]")
a,b = np.histogram(prop_dict["avg_brightness"][mask_N], bins = bright_bins)
bin_avg = (b[:-1]+b[1:])/2
ax14_0.errorbar(bin_avg,a,  yerr = np.sqrt(a), fmt = ".k")
for i in range(len(bin_avg)):
    ax14_0.hlines(y = a[i], xmin = b[i], xmax = b[i+1], color = "black")

ax14_0.set_xscale("log")
ax14_0.set_yscale("log")
xticks = np.array([1000, 2000, 5000, 10000])
yticks = np.array([1,2,5,0,20,50, 100, 200])
ax14_0.set_xticks(xticks)
ax14_0.set_yticks(yticks)
xticks_label = xticks.astype("str")
ax14_0.set_xticklabels(xticks_label)
yticks_label = yticks.astype("str")
ax14_0.set_yticklabels(yticks_label)
fig14_0.savefig(path_to_plots+"hist_brightness.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")



fig14, ax14 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
a,b,c = np.histogram2d(prop_dict["zenith"][mask_N], prop_dict["avg_brightness"][mask_N],
                       bins = [zenith_bins, bright_bins])
X,Y = np.meshgrid(zenith_bins, bright_bins)
im14 = ax14.pcolormesh(X,Y,a.T)
ax14.set_xlabel("Zenith [$^\circ$]")
ax14.set_ylabel("Brightness [MHz]")
ax14.set_yscale("log")
yticks = np.array([1000,2000, 3000,4000, 6000])
ax14.set_yticks(yticks)
yticks_label = yticks.astype("str")
old_label = ax14.get_yticklabels()
ax14.set_yticklabels([""]*len(old_label))
print(old_label) 
ax14.set_yticklabels(yticks_label)
#ax14.set_ylim((900, np.max(prop_dict["avg_brightness"])))
fig14.colorbar(im14, label = "Number of Trails [counts]")
fig14.savefig(path_to_plots+"hist_2D_zenith_brightness.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")


fig14_1, ax14_1 = plt.subplots(1, 1, figsize=figure_size, constrained_layout=True)
hist_zenith_bright = []
c,d = np.histogram(prop_dict["zenith"][index_first_run], bins = zenith_bins)
for i in range(len(zenith_bins)-1):
    mask_zenith_bin = np.logical_and(prop_dict["zenith"][mask_N]>=zenith_bins[i],
                            prop_dict["zenith"][mask_N]<zenith_bins[i+1])
    a,b = np.histogram(prop_dict["avg_brightness"][mask_N][mask_zenith_bin], bins =bright_bins)
    hist_zenith_bright.append(a/c[i])

X,Y = np.meshgrid(zenith_bins, bright_bins)
im14_1 = ax14_1.pcolormesh(X,Y,np.array(hist_zenith_bright).T)
ax14_1.set_xlabel("Zenith [$^\circ$]")
ax14_1.set_ylabel("Brightness [MHz]")
ax14_1.set_yscale("log")
yticks = np.array([1000,2000,3000,4000, 6000])
ax14_1.set_yticks(yticks)
yticks_label = yticks.astype("str")
ax14_1.set_yticklabels(yticks_label)
ax14_1.set_ylim((900, np.max(prop_dict["avg_brightness"])))
fig14_1.colorbar(im14_1, label = "Rate of Trails [counts/run]")
fig14_1.savefig(path_to_plots+"hist_2D_zenith_brightness_rate.jpg")
plt.close()

counter_of_plots+=1
print(counter_of_plots,"Plots created")






hfile.close()
