import os
home_path = os.path.expanduser("~/")
directory_hillas_params = "/lfs/l7/hess/users/spencers/realdata/Hillas0510_alldata/"
path_to_save_files = os.path.expanduser("~/")  +"hillas_parameter"
exec(open("functions_import_essentials.py").read())
import tables as tb

def Flatten(arr_list):
    a = [item for sublist in arr_list for item in sublist]
    a = np.array(a)
    return a

def read_hillas_file(nrun):
    print("Reading file "+str(nrun)+"_params.hdf5")
    try:
        x = tb.open_file(directory_hillas_params+str(nrun)+"_params.hdf5", "r")
    except:
        print("File "+str(nrun)+"_params.hdf5 is not readable")
        x = -1
    return x

nrun = int(sys.argv[1])
hfile = h5py.File(home_path+"trail_properties/run_properties.h5", "r")
props = {}
idx = np.where(np.array(hfile["nrun"])==nrun)[0]
for key in list(hfile.keys()):
    props[key]=np.array(hfile[key])[idx]

hfile.close()
start_time = props["start"]-props["t0"]
stop_time = props["duration"]+props["start"]-props["t0"]
track = props["N_track"]
print(start_time,stop_time)#, brightness, len(track))

x = read_hillas_file(nrun)
if x == -1:
    print("Run", nrun, "has no trails")
    exit()
    #

time = (Time(np.array(x.root.time_second[:][0],dtype='double'),
             format='unix')+ x.root.time_nanosecond[:] * u.ns).to_value(format = "unix")
time = time[0] - x.root.numdurt0[:][0][2] #numdurt0 gives time of run start in unix
twithdata = x.root.twithdata[:][0]

hillas_param_names = ['amps', 'lengths', 'widths', 'kurtosis', 'skewness', 'phis', 'alphas']
dict_hillas_parameter = {}
for N in range(len(track)):
    dict_hillas_parameter[N] = {}
    for str_hillas in hillas_param_names:
        dict_hillas_parameter[N][str_hillas] = {}
        for telId in np.arange(1,6):
            dict_hillas_parameter[N][str_hillas][telId] = {}
            # will later contain keys "during" and "around" with complete array of hillas params
            
for N in range(len(track)):
    if track[0] == -1:
        print("No tracks detected")
        break
    duration = stop_time[N]-start_time[N]
    idx_during = np.where(np.logical_and(time >= start_time[N],
                                         time <= stop_time[N]))[0]
    idx_before = np.where(np.logical_and(time >= start_time[N]-(60+duration/2),
                                         time <= start_time[N]-60))[0]
    idx_after = np.where(np.logical_and(time >= stop_time[N]+60,
                                         time <= stop_time[N]+60+duration/2))[0]
    all_idx = np.append(np.append(idx_before, idx_during), idx_after)
    
    #all_idx = np.where(np.logical_and(time >= start_time[N]-60,
    #                                     time <= stop_time[N]+60))[0]
    mask_before = all_idx <= np.max(idx_before)
    mask_after = all_idx >= np.min(idx_after)
    mask_during = np.logical_and(all_idx >=np.min(idx_during),
                                 all_idx <= np.max(idx_during))
    twd = twithdata[all_idx]
    params = {}
    params["amps"]     = x.root.hillasamps[:][0][all_idx]
    params["lengths"]  = x.root.hillaslengths[:][0][all_idx]
    params["widths"]   = x.root.hillaswidths[:][0][all_idx]
    params["kurtosis"] = x.root.hillaskurtosis[:][0][all_idx]
    params["skewness"] = x.root.hillasskewness[:][0][all_idx]
    #Phis and Alphas seem to be not dependent on track, leave therefore away for faster evaluation?
    params["phis"]     = x.root.hillasphis[:][0][all_idx]
    params["alphas"]   = x.root.hillasalphas[:][0][all_idx]
    paramnames = []

    ##########
    #Build a mask here looking for len of each list, look only at the ones with more than one entry
    #mask_hybrid1 = np.array([len(twd[i])>1 for i in range(len(twd))])
    #mask_hybrid2 = np.array([twd[mask_hybrid1][i] for i in twd[mask_hybrid1] if 5 in i)
    ##########
    for str_hillas in list(params.keys()):#Used to keep everything working if phis and alphas are excluded
        paramnames.append(str_hillas)
    twd_flat = {}
    twd_flat["before"] = Flatten(twd[mask_before])
    twd_flat["during"] = Flatten(twd[mask_during])  # Used to assign correct telId to event
    twd_flat["after"]  = Flatten(twd[mask_after])
    str_hillas = "lengths"
    mask_len_before = Flatten(params[str_hillas][mask_before])>0
    mask_len_during = Flatten(params[str_hillas][mask_during])>0 #Making a length cut removes negative amplitudes as well
    mask_len_after  = Flatten(params[str_hillas][mask_after])>0
    for str_hillas in paramnames:
        params[str_hillas+" before"] = Flatten(params[str_hillas][mask_before])[mask_len_before]
        params[str_hillas+" during"] = Flatten(params[str_hillas][mask_during])[mask_len_during]
        params[str_hillas+" after"]  = Flatten(params[str_hillas][mask_after])[mask_len_after]
    masks = {}
    for telId in np.arange(1,6):
        masks[telId] = {}
        masks[telId]["before"] = twd_flat["before"][mask_len_before] == telId
        masks[telId]["during"] = twd_flat["during"][mask_len_during] == telId
        masks[telId]["after"]  = twd_flat["after"][mask_len_after] == telId
    #
    for str_hillas in paramnames:
        for telId in np.arange(1,6):
            dict_hillas_parameter[N][str_hillas][telId]["during"] = params[str_hillas+" during"][masks[telId]["during"]]
            dict_hillas_parameter[N][str_hillas][telId]["around"] = np.append(params[str_hillas+" before"][masks[telId]["before"]],
                                                                              params[str_hillas+" after"][masks[telId]["after"]])
    local_time = time[all_idx]
    for telId in np.arange(1,6):
        telId_indices = [j for j in range(len(twd)) if telId in twd[j]]
        hist, bin_edges = np.histogram(local_time[telId_indices] , 
                                       bins = np.arange(int(np.min(local_time[telId_indices])) ,int(np.max(local_time[telId_indices])), 1))
        dict_hillas_parameter[N]["event rate"] = {}
        dict_hillas_parameter[N]["event rate"][telId] = {}
        dict_hillas_parameter[N]["event rate"][telId]["hist"] = hist
        dict_hillas_parameter[N]["event rate"][telId]["bin_edges"] = bin_edges
        dict_hillas_parameter[N]["event rate"][telId]["raw data"] = time[all_idx][telId_indices]

#with open(path_to_save_files+"/dict_hillas_params"+str(nrun)+".txt", "wb") as handle:
#    pickle.dump(dict_hillas_parameter, handle)
'''#still work in progress
hfile = h5py.File("testfile.h5", "w")
dict_group = hfile.create_group("N_trail")
for N,param in dict_hillas_parameter.items():
    subgroup = dict_group.create_group(param)
    for param, telId in dict_hillas_parameter[N].items():
        subsubgroup = subgroup.create_group(telId)
        for telId, pos in dict_hillas_parameter[N][params].items():
            subsubsubgroup = subsubgroup.create_group(pos)
            for pos, val in dict_hillas_parameter[N][params][telId].items():
                subsubsubgroup[pos] = val
 
'''
