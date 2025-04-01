import ROOT
import os
import sys
import tables
import numpy as np
import time
import pandas as pd

starttime=time.time()

HESSDST = '/lfs/l7/hess/users/marandon/CalibData/CT5/dst_NewScheme2/'

runnumber = int(sys.argv[1])
outdir = str(sys.argv[2])

outfilename = outdir+'/'+str(runnumber)+'_intensity.pkl'

# rootlogon.C equivalent
if ROOT.gSystem.Getenv("HESSUSER"):
    print('Using HESSUSER')
    ROOT.gSystem.SetIncludePath("-I$HESSUSER/include -I$HESSROOT/include")
    ROOT.gROOT.SetMacroPath(".:$HESSUSER:$HESSROOT")
else:
    print('Using HESSROOT')
    ROOT.gSystem.SetIncludePath("-I$HESSROOT/include")
    ROOT.gROOT.SetMacroPath(".:$HESSROOT")

def calc_filepath(runnumber,dstname=HESSDST):
    if (runnumber-runnumber%200)<100000:
        lowbound = '0'+str(runnumber-runnumber%200)
        highbound = '0'+str(int(lowbound) + 199)
        runnumber = '0'+str(runnumber)
    else:
        lowbound = str(runnumber-runnumber%200)
        highbound = str(int(lowbound) + 199)
        runnumber=str(runnumber)

    filepath = HESSDST+'run'+lowbound+'-'+highbound+'/run_'+runnumber+'_DST_001.root'
    return filepath


f = ROOT.TFile(calc_filepath(runnumber), 'READ')
 
sash_data = f.Get("run")
sash_data.GetEntry(0)
 
events = f.Get("DST")
 
hess = sash_data.GetHESSArray()
 
event_header = hess.Get(ROOT.Sash.EventHeader())
run_header = hess.Get(ROOT.Sash.RunHeader())
 
tels_in_run = run_header.GetTelsInRun()
tels_iterator = tels_in_run.begin()
participating_tels = []

while tels_iterator != tels_in_run.end():
    participating_tels.append(tels_iterator.GetId())
    tels_iterator += 1

run_number = run_header.GetRunNum()
run_duration = run_header.GetDuration()
run_timestamp = run_header.GetTimeStamp().GetUTC().GetUnixTime()

numdurt0=[run_number,run_duration,run_timestamp]

dataset_iterator = events.begin()


eventno=[]
pixelids=[]
intensities=[]
twithdata=[]
timestamps=[]
hillaslengths=[]
hillaswidths=[]
hillasamps=[]
hillasalphas=[]
hillasskewness=[]
hillaskurtosis=[]
hillasphis=[]
hillaslds=[]
time_nanosecond=[]
timestamp=[]
iterator=0

while dataset_iterator != events.end():
    eventno.append(iterator)
    iterator+=1
    ev_ids=[]
    ev_intensities=[]
    ev_twithdata=[]
    ev_hillaslength=[]
    ev_hillaswidth=[]
    ev_hillasamp=[]
    ev_hillasalpha=[]
    ev_hillaskurtosis=[]
    ev_hillasskewness=[]
    ev_hillasphi=[]
    ev_hillascog=[]
    ev_hillasld=[]
    ev_nanosec=[]

    dataset_iterator.Process(0,0)
    tels_with_data = event_header.GetTelWData()

    time_nanosecond.append(event_header.GetTimeStamp().fNanosec)
    tels_iterator = tels_with_data.begin()

    while tels_iterator != tels_with_data.end():
        ev_twithdata.append(tels_iterator.GetId())
        hillas = hess.Get(tels_iterator.GetId(),"Hillas0407", ROOT.Reco.HillasParameters())

        hwidth = hillas.GetWidth()
        hlength = hillas.GetLength()
        hamp = hillas.GetImageAmplitude()
        halpha = hillas.GetAlpha().GetDegrees()
        hskewness = hillas.GetSkewness()
        hkurtosis = hillas.GetKurtosis()
        hphi = hillas.GetPhi().GetDegrees()
        hcog = hillas.GetCenterOfGravity()
        hld = hillas.GetLocalDistance()

        ev_hillaslength.append(hlength)
        ev_hillaswidth.append(hwidth)
        ev_hillasamp.append(hamp)
        ev_hillasalpha.append(halpha)
        ev_hillasskewness.append(hskewness)
        ev_hillaskurtosis.append(hkurtosis)
        ev_hillascog.append(hcog)
        ev_hillasphi.append(hphi)
        ev_hillasld.append(hld)

        intensity_data = hess.Get(tels_iterator.GetId(),"Extended0407",ROOT.Sash.IntensityData())
        pixel_list = tels_iterator.__deref__().__deref__().GetConfig().GetPixelsInTel()
        pixel_list = intensity_data.GetNonZero()
        pixel_iterator = pixel_list.begin()
        idlist=[]
        intensitylist=[]

        while pixel_iterator != pixel_list.end():
            pixel_id = pixel_iterator.__deref__().__deref__().GetConfig().GetPixelID()
            idlist.append(pixel_id)
            pixel_intensity = intensity_data.GetIntensity(pixel_iterator.__deref__()).fIntensity
            intensitylist.append(pixel_intensity)
            pixel_iterator += 1

        ev_ids.append(idlist)
        ev_intensities.append(intensitylist)
        
        tels_iterator += 1

    if len(ev_twithdata)>0:
        pixelids.append(ev_ids)
        intensities.append(ev_intensities)
        twithdata.append(ev_twithdata)
        hillasamps.append(ev_hillasamp)
        hillaswidths.append(ev_hillaswidth)
        hillaslengths.append(ev_hillaslength)
        hillasalphas.append(ev_hillasalpha)
        hillasskewness.append(ev_hillasskewness)
        hillaskurtosis.append(ev_hillaskurtosis)
        hillasphis.append(ev_hillasphi)
        hillaslds.append(ev_hillasld)

'''
def writefile(h5file,arrname,arr):

    arr=np.asarray(arr,dtype=object)
    myobjects = h5file.create_vlarray(h5file.root, arrname, tables.ObjectAtom())
    myobjects.append(arr)

    return 0

h5file = tables.open_file(outfilename, 'w')

writefile(h5file,'numdurt0',numdurt0)
writefile(h5file,'pixelids',pixelids)
writefile(h5file,'intensities',intensities)
writefile(h5file,'twithdata',twithdata)
writefile(h5file,'hillasamps',hillasamps)
writefile(h5file,'hillaswidths',hillaswidths)
writefile(h5file,'hillaslengths',hillaslengths)
writefile(h5file,'hillasalphas',hillasalphas)
writefile(h5file,'hillasskewness',hillasskewness)
writefile(h5file,'hillaskurtosis',hillaskurtosis)
writefile(h5file,'hillasphis',hillasphis)
writefile(h5file,'time_nanosecond',time_nanosecond)

h5file.close()
'''

mydf=pd.DataFrame({'numdurt0':numdurt0,'pixelids':pixelids,'intensities':intensities,'twithdata':twithdata,'hillasamps':hillasamps,'hillasalphas':hillasalphas,'hillasskewness':hillasskewness,'hillaskurtosis':hillaskurtosis,'hillasphis':hillasphis,'time_nanosecond':time_nanosecond})

mydf.to_pickle(outfilename)
endtime=time.time()
print('Time Elapsed:', endtime-starttime)

f.Close()
