#!/usr/bin/env python
import os
import sys
import obspy
import getopt
import numpy as np
from obspy.taup import TauPyModel
from scipy import interpolate
from obspy.signal.trigger import recSTALTA
from obspy.signal.trigger import plotTrigger
import seispy
#import subprocess


try:
    opts,args = getopt.getopt(sys.argv[1:], "mI:C:F:B:O:vt")
except:
    print('Arguments are not found!')
#    sys.exit(1)

mccc = "./mccc"
dt = 0.025
ctime = [-10, 180]
verbose = ''
ismccc = 0
istaup=0
for op, value in opts:
    if op == "-I":
        indir = value
    elif op == "-C":
        cut_info = value
    elif op == "-F":
        filt_info = value
    elif op == "-B":
        bp_f1 = float(value.split('/')[0])
        bp_f2 = float(value.split('/')[1])
    elif op == "-O":
        outpath = value
    elif op == "-v":
        verbose = '-v'
    elif op == "-m":
        ismccc = 1
    elif op == "-t":
        istaup = 1
    else:
        sys.exit(1)

'''
delta_range = [30, 90]
model = TauPyModel(model="iasp91")
arrivals = model.get_travel_times(source_depth_in_km=600,
        distance_in_degree=delta_range[0],
        phase_list=["P"])
t1 = arrivals[0].time - 10
arrivals = model.get_travel_times(source_depth_in_km=10,
        distance_in_degree=delta_range[1],
        phase_list=["P"])
t2 = arrivals[0].time + 200
cut_info = '-3/'+str(t1)+'/'+str(t2)
'''
if istaup:
    print("Write P-wave traveltime into sac file...")
    os.system("taup_setsac -mod iasp91 -ph P-0 -evdpkm %s" % (os.path.join(indir,"*.SAC")))
if ismccc:
    print(mccc+" "+os.path.join(indir,"*.SAC")+" -D -C"+cut_info+" -F"+filt_info+" "+verbose)
    os.system(mccc+" "+os.path.join(indir,"*.SAC")+" -D -C"+cut_info+" -F"+filt_info+" "+verbose)
    
with open('tdel.dat','r') as fid:
    cc_info = [[cc.split('\t')[0], float(cc.split('\t')[1]), float(cc.split('\t')[2])] for cc in fid.readlines()]

data_sum = np.zeros(np.floor((ctime[1]-ctime[0])/dt))
st = obspy.core.Stream()
for info in cc_info:
    tr = obspy.read(info[0])[0]
    tr.detrend('linear')
    tr.detrend('constant')
    snrbegin = np.floor((tr.stats.sac.t0-20)/dt)
    snrend = np.floor((tr.stats.sac.t0+20)/dt)
    snro = np.floor(tr.stats.sac.t0/dt)
    snr = seispy.geo.snr(tr.data[snro:snrend], tr.data[snrbegin:snro])
    if snr < 15 or info[2]<0.6:
        continue
    nt1 = np.floor((tr.stats.sac.t0-info[1]+ctime[0])/dt)
    nt2 = np.floor((tr.stats.sac.t0-info[1]+ctime[1])/dt)
#    nt1 = np.floor((tr.stats.sac.t0+ctime[0])/dt)
#    nt2 = np.floor((tr.stats.sac.t0+ctime[1])/dt)

#    tr1 = tr.copy()
#    tr1.filter('bandpass',freqmin=0.05,freqmax=2)
#    data_sum += tr1.data[nt1:nt2]
    tr.filter('bandpass',freqmin=bp_f1, freqmax=bp_f2)
    newtr = tr.data[nt1:nt2]    
    tmp = obspy.Trace()
    tmp.data = newtr
    tmp.normalize()
    tmp.stats.delta = dt
    tmp.stats.network = tr.stats.network
    tmp.stats.station = tr.stats.station
    tmp.stats.channel = 'Z'
    tmp.stats.update(adict={'sac':{'evla':tr.stats.sac.evla,
                                  'evlo':tr.stats.sac.evlo,
                                  'evdp':tr.stats.sac.evdp,
                                  'stla':tr.stats.sac.stla,
                                  'stlo':tr.stats.sac.stlo,
                                  'b':ctime[0],
                                  'gcarc':tr.stats.sac.gcarc}})
    sacname = tr.stats.network+'.'+tr.stats.station+'.z'
    tmp.write(os.path.join(outpath,sacname),format='SAC')

#    st += tmp
    
'''
cft = recSTALTA(data_sum, int(0.5 / dt), int(10 / dt))

trigger = np.where(cft>1.5)[0][0]*dt
tmp = obspy.Trace()
tmp.data = data_sum
tmp.stats.delta = dt
plotTrigger(tmp, cft, 1.5, 0.5)
for tr in st:
#    cef =np.sum(template*tr.data)
#    print(tr.stats.network,tr.stats.station,cef)
#        tr.data = -tr.data
    tr.stats.sac.b = -trigger
    sacname = tr.stats.network+'.'+tr.stats.station+'.z'
    tr.write(os.path.join(outpath,sacname),format='SAC')
'''
