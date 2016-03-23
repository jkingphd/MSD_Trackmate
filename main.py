# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
import sys
from setup import *

## Comment out the following line to hardcode the filename
fname = sys.argv[1]
#fname = 'kpm2D15_Tracks.xml'

#%% 
## Set experimental parameters     
nMin = 50
eta = 0.000932
T = 273.15 + 23
scale = 0.28
dt = 1.0 # time interval in seconds
## These parameters should set to True if a drift velocity is present
## I recommend leaving y_drift as True for sedimentation experiments.
x_drift = False
y_drift = True

#%%
## Read datafile and scale x-, y-, and t- measurements
print 'Reading datafile and scaling measurements.'
tree = ET.parse(fname)
root = tree.getroot()
data = []
nTracks = np.int(root.attrib['nTracks'])
idx = 0
for j in range(nTracks):
    Track = j+1
    nSpots = np.int(root[j].attrib['nSpots'])
    for i in range(nSpots):
        time = np.float(root[j][i].attrib['t'])
        x = np.float(root[j][i].attrib['x'])
        y = np.float(root[j][i].attrib['y'])
        z = np.float(root[j][i].attrib['z'])
        data.append([x, y, z, time, Track])
        idx += 1
data = np.array(data)

data[:,3] *= dt
data[:,0:2] *= scale
data[:,2] *= 2*scale

#%%
## Step through individual tracks and perform msd analysis
## The length parameter is the length of the sliding window used to sample
## the run data. For length = 4, each run is split into N segments of length 4,
## with an overlap of 2.
## A length of 6 gives good results for 1 µm - 10 µm diameter particles
length = 6
particles = np.unique(data[:,-1])
f = open(fname.split('.')[0] + '_out.csv', 'w')
f.write('Track #, Length, Final_Length, x, y, t, dia_fc, dia_msd,\
        dia_msd_error, v_drift, density_fc, density_msd, density_msd_error\n')
for particle in particles:
    print 'Particle %d' % particle
    ## Find all positions for a given run.
    idx = np.where(data[:,-1] == particle)[0]
    dia_fc = data[idx, 2]
    x = data[idx, 0]
    y = data[idx, 1]
    t = data[idx, 3]
    n_orig = x.size
    ## Determine if particle is stuck or sticks during the run.
    ## If it sticks, truncate the run to the point in time where it sticks.
    a = stick_check(x, y, length = 6, threshold = 0.28)
    if a:
        x = x[:a[0]]
        y = y[:a[0]]
        t = t[:a[0]]
    x_0 = x[0]/scale
    y_0 = y[0]/scale
    t_0 = int(t[0]/dt)
    ## Check the length of the truncated/untruncated run. Only perford msd
    ## analysis on runs that equal or exceed nMin.
    n = x.size
    if n >= nMin:
        ## This is specifically for settling experiments. Calculate the
        ## average y-velocity over the entire run.
        v_drift = (y[-1] - y[0])/(t[-1] - t[0])
        ## Compensate for a steady drift velocity if flagged. This is very
        ## important. If drift velocity is present and not corrected for, msd
        ## analysis WILL fail.
        if x_drift == True:
            x = offset(x, t)
        if y_drift == True:
            y = offset(y, t)
        ## Split position and time arrays into segments and perform msd
        ## analysis on each segment, then calculate the mean result for the
        ## diameter and the standard error.
        xs = split_array(x, length, length/2)
        ys = split_array(y, length, length/2)
        ts = split_array(t, length, length/2)
        msd = []
        for i in range(xs.shape[0]):
            dx = (xs[i] - xs[i][0])[1:]
            dy = (ys[i] - ys[i][0])[1:]
            dT = (ts[i] - ts[i][0])[1:]
            dr2 = np.power(dx, 2) + np.power(dy, 2)
            msd.append(dr2/dT)
        msd = np.array(msd)
        means = []
        for i, row in enumerate(msd):
            means.append(np.mean(row))
        means = np.array(means)
        dia = calc_dia(np.mean(means/4.), eta, T)
        SE = StandardError(means/4.)
        dia_error = np.abs(calc_dia(np.mean(means/4.) + SE, eta, T) - dia)
        ## Calculate density from both the flowcam-reported diameter (the 
        ## minimum value of all reported diameters for the run) and the 
        ## msd-derived diameter
        density_msd = calc_density(dia, v_drift, eta)
        density_msd_error = np.abs(calc_density(dia + dia_error, v_drift, eta) - density_msd)
        dia_fc = np.min(dia_fc)
        density_fc = calc_density(dia_fc, v_drift, eta)
        f.write('%d, %d, %d, %0.3f, %0.3f, %d, %0.3f, %0.3f, %0.3f, %0.3f,\
                %0.3f, %0.3f, %0.3f\n' % (particle, n_orig, n, x_0, y_0,\
                t_0, dia_fc, dia, dia_error, v_drift, density_fc,\
                density_msd, density_msd_error))
    else:
        print 'Track too short.'
f.close()