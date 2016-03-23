# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats
from scipy import constants

def split_array(x, length, stagger):
    """
    split array x in chunks of length with stagger (amount of overlap)
    """
    L = [x[i:i+length] for i in range(0,len(x),stagger)]
    L = [xx for xx in L if len(xx)==length]
    return np.array(L)

def stick_check(x, y, length = 6, threshold = 0.1):
    xs = split_array(x, length, 1)
    ys = split_array(y, length, 1)
    idx = []
    for i in range(len(xs)):
        dx = xs[i] - xs[i][0]
        dy = ys[i] - ys[i][0]
        if np.all(np.abs(dx) < threshold) and np.all(np.abs(dy) < threshold):
            idx.append(i)
    return idx

def calc_dia(D, eta = 0.00089, T = 298.15):
    '''
    Calculate hydrodynamic diameter from Stokes-Einstein equation. 
    Needs diffusion coefficient, dynamic viscosity [N*s/m^2], and temperature [K].
    '''
    D *= 1E-12 # m^2/s
    r = constants.k*T/(6*constants.pi*eta*D)
    result = 2*r*1E6 # microns
    return result

def offset(y, t):
    v_avg = (y[-1] - y[0])/(t[-1] - t[0])
    result = y - v_avg*t
    return result

def StandardError(x, conf=0.95, ddof=1):
    """
    calculate the standard error over array

    *Parameters*

    x: array
     array of values to consider

    conf: float
     confidance interval to apply [default=0.95]

    ddof: int
     degree of freedom correction.  dof=size(x)-ddof

    *Returns*

    SE: float
     standard error.  Calculated with students t distribution
    """

    _n=np.size(x)
    if(_n<=1):
        _SE=0
    else:
        #ddof=0, standard dev, ddof=1, sample standard dev
        _SE=np.std(x,ddof=1)/np.sqrt(_n)*scipy.stats.t.isf(0.5*(1-conf),_n-1)
    return _SE
    
def calc_density(d, v_drift, eta = 0.00089):
    '''
    Returns delta rho.
    Needs a hydrodynamic diameter [µm], average velocity [µm/s], and dynamic viscosity [N*s/m^2].
    '''
    r = 0.5*d
    r *= 1E-6 # m
    v_drift *= 1E-6 # m/s
    result = (9*eta*v_drift)/(2*constants.g*r**2)
    return result