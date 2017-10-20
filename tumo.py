#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:41 2017

@author: simonvanovereem
"""

#Imports
import math as m
import numpy as np

#Constants
rho_w   = 0.01      #[day^-1]
rho_g   = 0.005     #[day^-1]
rho_avg = 0.5*(rho_w + rho_g) #[day^-1]
D       = 0.48e-2   #[cm^2/day]
cd      = 4.e4      #[cells/cm^2]
R       = 10.       #[cm]
R_g     = 2.        #[cm]
RI      = 8.        #[cm]
phi     = m.pi/8.   #[rad]
phil    = m.pi/4.   #[rad]

n       = 10.       #[-]
dr      = R/n       #[cm]
rc      = np.linspace(0.+0.5*dr,R-0.5*dr,n) #[cm]
dtheta  = phil/n    #[cm]
theta   = np.linspace(0.+0.5*dtheta,phil-0.5*dtheta,n)

for i in range(0,int(n)):
    for j in range(0,int(n)):
        if abs(rc[i] - R_g) < 0.5*dr or (abs(rc[i] - RI) < 0.5*dr and abs(theta[j] - phi) < 0.5*dtheta):
            rho = rho_avg
        elif rc[i] < R_g or (rc[i] > RI and theta[j] < phi):
            rho = rho_g
        else:
            rho = rho_w
        # Internal nodes
        Cc = -(2.*dr/(rc[i]*dtheta) + 
               (rc[i+1]*dtheta)/dr + 
               (rc[i-1]*dtheta)/dr + 
               rc[i]*rho*dr*dtheta)
        Cn = dr/(rc[i]*dtheta)
        Ce = ((rc[i]+0.5*dr)*dtheta)/dr
        Cs = dr/(rc[i]*dtheta)
        Cw = ((rc[i]-0.5*dr)*dtheta)/dr
        