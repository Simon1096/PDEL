#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:41 2017

@author: simonvanovereem
"""

#Imports
import math as m
import numpy as np
np.set_printoptions(threshold=np.nan)

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

n       = 3.       #[-]
dr      = R/n       #[cm]
rc      = np.linspace(0.+0.5*dr,R-0.5*dr,n) #[cm]
dtheta  = phil/n    #[cm]
theta   = np.linspace(0.+0.5*dtheta,phil-0.5*dtheta,n)
A       = np.zeros((int(n**2),int(n**2)))

def left():
    if j == 0:
        Cc = -(dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr +  
               rc[i]*rho*dr*dtheta)
        Cn = dr/(rc[i]*dtheta)
        Ce = ((rc[i+1]*0.5*dr)*dtheta)/dr
        Cs = 0
    elif j == int(n) - 1:
        Cc = -(dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr +  
               rc[i]*rho*dr*dtheta)
        Cn = 0
        Ce = ((rc[i+1]*0.5*dr)*dtheta)/dr
        Cs = dr/(rc[i-int(n)]*dtheta)
    else:
        Cc = -(2.*dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr +  
               rc[i]*rho*dr*dtheta)
        Cn = dr/(rc[i]*dtheta)
        Ce = ((rc[i+1]*0.5*dr)*dtheta)/dr
        Cs = dr/(rc[i]*dtheta)
    return Cc, Cn, Ce, Cs

def lower():
    if i == int(n) - 1:
        Cc = -(dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr + 
               rc[i]*rho*dr*dtheta)
        Cn = dr/(rc[i]*dtheta)
        Ce = 0
        Cw = ((rc[i-1]*0.5*dr)*dtheta)/dr
    else:
        Cc = -(2.*dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr + 
               (0.5*dr*dtheta)/dr + 
               rc[i]*rho*dr*dtheta)
        Cn = dr/(rc[i]*dtheta)
        Ce = ((rc[i+1]*0.5*dr)*dtheta)/dr
        Cw = ((rc[i-1]*0.5*dr)*dtheta)/dr
    return Cc, Cn, Ce, Cw

def right():
    if j == int(n) - 1:
        Cc = -(dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr + 
               rc[i]*rho*dr*dtheta)
        Cn = 0
        Cs = dr/(rc[i]*dtheta)
        Cw = ((rc[i-1]*0.5*dr)*dtheta)/dr
    else:
        Cc = -(2.*dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr + 
               rc[i]*rho*dr*dtheta)
        Cn = dr/(rc[i]*dtheta)
        Cs = dr/(rc[i]*dtheta)
        Cw = ((rc[i-1]*0.5*dr)*dtheta)/dr
    return Cc, Cn, Cs, Cw

def upper():
    if i != 0 and i != int(n) - 1:
        Cc = -(dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr + 
               (0.5*dr*dtheta)/dr + 
               rc[i]*rho*dr*dtheta)
        Ce = ((rc[i+1]*0.5*dr)*dtheta)/dr
        Cs = dr/(rc[i-int(n)]*dtheta)
        Cw = ((rc[i-1]*0.5*dr)*dtheta)/dr
    return Cc, Ce, Cs, Cw

def internal():
    Cc = -(2.*dr/(rc[i]*dtheta) + 
               (0.5*dr*dtheta)/dr + 
               (0.5*dr*dtheta)/dr + 
               rc[i]*rho*dr*dtheta)
    Cn = dr/(rc[i]*dtheta)
    Ce = ((rc[i+1]*0.5*dr)*dtheta)/dr
    Cs = dr/(rc[i-int(n)]*dtheta)
    Cw = ((rc[i-1]*0.5*dr)*dtheta)/dr
    return Cc, Cn, Ce, Cs, Cw
    

for j in range(0,int(n)):
    for i in range(0,int(n)):
        if abs(rc[i] - R_g) < 0.5*dr or (abs(rc[i] - RI) < 0.5*dr and abs(theta[j] - phi) < 0.5*dtheta):
            rho = rho_avg
        elif rc[i] < R_g or (rc[int(i%n)] > RI and theta[int(j%n)] < phi):
            rho = rho_g
        else:
            rho = rho_w
            
        if i == 0:
            (Cc, Cn, Ce, Cs) = left()
            A[int(j*n)][int(j*n)] = Cc

            if (int(j*n+n) < n**2):
              A[int(j*n)][int(j*n+n)] = Cn

            if (int(j*n-n) >= 0):
              A[int(j*n)][int(j*n-n)] = Cs

            A[int(j*n)][int(j*n+1)] = Ce
            
        elif j == 0 and i > 0:
            (Cc, Cn, Ce, Cw) = lower()

            A[int(j*n+i)][int(j*n+i)] = Cc

            if (int(j*n+i+n) < n**2):
              A[int(j*n+i)][int(j*n+i+n)] = Cn

            if (int(j*n+i+1) < n**2):
              A[int(j*n+i)][int(j*n+i+1)] = Ce

            A[int(j*n+i)][int(j*n+i-1)] = Cw


        elif i == int(n) - 1 and j != 0:
            (Cc, Cn, Cs, Cw) = right()
            A[int(j*n+i)][int(j*n+i)] = Cc

            if (int(j*n+i+n) < n**2):
              A[int(j*n+i)][int(j*n+i+n)] = Cn

            if (int(j*n-n+i) >= 0):
              A[int(j*n+i)][int(j*n-n+i)] = Cs

            A[int(j*n+i)][int(j*n+i-1)] = Cw
            
        elif j == int(n) - 1 and i != 0 and i != int(n) - 1:
            (Cc, Ce, Cs, Cw) = upper()
            A[int(j*n+i)][int(j*n+i)] = Cc
            A[int(j*n+i)][int(j*n+i-n)] = Cs
            A[int(j*n+i)][int(j*n+i+1)] = Ce
            A[int(j*n+i)][int(j*n+i-1)] = Cw
        
        else:
            (Cc, Cn, Ce, Cs, Cw) = internal()
            A[int(j*n+i)][int(j*n+i)] = Cc
            A[int(j*n+i)][int(j*n+i+n)] = Cn
            A[int(j*n+i)][int(j*n+i-n)] = Cs
            A[int(j*n+i)][int(j*n+i+1)] = Ce
            A[int(j*n+i)][int(j*n+i-1)] = Cw

# TODO: 
# - Divide everything in the matrix by Rc
# - Choose time integration method
