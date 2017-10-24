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
Un      = np.array(cd*np.exp(-(rc)**2))
U       = np.zeros(int(n**2))

def internal():
  if i == 0:
    Cc = -(dr/(rc[i]*dtheta) + 
      (0.5*dr*dtheta)/dr +  
      rc[i]*rho*dr*dtheta)
  elif i == n - 1:
    Cc = -(dr/(rc[i]*dtheta) + 
      (0.5*dr*dtheta)/dr + 
      rc[i]*rho*dr*dtheta)
  elif j == 0:
    Cc = -(dr/(rc[i]*dtheta) + 
      (0.5*dr*dtheta)/dr + 
      rc[i]*rho*dr*dtheta)
  elif j == n - 1:
    Cc = -(dr/(rc[i]*dtheta) + 
      (0.5*dr*dtheta)/dr + 
      (0.5*dr*dtheta)/dr + 
      rc[i]*rho*dr*dtheta)
  else:
    Cc = -(2.*dr/(rc[i]*dtheta) + 
      (0.5*dr*dtheta)/dr + 
      (0.5*dr*dtheta)/dr + 
      rc[i]*rho*dr*dtheta)

  Cn = dr/(rc[i]*dtheta) if j < n else 0
  Cs = dr/(rc[i-int(n)]*dtheta) if j != 0 else 0
  Ce = ((rc[i+1]*0.5*dr)*dtheta)/dr if i + 1 < n else 0
  Cw = ((rc[i-1]*0.5*dr)*dtheta)/dr if i != 0 else 0
  return Cc/(rc[i]*dr*dtheta), Cn/(rc[i]*dr*dtheta), Ce/(rc[i]*dr*dtheta), Cs/(rc[i]*dr*dtheta), Cw/(rc[i]*dr*dtheta)
  
for j in range(0,int(n)):
    for i in range(0,int(n)):
        if abs(rc[i] - R_g) < 0.5*dr or (abs(rc[i] - RI) < 0.5*dr and abs(theta[j] - phi) < 0.5*dtheta):
            rho = rho_avg
        elif rc[i] < R_g or (rc[int(i%n)] > RI and theta[int(j%n)] < phi):
            rho = rho_g
        else:
            rho = rho_w
            
        (Cc, Cn, Ce, Cs, Cw) = internal()
        A[int(j*n+i)][int(j*n+i)] = Cc

        if j != n - 1:
          A[int(j*n+i)][int(j*n+i+n)] = Cn
        if j != 0:
          A[int(j*n+i)][int(j*n+i-n)] = Cs
        if i != n - 1:
          A[int(j*n+i)][int(j*n+i+1)] = Ce
        if i != 0:
          A[int(j*n+i)][int(j*n+i-1)] = Cw
            
    np.put(U,np.arange(int(j*n),int((j+1)*n)), Un)

def time():
    global U
    I   = np.identity(int(n**2))
    dt  =  0.01
    amp_mat = np.linalg.inv(I-dt*A)
    for i in range(0,1357):
        Unew = np.dot(amp_mat,U)
        U = Unew
    print(sum(i > cd for i in U))
    return U



