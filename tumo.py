#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:41 2017

@author: simonvanovereem
@minion: jasonxie
"""

#Imports
import math as m
import numpy as np
np.set_printoptions(threshold=np.nan)

#Constants
rho_w	= 0.01		#[day^-1]
rho_g	= 0.005	 #[day^-1]
rho_avg = 0.5*(rho_w + rho_g) #[day^-1]
D		= 0.48e-2	 #[cm^2/day]
cd		= 4.e4		#[cells/cm^2]
R		= 10.		 #[cm]
R_g		= 2.		#[cm]
RI		= 8.		#[cm]
phi		= m.pi/8.	 #[rad]
phil	= m.pi/4.	 #[rad]

N 		= [10]
dr 		= 0
dtheta 	= 0
rc 		= []
theta 	= []

def internal(i, j, n):
	if i == 0:
		# Left side
		if j == 0 or j == n - 1:
			Cc = -D * (dr / (rc[i]*dtheta) + 
				(rc[i]+0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta
		else:
			Cc = -D * (2.*dr / (rc[i]*dtheta) + 
				(rc[i]+0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta
	elif i == n - 1:
		# Right side
		if j == n - 1:
			Cc = -D * (dr / (rc[i]*dtheta) + 
				(rc[i]-0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta
		else:
			Cc = -D * (2.*dr / (rc[i]*dtheta) + 
				(rc[i]-0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta
	elif j == 0:
		# Bottom side
		if i == n - 1:
			Cc = -D * (dr / (rc[i]*dtheta) + 
				(rc[i]-0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta
		else:
			Cc = -D * (dr/(rc[i]*dtheta) + 
				(rc[i]+0.5*dr) * dtheta / dr + 
				(rc[i]-0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta
	elif j == n - 1:
		# Top side
		Cc = -D * (dr / (rc[i]*dtheta) + 
			(rc[i]+0.5*dr) * dtheta / dr + 
			(rc[i]-0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta
	else:
		# Internal nodes
		Cc = -D * (2.*dr / (rc[i]*dtheta) + 
			(rc[i]+0.5*dr) * dtheta / dr + 
			(rc[i]-0.5*dr) * dtheta / dr) + rc[i]*rho*dr*dtheta

	# R_x = 0.5 * dr
	Cn = D*dr / (rc[i]*dtheta) if j < n - 1 else 0
	Cs = D*dr / (rc[i]*dtheta) if j != 0 else 0
	Ce = D * ((rc[i]+0.5*dr)*dtheta) / dr if i < n-1 else 0
	Cw = D * ((rc[i]-0.5*dr)*dtheta) / dr if i != 0 else 0
	return Cc/(rc[i]*dr*dtheta), Cn/(rc[i]*dr*dtheta), Ce/(rc[i]*dr*dtheta), Cs/(rc[i]*dr*dtheta), Cw/(rc[i]*dr*dtheta)

def initMatrix(n):
	global rho, rc, theta, dtheta, dr

	dr		= R / n
	rc		= np.linspace(0.+0.5*dr,R-0.5*dr,n)
	dtheta	= phil / n
	theta	= np.linspace(0.+0.5*dtheta,phil-0.5*dtheta,n)

	A		= np.zeros((n**2,n**2))
	Un		= np.array(cd*np.exp(-(rc)**2))
	U		= np.zeros(n**2)
	for j in range(0,n):
		for i in range(0,n):
			if abs(rc[i] - R_g) < 0.5*dr or (abs(rc[i] - RI) < 0.5*dr and abs(theta[j] - phi) < 0.5*dtheta):
				rho = rho_avg
			elif rc[i] < R_g or (rc[i] > RI and theta[j] < phi):
				rho = rho_g
			else:
				rho = rho_w
				
			Cc, Cn, Ce, Cs, Cw = internal(i, j, n)
			A[j*n+i][j*n+i] = Cc

			if j != n - 1:
				A[j*n+i][j*n+i+n] = Cn

			if j != 0:
				A[j*n+i][j*n+i-n] = Cs

			if i != n - 1:
				A[j*n+i][j*n+i+1] = Ce

			if i != 0:
				A[j*n+i][j*n+i-1] = Cw
				
		np.put(U,np.arange(j*n,(j+1)*n), Un)
	return A, U

def time(A, U, n):
	I = np.identity(n**2)
	T = [1, 0.5, 0.1]
	for dt in T:
		Ut = U
		amp_mat = np.linalg.inv(I-dt*A)
		iterations = 0;
		while sum(c > cd for c in Ut) < 0.2 * len(U):
			Unew = np.dot(amp_mat,Ut)
			Ut = Unew
			iterations += 1
		print("deltaT: %.3f\nNumber of steps taken: %d" % (dt, iterations))
		tod = iterations * dt
		todStr = ""
		if tod > 365.25:
			count = int(tod / 365.25)
			todStr += "%d year" % count
			tod -= count * 365.25
			if count > 1:
				todStr += "s"
		if tod > 7:
			count = int(tod / 7)
			todStr += " %d week" % count
			tod -= count * 7
			if count > 1:
				todStr += "s"
		if tod > 0:
			todStr += " %.3f days" % tod
			if tod > 1:
				todStr += "s"

		print("Time till death: " + todStr)
	print("======================================")
	return U

print("======================================")
for n in N:
	A, U = initMatrix(n)
	print("Number of nodes in each direction: %d" % n)
	time(A, U, n)