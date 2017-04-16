#!/usr/bin/env python

# Simulation of a chaotic three-varialble Belousov-Zhabotinsky reaction in a CSTR
# Based on the model from Gyorgyi & Field (1992)
# By Steven Selverston
# For: Prof. Jim Crutchfield
# PHY 150/250
# 6/9/2009


import scipy as sp
import numpy as np
from scipy import integrate
import pylab as pl
from enthought.mayavi import mlab
from enthought.mayavi.mlab import axes,plot3d,points3d


print 'Welcome to the Chaotic Belousov-Zhabotinsky Simulator! \n '
print 'Press the "Return" button to choose the Default values. \n '

plotrange = raw_input('Enter the Time Range (Dimensionless):          (Default = 5)' '\n')
if plotrange == '': plotrange = 5
else: plotrange = float(plotrange)


	
# All the parameter values are from Gyorgyi & Field 1992

# Get user input to build a custom simulation

rangekf = raw_input('High kf range (1) or Low kf range (2)?:        (Default = High kf)' '\n')
if rangekf == '': rangekf = 1
else: rangekf = int(rangekf)

plotprofile = raw_input('Plot the 2D Profile? Yes (1)  No (2):          (Default = Yes)' '\n')
if plotprofile == '': plotprofile = 1
else: plotprofile = int(plotprofile)

# This gets the user input on whether or not to plot the attractor as well
plotattractor = raw_input('Plot the 3D Attractor? Yes (1)  No (2):        (Default = No)' '\n')
if plotattractor == '': plotattractor = 2
else: plotattractor = int(plotattractor)



if rangekf == 1:
	
	stepsize = raw_input('Enter the Step Size:                           (Default = 0.01)' '\n')
	if stepsize == '': stepsize = 0.01
	else: stepsize = float(stepsize)
	
	kf = raw_input('Enter the kf:                                  (Default = 0.00216)' '\n')
	if kf == '': kf= 0.00216
	else: kf = float(kf)
	
	A = 0.14   # BrO3-
	M = 0.3  # MA
	H = 0.26   # H+
	C = 0.001  # Ce-total    high kf
	alpha = 333.3   # high kf
	beta = 0.2609
	uo = np.array([.446751,5.275282,.393890])   #bz2 high kf
	
# This next if statement just assigns different parameters and initial guesses (for Low flowrate conditions)
if rangekf == 2:
	stepsize = raw_input('Enter the Step Size:  (Default = 0.005)' '\n')
	if stepsize == '': stepsize = 0.005
	else: stepsize = float(stepsize)
	
	kf = raw_input('Enter the kf:  (Default = 0.0003902)' '\n')
	if kf == '': kf = 0.0003902
	else: kf = float(kf)
	A = 0.1   # BrO3-
	M = 0.25  # MA
	H = 0.26   # H+
	C = 0.000833  # Ce-total    high kf
	alpha = 666.7   # high kf
	beta = 0.3478
	uo = np.array([.0468627,.89870,.846515])  # bz2 low kf

# let the user know that the program is running
print 'running...'

iters = np.arange(0.0,plotrange,stepsize)	

k1 = 4000000
k2 = 2.0
k3 = 3000
k4 = 55.2
k5 = 7000
k6 = 0.09
k7 = 0.23

To = (10*k2*A*H*C)**(-1)
Vo = 4*A*H*C/M**2
Xo = k2*A*H**2/k5
Yo = 4*k2*A*H**2/k5
Zo = C*A/40*M
	
# The equations are from Gyorgyi & Field 1992, although they are in a slightly different form

# the definition of BZ writes the equations in a vector form for use in SciPy's integrate module odeint

# Here x represents dx/dt, y represents dy/dt, and z represents dz/dt.

def BZ(u,t):
	x = To*(-k1*H*Yo*u[0]*((alpha*k6*Zo*Vo*u[1]*u[2]/(k1*H*Xo*u[0] + k2*A*H**2 + kf))/Yo) + k2*A*H**2*Yo/Xo*((alpha*k6*Zo*Vo*u[1]*u[2]/(k1*H*Xo*u[0] + k2*A*H**2 + kf))/Yo) - 2*k3*Xo*u[0]**2 + 0.5*k4*A**0.5*H**1.5*Xo**(-0.5)*(C-Zo*u[1])*u[0]**(0.5) - 0.5*k5*Zo*u[0]*u[1] - kf*u[0])
	
	y = To*(k4*A**0.5*H**1.5*Xo**(0.5)*(C/Zo-u[1])*u[0]**0.5 - k5*Xo*u[0]*u[1] - alpha*k6*Vo*u[1]*u[2] - beta*k7*M*u[1] - kf*u[1])
	
	z = To*(2*k1*H*Xo*Yo/Vo*u[0]*((alpha*k6*Zo*Vo*u[1]*u[2]/(k1*H*Xo*u[0] + k2*A*H**2 + kf))/Yo) + k2*A*H**2*Yo/Vo*((alpha*k6*Zo*Vo*u[1]*u[2]/(k1*H*Xo*u[0] + k2*A*H**2 + kf))/Yo) + k3*Xo**2/Vo*u[0]**2 - alpha*k6*Zo*u[1]*u[2] - kf*u[2])
	
	return np.array([x,y,z], dtype = float)


# The equations are integrated in vector form using odeint

results = integrate.odeint(BZ,uo,iters)

# The 2D plotting parameters are assigned for Pylab/matplotlib

# The plots are designed to adjust their axes according to the maximum and minimum values after taking the logarithms

if plotprofile ==1:
	pl.xlabel('Dimensionless Time   ') # set x-axis label
	pl.ylabel('Dimensionless Concentrations (log10)') # set y-axis label
	pl.title('\n' + '\n' + 'Chaotic Oregonator for CSTR at \n'+ 'kf = ' + str(kf)+' alpha = '+str(alpha)+' beta = ' +str(beta)+' C = ' + str(C)  + '\n' + ' stepsize = ' +str(stepsize) + ' iterations = ' + str(len(iters)) + '\n'+'Blue = HBrO2, Red = Ce(IV), Green = BrMA') # set plot title.plot(iters,np.log10(results))
	pl.axis([0,plotrange,np.min(np.log10(results)),np.max(np.log10(results))]) 
	pl.plot(iters,np.log10(results[:,0]),iters,np.log10(results[:,1]),iters,np.log10(results[:,2]))
	pl.show()

# The 3D plotting parameters are assigned for MayaVi

# these if statements just adjust the point sizes based on the number of iterations
if len(iters) <= 10:
	if len(iters) > 0:
		pointsize = .005
if len(iters) <=1000:
	if len(iters) > 10:
		pointsize = .004
if len(iters) <=2000:
	if len(iters) > 1000:
		pointsize = .002
if len(iters) <=10000:
	if len(iters) > 2000:
		pointsize = .0008
if len(iters) <= 30000:
	if len(iters) > 10000:
		pointsize = .0005
if len(iters) <= 50000:
	if len(iters) > 30000:
		pointsize = .0003
if len(iters) > 50000:
	pointsize = .0001


# Again, the axes are adjusting to maximum and minimum values from the array
if plotattractor == 1:
	points3d(np.log10(results[:,0]),np.log10(results[:,1]),np.log10(results[:,2]),iters,color=(0.,0.,0.),scale_factor = pointsize)
	axes(color=(1.,0.,0.),extent=[min(np.log10(results[:,0])),max(np.log10(results[:,0])),min(np.log10(results[:,1])),max(np.log10(results[:,1])),min(np.log10(results[:,2])),max(np.log10(results[:,2]))],ranges=[min(np.log10(results[:,0])),max(np.log10(results[:,0])),min(np.log10(results[:,1])),max(np.log10(results[:,1])),min(np.log10(results[:,2])),max(np.log10(results[:,2]))])
