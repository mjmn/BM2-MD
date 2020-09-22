#!/home/mniesen/Software/Anaconda/anaconda3/bin/python

import MDAnalysis as mda
import numpy as np
import numpy.linalg as npla
import sys
import math

### Functions used

def channelaxis(dirv,cv):
    '''
    Determine a normalized vector along the channel axis, baased on a direction and a center point.
    '''
    r = np.dot(dirv,cv)
    if r<0:
        return -1*cv/npla.norm(cv)
    else:
        return cv/npla.norm(cv)


def orderFunction(v1,v2):
    '''
    Calculate the contribution of a water to the S_{HH} order parameter.
    Assumes the two input vectors are normalized.
    '''
    cosTheta2 = np.dot(v1,v2)**2
    return 3*cosTheta2-1
    

# Hardcoded parameters
dx = 4.0 # Spacing for bins used in analysis, 1 Angstrom
dt = 100000.0 # Spacing used for time-averaging, 100ns
startT = 3001 # First frame to analyze, we left out the first 30ns
minx = -34.0 # Minimum channel coordinate
maxx = 34.0 # Maximum channel coordinate
## THE FOLLOWING SETTINGS ARE FOR CALCULATING A SINGLE CROSS-CHANNEL QUANTITY
#dx = 40.0
#minx = -12.0
#maxx = 12.0
## END OF CROSS-CHANNEL AVERAGING
cylR2 = 225.0 # Maximum distance, squared, from the channel axis for waters to be included in the analysis
dcos = 0.02 # Spacing used in collecting cosine theta distributions

# Define user inputs
if len(sys.argv)<2 or ("-h" in sys.argv):
    print("Use: ./water_order_analysis.py [TPR filename] [TRR filename]")
    exit()

TPR = sys.argv[1]
TRR = sys.argv[2]

# Initialize
u = mda.Universe(TPR,TRR)
bto = -1 # Initial time-bin
# Derrived output size
timestep = u.trajectory[0].dt
sx = math.ceil((maxx-minx)/dx) # number of bins along channel coordinate
st = math.floor(len(u.trajectory)*timestep/dt) # Number of bins in time
sc = math.ceil(2.0/dcos) # Number of bins in cosine distributions
# Make sure the time-blocks requested are not too long
if len(u.trajectory)*timestep < dt:
    print("Error; requested time-block is longer than the entire trajectory")
    exit()
# Pre-create zero arrays for analysis output
density = np.zeros((sx,st))
#densityO = np.zeros((sx,st))
HHorder = np.zeros((sx,st))
OHorder = np.zeros((sx,st))
cosOH = np.zeros((sc,sx),dtype=int)
cosHH = np.zeros((sc,sx),dtype=int)

# Group selections
center = u.select_atoms('protein and resname ALA and (resid 16 or resid 49 or resid 82 or resid 115) and name CA') # ALA 17 Calphas
direction = u.select_atoms('protein and resname TRP and (resid 22 or resid 55 or resid 88 or resid 121) and name CA') # TRP 23 Calphas
calphas = u.select_atoms('protein and name CA') # All Calphas
water_O = u.select_atoms('moltype TIP3 and type OT') # Water oxygens
water_H = u.select_atoms('moltype TIP3 and type HT') # Water hydrogens

# Collect data
for ts in u.trajectory[startT:]:
    #print(ts.time)
    # Current time-bin
    bt = int((ts.time - startT*timestep)/dt)
    # for avoiding incomplete time-blocks
    if bt>=st:#st:
        break
    # We will assign water depending on where along the channel axis it's oxygen is located
    cc = center.centroid() # Center of the channel
    dc = direction.centroid() - cc # Direction of the channel axis vector
    p1, p2, p3 = calphas.principal_axes() # Principal axes of all Calpha atoms
    caxis = channelaxis(dc,p3) # Calculate normalized vector along channel axis
    # Determine where each water is, and calculate it's contribtutions to the density/OHorder and HHorder
    for wi in range(len(water_O)):
        posO = water_O[wi].position - cc
        posOx = np.dot(posO,caxis)
        posOyz2 = npla.norm(posO)**2 - abs(posOx)**2
        if (posOx > minx) and (posOx < maxx) and (posOyz2 < cylR2):
            bx = int((posOx-minx)/dx)
            density[bx][bt] += 1
            # Add HHorder contribution
            hhvec = water_H[2*wi+1].position - water_H[2*wi].position # Vector between the hydrogens
            hhvec /= npla.norm(hhvec) # Normalize
            val = orderFunction(hhvec,caxis) # Calculate contribution of this water
            HHorder[bx][bt] += val # Add contribution, note we didn't divide by 2 yet
            cHH = np.dot(hhvec,caxis)
            cosHH[min([int((cHH+1.0)/dcos),(sc-1)])][bx] += 1
            # Add OHorder contribution
            ohvec = water_H[2*wi].position - water_O[wi].position # Vector between the hydrogens
            ohvec /= npla.norm(ohvec) # Normalize
            val = orderFunction(ohvec,caxis) # Calculate contribution of this water
            OHorder[bx][bt] += val # Add contribution, note we didn't divide by 2 yet
            cOH = np.dot(ohvec,caxis)
            cosOH[min([int((cOH+1.0)/dcos),(sc-1)])][bx] += 1
            # Add OHorder contribution for the second vector
            ohvec = water_H[2*wi+1].position - water_O[wi].position # Vector between the hydrogens
            ohvec /= npla.norm(ohvec) # Normalize
            val = orderFunction(ohvec,caxis) # Calculate contribution of this water
            OHorder[bx][bt] += val # Add contribution, note we didn't divide by 2 yet
            cOH = np.dot(ohvec,caxis)
            cosOH[min([int((cOH+1.0)/dcos),(sc-1)])][bx] += 1

# Process data and save
HHorder = np.divide(HHorder,2.0*density,out=np.zeros_like(HHorder),where=(density>0.5)) # Average, avoid zero division
OHorder = np.divide(OHorder,4.0*density,out=np.zeros_like(OHorder),where=(density>0.5)) # Average, avoid zero division
density = density/(dt/u.trajectory[0].dt) # Change it to the instantenous number density, instead of the density summed over time
# Write
np.savetxt('hhorder.txt',HHorder,fmt='%12.6f')
np.savetxt('ohorder.txt',OHorder,fmt='%12.6f')
np.savetxt('density.txt',density,fmt='%12.6f')
np.savetxt('cosine_thetaOH.txt',cosOH,fmt='%i')
np.savetxt('cosine_thetaHH.txt',cosHH,fmt='%i')
