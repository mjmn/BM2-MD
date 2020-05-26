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


# Hardcoded parameters
dt = 1000.0 # Spacing used for time-averaging, 1ns
minx = -20.0 # Minimum channel coordinate
maxx = 20.0 # Maximum channel coordinate

# Define user inputs
if len(sys.argv)<2 or ("-h" in sys.argv):
    print("Use: ./water_transport.py [TPR filename] [TRR filename]")
    exit()

TPR = sys.argv[1]
TRR = sys.argv[2]

# Initialize
u = mda.Universe(TPR,TRR)
bto = -1 # Initial time-bin
# Derrived output size
st = math.floor(len(u.trajectory)*u.trajectory[0].dt/dt) # Number of bins in time
# Make sure the time-blocks requested are not too long
if len(u.trajectory)*u.trajectory[0].dt < dt:
    print("Error; requested time-block is longer than the entire trajectory")
    exit()
# Pre-create zero arrays for analysis output
transported = np.zeros((st,3)) # Time, Transport to N-term, Transport to C-term
passageTimesN = [] # Empty list to contain the passage times in the N-terminal direction
passageTimesC = [] # Empty list to contain the passage times in the C-terminal direction

# Group selections
center = u.select_atoms('protein and resname ALA and (resid 16 or resid 49 or resid 82 or resid 115) and name CA') # ALA 17 Calphas
direction = u.select_atoms('protein and resname TRP and (resid 22 or resid 55 or resid 88 or resid 121) and name CA') # TRP 23 Calphas
calphas = u.select_atoms('protein and name CA') # All Calphas
water_O = u.select_atoms('moltype TIP3 and type OT') # Water oxygens

## Initialize a vector that tracks the location and time of entering the channel for water molecules
waterLoc = np.zeros((len(water_O),1),dtype=int)
waterTimeStamp = np.zeros((len(water_O),1),dtype=float)

# Collect data
for ts in u.trajectory:
    # Current time-bin
    bt = int(ts.time/dt)
    # for avoiding incomplete time-blocks
    if bt>=st:#st:
        break
    transported[bt][0] = bt*dt
    # We will assign water depending on where along the channel axis it's oxygen is located
    cc = center.centroid() # Center of the channel
    dc = direction.centroid() - cc # Direction of the channel axis vector
    p1, p2, p3 = calphas.principal_axes() # Principal axes of all Calpha atoms
    caxis = channelaxis(dc,p3) # Calculate normalized vector along channel axis
    # For each water see if it passed a boundary
    for wi in range(len(water_O)):
        posO = water_O[wi].position - cc
        posOx = np.dot(posO,caxis)
        if (posOx<minx): # This water is at the N-terminal half of the channel
            if waterLoc[wi]==4: # It started from the C-terminal half and went in the channel
                transported[bt][1] += 1 # Add a count
                passageTimesN.append(ts.time-waterTimeStamp[wi]) # Collect data on how long it took this water to transport
            waterLoc[wi] = 1 # Label this water as being on the N-terminal half
        elif (posOx>maxx): # This water is at the C-terminal half of the channel
            if waterLoc[wi]==3: # It started N-terminal and entered the channel
                transported[bt][2] += 1 # add a count
                passageTimesC.append(ts.time-waterTimeStamp[wi]) # Collect data on how long it took this water to transport
            waterLoc[wi] = 2 # Label this water as being on the C-terminal half
        else: # Water is in the channel
            if waterLoc[wi]==1:
                waterLoc[wi] = 3 # In the channel, entered from N-term
                waterTimeStamp[wi] = ts.time # Store the time at which the water entered the channel
            elif waterLoc[wi]==2:
                waterLoc[wi] = 4 # In the channel, entered from C-term
                waterTimeStamp[wi] = ts.time # Store the time at which the water entered the channel
# Process data and save
# Convert passagetimes from list to np array and save
np.savetxt('passageTimesN.txt',np.array(passageTimesN),fmt='%12.2f')
np.savetxt('passageTimesC.txt',np.array(passageTimesC),fmt='%12.2f')
# Write
np.savetxt('transported.txt',transported,fmt='%i')
