#!/home/mniesen/Software/Anaconda/anaconda3/bin/python

import MDAnalysis as mda
import numpy as np
import numpy.linalg as npla
import sys
import math

### Functions used

def calc_normalized_dipole(O, H1, H2):
    '''
    Given positions for an Oxygen and two hydrogens, return the normalized dipole vector.
    NOTE: I don't check for periodic boundaries here, my trajectories are saved to keep whole molecules.
    '''
    dvec = (H1-O)+(H2-O) # Direction doesn't matter for correlation function, just note the convention if using this function elsewhere
    return dvec/npla.norm(dvec)

def calcMSD(r1, r2, box):
    '''
    Return the mean-squared displacement given two position vectors.
    '''
    dr = np.abs(r2-r1)
    for i in range(len(dr)):
        if dr[i]>box[i]/2.0: # Bigger displacement than half a box-length shouldn't happen on a 500ps timescale
            dr[i] = box[i]-dr[i]
    return npla.norm(dr)**2

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

def assignToBin(rw, dz, origin, zaxis, minz, maxz, radial_cutoff2):
    '''
    Assign a particle to a bin along the channel axis
    '''
    posO = rw - origin
    posOx = np.dot(posO, zaxis)
    posOyz2 = npla.norm(posO)**2 - abs(posOx)**2
    if (posOx > minz) and (posOx < maxz) and (posOyz2 < radial_cutoff2):
        bx = int((posOx-minz)/dz)
    else:
        bx = -1 # It isn't in the binned-range
    return bx
    

# Hardcoded parameters
maxTau = 20 # At most 50 timesteps (10ps each)
dx = 4.0 # Spacing for bins used in analysis, 1 Angstrom
minx = -34.0 # Minimum channel coordinate
maxx = 34.0 # Maximum channel coordinate
## THE FOLLOWING SETTINGS ARE FOR CALCULATING A SINGLE CROSS-CHANNEL QUANTITY
#dx = 24.0
#minx = -12.0
#maxx = 12.0
## END OF CROSS-CHANNEL AVERAGING
cylR2 = 225.0 # Maximum distance, squared, from the channel axis for waters to be included in the analysis

# Define user inputs
if len(sys.argv)<2 or ("-h" in sys.argv):
    print("Use: ./water_order_analysis.py [TPR filename] [TRR filename] {startT} {endT}\nArguments in {} are optional.\n")
    exit()

TPR = sys.argv[1]
TRR = sys.argv[2]
if len(sys.argv)>2:
    startT = int(sys.argv[3])
    endT = int(sys.argv[4])
else:
    startT = 0
    endT = -1

# Initialize
u = mda.Universe(TPR,TRR)
# Derrived output size
sx = math.ceil((maxx-minx)/dx) # number of bins along channel coordinate
# Pre-create zero arrays for analysis output
crot_mat = np.zeros((sx, maxTau+1, 3), dtype=float)

# Group selections
center = u.select_atoms('protein and resname ALA and (resid 16 or resid 49 or resid 82 or resid 115) and name CA') # ALA 17 Calphas
direction = u.select_atoms('protein and resname TRP and (resid 22 or resid 55 or resid 88 or resid 121) and name CA') # TRP 23 Calphas
calphas = u.select_atoms('protein and name CA') # All Calphas
water_O = u.select_atoms('moltype TIP3 and type OT') # Water oxygens
water_H = u.select_atoms('moltype TIP3 and type HT') # Water hydrogens

# Empty lists for tracking historical water dipole and bin assignment
dipole_history = [ [] for i in range(len(water_O)) ]
bin_history = [ [] for i in range(len(water_O)) ]

# Collect data
for ts in u.trajectory[startT:endT]:
    box = ts.dimensions[0:3] # Used to fixed for periodic boundary conditions in MSD calculation
    print(box)
    # We will assign water depending on where along the channel axis it's oxygen is located
    cc = center.centroid() # Center of the channel
    dc = direction.centroid() - cc # Direction of the channel axis vector
    p1, p2, p3 = calphas.principal_axes() # Principal axes of all Calpha atoms
    caxis = channelaxis(dc,p3) # Calculate normalized vector along channel axis
    # If we're in a new time-bin, determine where along the channel axis each water is
    for wi in range(len(water_O)):
        rw = water_O[wi].position # Current position of the water-molecule
        # Check which channel bin the water is in at the moment
        bw = assignToBin(rw, dx, cc, caxis, minx, maxx, cylR2)
        #print("Current: ",rw, bw, "\n")
        # Calculate the normalized dipole vector for this water
        rh1 = water_H[wi*2].position
        rh2 = water_H[wi*2+1].position
        n_dipole = calc_normalized_dipole(rw, rh1, rh2)
        # Store the information for this water
        dipole_history[wi].insert(0,n_dipole)
        bin_history[wi].insert(0,bw)
        # Now add contributions to MSD
        for tau in range(len(dipole_history[wi])):
            # Get the position and bin data for this water, tau timesteps ago
            if tau ==  maxTau: # Forget data that is too far in the past
                n_dipoleT = dipole_history[wi].pop()
                bwT = bin_history[wi].pop()
            else:
                n_dipoleT = dipole_history[wi][tau]
                bwT = bin_history[wi][tau]
            #print("Previous: ", rwT, bwT, "\n")
            if bwT!=-1: # If the water started out inside the channel, else there is no contribution
                S = 0.5*orderFunction(n_dipole, n_dipoleT) # Calculate (3/2)*mu(t)*mu(t+tau)-1
                Crot = np.dot(n_dipole, n_dipoleT) # Simple dot-product
             #   print("MSD: ", msd, "\n")
                crot_mat[bwT, tau, 0] += S
                crot_mat[bwT, tau, 1] += Crot
                crot_mat[bwT, tau, 2] += 1

# Process data and save
#msd_mat[:,:,0] = np.divide(msd_mat[:,:,0], msd_mat[:,:,1],out=np.zeros_like(msd_mat[:,:,0]),where=(msd_mat[:,:,1]>0.5)) # Average, avoid zero division
crot_mat[:,:,0] /= crot_mat[:,:,2]
crot_mat[:,:,1] /= crot_mat[:,:,2]
# Write
fno = 'Crot_vs_channelaxis_%i-%i.txt' % (startT, endT)
fo = open(fno,'ab') # File for writing
for b in range(sx):
    np.savetxt(fo, crot_mat[b, :, :], fmt='%12.3f')
fo.close()
