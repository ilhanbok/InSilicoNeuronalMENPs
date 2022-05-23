'''
QuantifyMag_MPIvoxel_L23.py
Description: Runs a biophysical simulation on Allen Brain Atlas morphologies to
             compute surrounding voxel-based magnetization.
Usage: Download and patch desired morphologies from ABA, then run the batch file `loopMorpho.bat`
Author(s): Ilhan Bok, Xiaofei Qu
Last Modified: Jan. 19, 2022
'''

import numpy as np
import sys
import matplotlib.pyplot as plt
import LFPy
import neuron
import math

def rndu(x, delta):
    return int(math.ceil(x / delta)) * delta
    
def rndd(x, delta):
    return int(math.floor(x / delta)) * delta
    
def rnd(x, delta):
    return int(round(x / delta)) * delta

#compile mod files every time, because of incompatibility with Mainen96 files:
neuron.nrn_dll_loaded.append('nrnmech.dll')

cellname = sys.argv[1]
# define cell parameters used as input to cell-class
cellParameters = {
   'morphology'    : cellname + '.asc',
    'templatefile'  : ['templatePC.hoc',
                      ],
    'templatename'  : 'cADpyr231_L6_TPC_L4_0cb1e9aa6b',
    'templateargs'  : cellname,
    'passive' : False,
    'nsegs_method' : None,
    'dt' : 2**-6,
    'tstart' : 0,
    'tstop' : 50,
    'v_init' : -60,
    'celsius': 34,
    'pt3d' : True,
}
# delete old sections from NEURON namespace
LFPy.cell.neuron.h("forall delete_section()")

# Initialize cell instance, using the LFPy.Cell class
cell = LFPy.TemplateCell(**cellParameters)
cell.set_rotation(x=4.729, y=-3.166, z=0)

# Override passive reversal potential, AP is generated
for sec in cell.allseclist:
    for seg in sec:
        seg.e_pas = 100

# Interval for recording electrode placement
delta = 20
num_radii = 1
cellx = cell.x.flatten()
celly = cell.y.flatten()
cellz = cell.z.flatten()

xmid = rnd((rndd(min(cellx),delta)-2*delta + rndu(max(cellx),delta)+2*delta)/2,delta)
xmin = rndd(min(cellx),delta)-2*delta
xmax = rndu(max(cellx),delta)+2*delta+1
ymin = rndd(min(celly),delta)-2*delta
ymax = rndu(max(celly),delta)+2*delta+1
zmin = rndd(min(cellz),delta)-2*delta
zmax = rndu(max(cellz),delta)+2*delta+1

# Generate the grid in xz-plane over which we calculate local field potentials
X, Y, Z = np.mgrid[xmin:xmax:delta,
                   ymin:ymax:delta,
                   zmin:zmax:delta]

# Print minimum and maximum x, y, and z for scaling by MATLAB
f = open('CellSizes.out', 'a')
f.write(str(cellname) + ' ' + str(xmin) + ' ' + str(xmax) + ' ' + str(ymin) + ' ' + str(ymax)
    + ' ' + str(zmin) + ' ' + str(zmax) + '\n')
f.close()

# define parameters for extracellular recording electrode, using optional method
electrodeParameters = {
    'sigma' : 0.3,          # extracellular conductivity
    'x' : X.flatten(),      # x,y,z-coordinates of contacts
    'y' : Y.flatten(),
    'z' : Z.flatten(),
    'method' : 'root_as_point',  #sphere source soma segment
    'N' : np.array([[0, 1, 0]]*X.size), #surface normals
    'r' : 2.5,              # contact site radius
    'n' : 20,               # datapoints for averaging
}


# create extracellular electrode object for LFPs on grid
electrode = LFPy.RecExtElectrode(cell=cell, **electrodeParameters)

# perform NEURON simulation
# Simulated results saved as attribute `data` in the RecExtElectrode instance
cell.simulate(probes=[electrode])

from example_suppl_2D import plot_mag_dist, plot_mag_voxel, find_tshow_global
t_show_global = find_tshow_global(cell) # note that this is for times less than 20 ms (to filter out saturation)

# Create magnetization voxel plots
vector_str = 'XYZ' #Include static field directions in {X,Y,Z} to compute the results for
for currdir in vector_str:
    print('Computing for ' + currdir)
    dirvec = (1,0,0)
    # Map direction to vector
    if currdir == "Y":
        dirvec = (0,1,0)
    elif currdir == "Z":
        dirvec = (0,0,1)
    fig = plot_mag_voxel(cell, electrode, X, Y, Z, t_show_global, dirvec, xmin, xmax, zmin, zmax)
    fig.savefig('PC_LFPy(' + cellname + ')' + currdir + '.svg', dpi=350)
    mag_dict = plot_mag_dist(cell, electrode, X, Y, Z, t_show_global,dirvec,
                    list(range(xmin,xmax,delta)),
                    list(range(ymin,ymax,delta)),
                    list(range(zmin,zmax,delta)))
    f = open( 'PC_VAR_QuantifyMagDist(' + cellname + ')' + currdir + '.out', 'w')
    f.write(repr(mag_dict))
    f.close()