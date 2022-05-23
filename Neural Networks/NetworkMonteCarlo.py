'''
NetworkMonteCarlo.py
Description: Read or generate network burst activity, and plot magnetization
    heatmap and Monte Carlo data
Usage: Specify whether reading NetPyNE simulation from file, then run script
Author(s): Ilhan Bok
Last Modified: Mar. 1, 2022
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from statistics import stdev

import json, math
import matplotlib.font_manager as fm
import scipy.interpolate as sciterp
import matplotlib.cbook as cbook
import sys
    
# https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
import numpy as np

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def map_magnitude(v1,v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)*math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)/50*0.02*1e6
    
def plot_mag(direction,LFP, X, Y, Z, slice_idx, dt, lim):
    '''
    plot the morphology and LFP contours, synaptic current and soma trace
    '''
    # # contour lines:
    n_contours = 200
    n_contours_black = 40
    
    print('Here are X, Y, and Z')
    print(X)
    print(Y)
    print(Z)
    
    # Take negative gradient of LFP to find electric field
    q = np.negative(np.gradient(LFP/1000, 5e-6))
    
    print('Applying linear magnetization algorithm...\n')
    
    # This is the extracellular potential, reshaped to the X, Z mesh
    #MAG = np.arcsinh(electrode.data[:, tidx]).reshape(X.shape)
    MAG = np.zeros((len(q[0]),len(q[0][0]),len(q[0][0][0])));
    for o in range(0,len(q[0])):
        for i in range(0,len(q[0][0])):
            for ii in range(0,len(q[0][0][0])):
                original_cos = math.cos(angle_between(direction,(q[0][o][i][ii],q[1][o][i][ii],q[2][o][i][ii])));
                MAG[o][i][ii] = (original_cos)*math.sqrt(q[0][o][i][ii]**2 + q[1][o][i][ii]**2 + q[2][o][i][ii]**2)/50*0.02

    # figure object
    plt.figure(figsize=(12, 8))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                        wspace=0.4, hspace=0.4)

    # Plot LFP around the cell with in color and with equipotential lines
    ax1 = plt.subplot(121, aspect='equal', frameon=False)

    slice_num = 115
    GLOBAL_SLICE = MAG[:,slice_num,:]
    print('MAG shape')
    print(str(MAG.shape))
    print(str(GLOBAL_SLICE.shape))
    m_range = lim
    # contour lines
    ct1 = ax1.contourf(X, Z, GLOBAL_SLICE, n_contours)
    ct1.set_clim((np.min(m_range), np.max(m_range)))
    ax1.contour(X, Z, GLOBAL_SLICE, n_contours_black, colors='k')
    print(np.min(GLOBAL_SLICE))
    print(np.max(GLOBAL_SLICE))

    # Figure formatting and labels
    plt.suptitle('example 1', fontsize=14)

    ax1.set_title('Change of Magnetization at t=' + str(slice_idx) + ' ms', fontsize=12)
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    
    # Add colorbar
    
    print('Finalizing plot...\n')
    sm = plt.cm.ScalarMappable(cmap=matplotlib.cm.get_cmap(name='jet'),
        norm=plt.Normalize(vmin=np.min(m_range), vmax=np.max(m_range)))
    plt.colorbar(sm,ax=ax1,format='%.0e')

    return plt

def plot_mag_timeflow_MONTECARLO(direction,LFP, X, Y, Z, slice):
    np.set_printoptions(threshold=sys.maxsize)
    plt.figure()
    '''
    plot the morphology and LFP contours, synaptic current and soma trace
    '''
    ''''''
    read_ds = True
    ''''''
    speed = 65 # (cm/s)
    D = 15 # (um^2/s)
    dt_orig = 0.1 # (ms)
    y_length = 700 # (um)
    numMC = 40e3 # number of total monte carlo runs, 40,000 should on average cover each square micron 
    ####### Calculate and rescale the parameters #######
    numMC = int(numMC) # this is needed to convert from float
    speed = speed * 10 # (um/ms)
    D = D / 1e12 # (m^2/s)
    dt = dt_orig / 1000 # (s)
    dv = math.sqrt(6*D*dt) # (m)
    dv = dv * 1e6 # (um) <- our unit
    # Center the particle maximally around the network burst peak
    tstart = slice - y_length / speed / 2
    num_of_steps = math.floor((y_length / speed) / dt_orig)
    if read_ds:
        with open('ds'+str(slice)+'.out','r') as f:
            ds = eval(f.read().replace('array', 'np.array').replace('nan', 'np.nan'))
            dts = ds[0]
            dxs = ds[1]
            dys = ds[2]
            dzs = ds[3]
            mag_4d = ds[4]
    else:

        # Do a 4-D interpolation on every particle
        # First, get the 4-dimensional LFP profile
        arr_3d = (np.array(LFP[0])).reshape((len(X),len(Y),len(Z))).transpose()
        q = np.negative(np.gradient(arr_3d/1000, 5e-6))
        mag_4d = np.zeros((150,
                    len(q[0]),
                    len(q[0][0]),
                    len(q[0][0][0])
                    ))
        for t in range(0,150):
            arr_3d = (np.array(LFP[t])).reshape((len(X),len(Y),len(Z))).transpose()
            q = np.negative(np.gradient(arr_3d/1000, 5e-6))
            mag_3d = np.zeros((
                    len(q[0]),
                    len(q[0][0]),
                    len(q[0][0][0])))
            for o in range(0,len(q[0])):
                for i in range(0,len(q[0][0])):
                    for ii in range(0,len(q[0][0][0])):
                        mag_3d[o][i][ii] = map_magnitude((1,0,0),(q[0][o][i][ii],q[1][o][i][ii],q[2][o][i][ii]))
            mag_4d[t] = mag_3d
        print('sizes:')
        print(mag_4d.shape)
        print(mag_3d.shape)
        # Idea: interpolate for all monte carlo simulations at once
        dxs = np.array([])
        dys = np.array([])
        dzs = np.array([])
        for i in range(0,numMC):
            dxs = np.append(dxs,np.cumsum(dv*(1-2*np.random.rand(1,num_of_steps)))+np.random.rand()*200)
            dys = np.append(dys,np.cumsum(dv*(1-2*np.random.rand(1,num_of_steps))+speed*dt_orig))
            dzs = np.append(dzs,np.cumsum(dv*(1-2*np.random.rand(1,num_of_steps)))+np.random.rand()*200)
        dxs[dxs >= 200] = np.nan
        dzs[dzs >= 200] = np.nan
        print(dxs)
        dts = np.tile(np.cumsum(np.full(num_of_steps,dt_orig))-dt_orig+tstart,numMC)
        with open('ds'+str(slice)+'.out','w') as f:
            f.write(repr((dts, dxs, dys, dzs, mag_4d)))
    ####### Now perform 4-D interpolation of data #######
    read_ls = True
    if read_ls:
        with open('ls'+str(slice)+'.out','r') as f:
            ls = eval(f.read().replace('array', 'np.array').replace('nan','np.nan'))
    else:
        ls = sciterp.interpn((list(range(0,150)),X,Y,Z),mag_4d,np.transpose(np.array([dts,dxs,dys,dzs])),bounds_error=False)
        with open('ls'+str(slice)+'.out','w') as f:
            f.write(repr(ls))
    for i in range(0,numMC):
        plt.plot(list(range(0,num_of_steps)),ls[i*num_of_steps:(i+1)*num_of_steps],color='black', alpha=0.01,zorder=1)
    standevs = []
    lsr = ls.reshape(-1,10)
    print(str(lsr.shape))
    truemean = np.nanmean(lsr,axis=0)
    for ti in range(0,num_of_steps):
        standevs.append(np.nanstd(lsr[:,ti]));
    print(standevs)
    with open('truemean'+str(slice)+'.out','w') as f:
        f.write(repr(truemean))
    with open('standevs'+str(slice)+'.out','w') as f:
        f.write(repr(standevs))
    plt.plot(list(range(0,num_of_steps)),truemean, linewidth=3.0, color='red',zorder=3)
    plt.fill_between(list(range(0,num_of_steps)), truemean-standevs, truemean+standevs,
        alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848', zorder=2)
    plt.ylim([-1e4, 1e4])
    plt.yscale('symlog', linthreshy=1e2)
    plt.xlim([0, 9])
    plt.xticks(list(range(0,10)), np.round(dts[0:10],2))
    return plt

######################
from netpyne import specs, sim, analysis

read_from_file = False
filename = 'NetLFPDataOverrideMONTECARLO_CORRECT.txt'

X = range(0, 201, 5)
Y = range(0, 701, 5)
Z = range(0, 201, 5)

''' Specify desired applied magnetic field direction here {X,Y,Z} '''
currdir = "X"

dirvec = (1,0,0)
# Map the desired direction computation to the appropriate vector
if currdir == "Y":
    dirvec = (0,1,0)
elif currdir == "Z":
    dirvec = (0,0,1)

# Hopefully I've already run this script and can skip the recomputation of LFPs
if read_from_file:
    json_read = open(filename, 'r')
    LFP = np.array(json.load(json_read))
    print('Congrats! Able to load LFP from file.')
else:

    # Network parameters
    netParams = specs.NetParams()  # object of class NetParams to store the network parameters

    netParams.sizeX = 200 # x-dimension (horizontal length) size in um
    netParams.sizeY = 700 # y-dimension (vertical height or cortical depth) size in um
    netParams.sizeZ = 200 # z-dimension (horizontal length) size in um
    netParams.propVelocity = 100.0 # propagation velocity (um/ms)
    netParams.probLengthConst = 150.0 # length constant for conn probability (um)

    ## Population parameters
    netParams.popParams['E2'] = {'cellType': 'E', 'numCells': 495-400, 'yRange': [0,150], 'cellModel': 'HH'}
    netParams.popParams['I2'] = {'cellType': 'I', 'numCells': 495+400, 'yRange': [0,150], 'cellModel': 'HH'}
    netParams.popParams['E3'] = {'cellType': 'E', 'numCells': 588-400, 'yRange': [150,500], 'cellModel': 'HH'}
    netParams.popParams['I3'] = {'cellType': 'I', 'numCells': 588+400, 'yRange': [150,500], 'cellModel': 'HH'}
    netParams.popParams['E4'] = {'cellType': 'E', 'numCells': 570-400, 'yRange': [500,700], 'cellModel': 'HH'}
    netParams.popParams['I4'] = {'cellType': 'I', 'numCells': 570+400, 'yRange': [500,700], 'cellModel': 'HH'}

    ## Cell property rules
    netParams.loadCellParamsRule(label='CellRule', fileName='IT2_reduced_cellParams_mod.json')
    netParams.cellParams['CellRule']['conds'] = {'cellType': ['E','I']}

    ## Synaptic mechanism parameters
    netParams.synMechParams['exc'] = {'mod': 'Exp2Syn', 'tau1': 0.8, 'tau2': 5.3, 'e': 0}  # NMDA synaptic mechanism
    netParams.synMechParams['inh'] = {'mod': 'Exp2Syn', 'tau1': 0.6, 'tau2': 8.5, 'e': -75}  # GABA synaptic mechanism

    ## Cell connectivity rules
    netParams.connParams['E->all'] = {
      'preConds': {'cellType': 'E'}, 'postConds': {'y': [100,1000]},  #  E -> all (100-1000 um)
      'probability': 0.1,                  # probability of connection
      'weight': '5*5.0*post_ynorm',         # synaptic weight
      'delay': 'dist_3D/propVelocity',      # transmission delay (ms)
      'synMech': 'exc'}                     # synaptic mechanism

    netParams.connParams['I->E'] = {
      'preConds': {'cellType': 'I'}, 'postConds': {'pop': ['E2','E4','E5']},       #  I -> E
      'probability': '0.4*exp(-dist_3D/probLengthConst)',   # probability of connection
      'weight': 5*1.0,                                      # synaptic weight
      'delay': 'dist_3D/propVelocity',                      # transmission delay (ms)
      'synMech': 'inh'}                                     # synaptic mechanism


    # Simulation configuration
    simConfig = specs.SimConfig()        # object of class SimConfig to store simulation configuration
    simConfig.duration = 151           # Duration of the simulation, in ms
    simConfig.dt = 0.1                # Internal integration timestep to use
    simConfig.verbose = False            # Show detailed messages
    simConfig.recordStep = 1             # Step size in ms to save data (eg. V traces, LFP, etc)
    simConfig.filename = 'net_lfp'   # Set file output name

    simConfig.recordLFP = [[x, y, z] for z in Z for y in Y for x in X]
    simConfig.analysis['plotRaster'] = {'orderBy': 'y', 'orderInverse': True, 'saveFig':'raster_ydir_new200.eps', 'figSize': (9,3), 'timeRange':[0,150]}

    # Create network and run simulation
    sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig)

    print('Initialization completed.\n')

    edata = analysis.lfp.plotLFP(plots=[])
    print('LFP plotted.\n')
    LFP = edata[1].get('LFP')
    with open('./' + filename, 'w') as json_file:
      json.dump(LFP.tolist(), json_file)

# plot LFPs and find electric field
slice_idx = 130;

print('Finding arr_3d')
arr_1d = LFP[slice_idx];
arr_3d = np.array(arr_1d).reshape((len(X),len(Y),len(Z))).transpose();
print('DONE (finding arr_3d)')

print('Plotting magnetization change...')
''' Specify the desired time slices (in ms) here '''
slices = [37, 69, 101, 129]
''' Specify which directions (some string permutation of X,Y,Z) here '''
dirsel = 'X'
for currdir in dirsel:
    print('Computing for ' + currdir)
    dirvec = (1,0,0)
    if currdir == "Y":
        dirvec = (0,1,0)
    elif currdir == "Z":
        dirvec = (0,0,1)
    for slice in slices:
        plt = plot_mag_timeflow_MONTECARLO(dirvec, LFP, list(X), list(Y), list(Z), slice)
        plt.savefig('BrownianMonteCarlo'+currdir+'_slice=' + str(slice)+'.png', format="png", dpi=600)
print('All done.\n')