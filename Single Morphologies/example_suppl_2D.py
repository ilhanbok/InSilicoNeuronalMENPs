#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Some plottin' functions used by the example scripts

Copyright (C) 2017 Computational Neuroscience Group, NMBU.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""
'''
example_suppl_2D.py
Description: Supplementary code for computing surrounding voxel-based magnetization
        for Allen Brain Atlas morphologies.
Usage: Download and patch desired morphologies from ABA, then run the batch file `loopMorpho.bat`
Author(s): Ilhan Bok, Xiaofei Qu
Last Modified: Jan. 19, 2022
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import neuron
import math

from mpl_toolkits.axes_grid1.anchored_artists import (AnchoredSizeBar)
import matplotlib.font_manager as fm

import matplotlib.cbook as cbook
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib_scalebar.scalebar import IMPERIAL_LENGTH

plt.rcParams.update({
    'axes.xmargin': 0.0,
    'axes.ymargin': 0.0,
})
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

import numpy as np 
import matplotlib.pyplot as plt

def nearest(x):
    return round(x / 20.0) * 20

def plot_mag_dist(cell, electrode, X, Y, Z, t_show_global, target_direction, xlist, ylist, zlist):
    '''
    plot the morphology and LFP contours, synaptic current and soma trace
    '''
    # some plot parameters
    tidx = np.where(cell.tvec == t_show_global)
    print(tidx)
    # contour lines:
    n_contours = 200
    n_contours_black = 40

    # This is the extracellular potential, reshaped to the X, Z mesh
    LFP = np.arcsinh(electrode.data[:, tidx]).reshape(X.shape)
    
    q = np.negative(np.gradient(LFP/1000, 20e-6))

    print('I\'m currently finding the magnetization distribution.')
    
    # This is the extracellular potential, reshaped to the X, Z mesh
    MAG = np.zeros((len(q[0]),len(q[0][0]),len(q[0][0][0])))
    for o in range(0,len(q[0])):
        for i in range(0,len(q[0][0])):
            for ii in range(0,len(q[0][0][0])):
                original_cos = math.cos(angle_between(target_direction,(q[0][o][i][ii],q[1][o][i][ii],q[2][o][i][ii])));
                MAG[o][i][ii] = (original_cos)*math.sqrt(q[0][o][i][ii]**2 + q[1][o][i][ii]**2 + q[2][o][i][ii]**2)/50*0.02
    
    delta = 20 # Electrode spacing (for consistency)
    scan_radii = delta # This is for now. One option is to make the radius dependent on neurite diameter
                    # Or we can just keep it at the electrode spacing: that also makes sense
    # Dictionary to store which voxels belong to which type of segment.
    # -> Currently there's a bit of overlap but this can be changed
    voxel_dict = {}
    # plot_morphology (synapses are not present for Allen Brain Atlas morphologies)
    for sec in neuron.h.allsec():
        category = sec.name().split('.')[-1].split('[')[0]
        if not category in list(voxel_dict):
            voxel_dict[category]  = []
        idx = cell.get_idx(sec.name())
        all_xs = np.r_[cell.x[idx, 0], cell.x[idx, 1][-1]]
        all_ys = np.r_[cell.y[idx, 0], cell.y[idx, 1][-1]]
        all_zs = np.r_[cell.z[idx, 0], cell.z[idx, 1][-1]]
        for i in range(0,len(all_xs)):
            segxyz = (all_xs[i],all_ys[i],all_zs[i])
            refxyz = (nearest(all_xs[i]),nearest(all_ys[i]),nearest(all_zs[i]))
            for dx in range(-1,2):
                for dy in range(-1,2):
                    for dz in range(-1,2):
                        if math.dist(segxyz,(refxyz[0]+dx*delta,refxyz[1]+dy*delta,refxyz[2]+dz*delta)) < scan_radii:
                            try:
                                curr_voxel = (xlist.index(refxyz[0]+dx*delta),
                                              ylist.index(refxyz[1]+dy*delta),
                                              zlist.index(refxyz[2]+dz*delta))
                                if curr_voxel not in voxel_dict[category]:
                                    voxel_dict[category].append(curr_voxel)
                            except ValueError:
                                print('Voxel wasn\'t inside our sampling zone (not a problem)')
                                # (i.e. just omit the voxel... :)
        #all_ds = np.r_[cell.d[idx, 0], cell.d[idx, 1][-1]] # use this for variable diameter voxel grabs
    
    print('keys:')
    print(list(voxel_dict))
    
    # Now iterate through all voxels and find all corresponding values
    mag_dict = { i : [] for i in list(voxel_dict) }
    for k in list(voxel_dict):
        for v in voxel_dict[k]:
            mag_dict[k].append(MAG[v[0]][v[1]][v[2]])
    return mag_dict

def plot_mag_voxel(cell, electrode, X, Y, Z, t_show_global, target_direction, xmin, xmax, zmin, zmax):
    '''
    plot the morphology and LFP contours, synaptic current and soma trace
    '''
    print('Here are all axes parameters')
    print(xmin)
    print(xmax)
    print(zmin)
    print(zmax)
    print('-----------------------------')
    # some plot parameters
    tidx = np.where(cell.tvec == t_show_global)
    print(tidx)
    # contour lines:
    n_contours = 200
    n_contours_black = 40

    # This is the extracellular potential, reshaped to the X, Z mesh
    LFP = np.arcsinh(electrode.data[:, tidx]).reshape(X.shape)
    
    q = np.negative(np.gradient(LFP/1000, 20e-6))
    
    print(q)
    print(q[0])
    

    # figure object
    fig = plt.figure(figsize=(12, 8))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None,
                        wspace=0.4, hspace=0.4)

    # Plot LFP around the cell with in color and with equipotential lines
    ax1 = fig.add_subplot(121, aspect='equal', frameon=False)
    #ax1 = fig.add_subplot(121, projection = '3d')

    # plot_morphology(plot_synapses=True)
    for sec in neuron.h.allsec():
        idx = cell.get_idx(sec.name())
        ax1.plot(np.r_[cell.x[idx, 0], cell.x[idx, 1][-1]],
                 np.r_[cell.z[idx, 0], cell.z[idx, 1][-1]],
                 color='orange')
    for i in range(len(cell.synapses)):
        ax1.plot([cell.synapses[i].x], [cell.synapses[i].z], '.',
                 markersize=10)

    # Figure formatting and labels
    fig.suptitle('example 1', fontsize=14)

    ax1.set_title('Electric Field Vectors at t=' + str(t_show_global) + ' ms', fontsize=12)
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    
    # This is the extracellular potential, reshaped to the X, Z mesh
    MAG = np.zeros((len(q[0]),len(q[0][0]),len(q[0][0][0])));
    for o in range(0,len(q[0])):
        for i in range(0,len(q[0][0])):
            for ii in range(0,len(q[0][0][0])):
                original_cos = math.cos(angle_between(target_direction,(q[0][o][i][ii],q[1][o][i][ii],q[2][o][i][ii])));
                MAG[o][i][ii] = (original_cos)*math.sqrt(q[0][o][i][ii]**2 + q[1][o][i][ii]**2 + q[2][o][i][ii]**2)/50*0.02

    # figure object
    fig = plt.figure(figsize=(12, 8))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None,
                        wspace=0.4, hspace=0.4)

    # Plot LFP around the cell with in color and with equipotential lines
    ax1 = fig.add_subplot(121, aspect='equal', frameon=False)

    # plot_morphology(plot_synapses=True)
    for sec in neuron.h.allsec():
        idx = cell.get_idx(sec.name())
        ax1.plot(np.r_[cell.x[idx, 0], cell.x[idx, 1][-1]],
                 np.r_[cell.z[idx, 0], cell.z[idx, 1][-1]],
                 color='k')
    for i in range(len(cell.synapses)):
        ax1.plot([cell.synapses[i].x], [cell.synapses[i].z], '.',
                 markersize=10)

    # contour lines
    m_range = [-1e-5, 1e-5]
    X, Z = np.mgrid[xmin:xmax:20, zmin:zmax:20]
    ct1 = ax1.pcolormesh(X, Z, MAG[:,1,:], vmin=np.min(m_range), vmax=np.max(m_range))

    # Plot soma potential
    ax3 = fig.add_subplot(224)
    ax3.plot(cell.tvec, cell.somav)

    # Figure formatting and labels
    fig.suptitle('example 1', fontsize=14)

    ax1.set_title('Change of Magnetization at t=' + str(t_show_global) + ' ms', fontsize=12)
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax1.set_yticks([])
    ax1.set_yticklabels([])

    ax3.set_title('somatic membrane potential', fontsize=12)
    ax3.set_ylabel('(mV)')
    ax3.set_xlabel('time (ms)')
    
    # Add colorbar
    
    print('Finalizing plot...\n')
    sm = plt.cm.ScalarMappable(cmap=matplotlib.cm.get_cmap(name='viridis'),
        norm=plt.Normalize(vmin=np.min(m_range), vmax=np.max(m_range)))
    plt.colorbar(sm,ax=ax1,format='%.0e')

    return fig

def find_tshow_global(cell):
    N = 20
    indexN = np.where(cell.tvec == 20)[0][0]
    filtsomav = cell.somav[0:indexN]
    return cell.tvec[np.where(filtsomav == np.max(filtsomav))[0][0]]