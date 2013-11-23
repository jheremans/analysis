import os, sys, time
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D

from analysis.lib.fitting import fit


# TODO implement formatting options
def plot_fit1d(res, fit_xvals=None,fit_num_points=100,ax=None, ret=None, **kw):
    '''
    function to plot a fitresult as returned by analysis.lib.fitting.fit.fit1d().
    '''

    # know keywords:
    print_info = kw.pop('print_info', True)
    info_xy = kw.pop('info_xy', 'auto')
    plot_data = kw.pop('plot_data', True)
    
    if res == None:
        return False
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if plot_data:
        if 'yerr' in res.keys():
            ax.errorbar(res['x'],res['y'],fmt='o',yerr=res['yerr'])
        else:    
            ax.plot(res['x'], res['y'], 'o')
    if fit_xvals == None:
        fit_xvals=np.linspace(res['x'][0],res['x'][-1],fit_num_points)
    ax.plot(fit_xvals, res['fitfunc'](fit_xvals), '-', lw=2)

    if print_info:
        params_str = res['fitfunc_str'] + '\n' + fit.str_fit_params(res)
        
        if info_xy == 'auto':
            info_x = ax.get_xlim()[0] + (ax.get_xlim()[-1]-ax.get_xlim()[0])*0.02
            info_y = ax.get_ylim()[0] + (ax.get_ylim()[-1]-ax.get_ylim()[0])*0.02
        else:
            info_x = info_xy[0]
            info_y = info_xy[1]

        ax.text(info_x, info_y, params_str, size='x-small',
                color='k', ha='left', va='bottom', 
                bbox=dict(facecolor='white', alpha=0.5))
    
    if ret == None:
        return
    if ret == 'ax':
        return ax

    return

#### Density matrix and state plotting
def dm(dm, **kw):
    '''
    Plot the real and imaginary parts of a density matrix as 3d bar charts.
    '''

    dmlabel = kw.pop('dmlabel', r'$\rho$')

    dmr = dm.real
    dmi = dm.imag

    fig = plt.figure(figsize=(12,4))
    axr = fig.add_subplot(121, projection='3d')
    axi = fig.add_subplot(122, projection='3d')

    lx= len(dm[0])
    ly= len(dm[:,0])
    xticks = [(r'$| %d \rangle$' % i) for i in range(lx)]
    yticks = [(r'$\langle %d |$' % i) for i in range(ly)]

    matrix(dmr, axr, title=r'$\Re$'+dmlabel, 
            xticks=xticks, yticks=yticks, **kw)
    matrix(dmi, axi, title=r'$\Im$'+dmlabel, 
            xticks=xticks, yticks=yticks, **kw)

def matrix(mat, ax=None, **kw):
    '''
    plot a matrix as 3d bar chart.
    '''
    xlabel = kw.pop('xlabel', '')
    ylabel = kw.pop('ylabel', '')
    xticks = kw.pop('xticks', [])
    yticks = kw.pop('yticks', [])
    title = kw.pop('title', '')
    zlabel = kw.pop('zlabel', '')
    ret = kw.pop('ret', None)

    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
    lx= len(mat[0])            # Work out matrix dimensions
    ly= len(mat[:,0])
    xpos = np.arange(0,lx,1)    # Set up a mesh of positions
    ypos = np.arange(0,ly,1)
    xpos, ypos = np.meshgrid(xpos+0.25, ypos+0.25)

    xpos = xpos.flatten()   # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)

    dx = 0.75 * np.ones_like(zpos)
    dy = dx.copy()
    dz = mat.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, alpha=0.5, color='RoyalBlue', 
            zsort='average')

    if xticks == []:
        xticks = [str(i) for i in np.arange(lx)]
    if yticks == []:
        yticks = [str(i) for i in np.arange(ly)]
    
    ax.w_xaxis.set_ticks(np.arange(lx) + 0.675)
    ax.w_yaxis.set_ticks(np.arange(ly) + 0.675)
    ax.w_xaxis.set_ticklabels(xticks)
    ax.w_yaxis.set_ticklabels(yticks)
    ax.set_zlim3d([-1,1])
    ax.set_zlabel(zlabel)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if ret == 'ax':
        return ax

    return

def msmnt_outcomes(outcomes, u_outcomes=None,
        renormalize=True, ax=None, **kw):
    '''
    plot a bar plot with measurement outcomes.

    arguments:
    - renormalize : bool (True)
      if True, assume data range is 0..1; renormalize to -1..1
      else leave as is (and assume range is 0..1 per default)

    known kws:
    - fc : color ('RoyalBlue')
      face color of the data bars
    - xticks : array (None)
      xaxis tick labels for the outcome bars
    - reference : array (None)
      reference values, to be shown for comparison (renormalized in the same way as
      the data)
    '''

    fc = kw.pop('fc', 'RoyalBlue')
    xticks = kw.pop('xticks', None)
    reference = kw.pop('reference', None)
    
    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    # get the correct dimensions
    l = len(outcomes)

    # renormalization
    if renormalize:
        outcomes = outcomes*2 - 1.
        if u_outcomes != None:
            yerr = 2*u_outcomes.copy()
        if reference != None:
            ref = reference.copy()*2 - 1.
    else:
        if u_outcomes != None:
            yerr = u_outcomes.copy()
        if reference != None:
            ref = reference.copy()

    # plot the data
    pos = np.arange(l)
    w = 0.8

    if reference != None:
        ax.bar(pos, ref, w, color='w')
    
    kw = {}
    if u_outcomes != None:
        kw['yerr'] = yerr
        kw['ecolor'] = 'k'
        kw['capsize'] = 5

    ax.bar(pos, outcomes, w, color=fc, **kw)

    # format the plot
    if not renormalize:
        yticks = [0,0.25,0.5,0.75,1]
        ylim = (-0.1,1.1)
        ax.add_line(Line2D([-0.1, l-0.1], [0,0], color='k'))
    else:
        yticks = [-1,-0.5,0,0.5,1]
        ylim = (-1.1,1.1)
        ax.add_line(Line2D([-0.1, l-0.1], [0,0], color='k'))

    if xticks == None:
        xticks = [str(i) for i in range(l)]

    ax.set_xticks(pos + 0.4)
    ax.set_xlim(-0.1, l-0.1)
    ax.set_xticklabels(xticks)

    ax.set_yticks(yticks)
    ax.set_ylim(ylim)




    


