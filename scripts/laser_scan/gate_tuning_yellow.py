import numpy as np
from scipy.interpolate import griddata

def get_data(fn):
    d = np.loadtxt(fn)
    return d

def get_yaxis_vals(d, ycol=3):
    return np.unique(d[:,ycol])

def make_grid_data(d, ycol=3, xpts=200):
    yvals = d[:,ycol]
    xvals = d[:,1]
    zvals = d[:,2]

    xi = linspace(min(xvals), max(xvals), xpts)
    yi = get_yaxis_vals(d, ycol=ycol)
    yi.sort()
    zi = griddata((xvals, yvals), zvals, (xi[None,:], yi[:,None]), method='nearest')
    zi /= zi.max()

    return xi, yi, zi

def plot_grid_data(d, ycol=3):
    xi, yi, zi = make_grid_data(d, ycol=ycol)
    deltay = yi[1] - yi[0]

    fig, ax = subplots(1,1, figsize=(4,4))
    ax.imshow(zi, origin='lower', extent=(min(xi),max(xi),yi.min()-deltay/2.,yi.max()+deltay/2.), 
        aspect='auto', cmap='bwr', interpolation='nearest')

    ax.set_yticks(yi[::4])
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Scan index')

    return fig, ax

def plot_offset_lines(d, ycol=3, of=5000):
    yvals = d[:,ycol]
    xvals = d[:,1]
    zvals = d[:,2]

    yi = get_yaxis_vals(d, ycol=ycol)
    yi.sort()

    fig, ax = subplots(1,1, figsize=(4,8))

    for i,y in enumerate(yi):
        idx = np.where(yvals==yi[i])[0]
        ax.plot(xvals[idx], zvals[idx] + i*of, 'b-')

    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Counts (Hz)')

    return fig, ax

FN = r"D:\measuring\data\20130902\170211_LaserScansYellowRepump_LT1_Yellow\170211_LaserScansYellowRepump_LT1_Yellow.dat"

# gate + 500 nW red; red scan
# FN = r'D:\measuring\data\20130901\213043_LaserScansYellowRepump_LT1_Red\213043_LaserScansYellowRepump_LT1_Red.dat'
# yellow scan
# FN = r'D:\measuring\data\20130901\213043_LaserScansYellowRepump_LT1_Yellow\213043_LaserScansYellowRepump_LT1_Yellow.dat'

# normal red data - no additional laser on
# FN = r'D:\measuring\data\20130901\221059_LaserScansYellowRepump_LT1_Red\221059_LaserScansYellowRepump_LT1_Red.dat'
# corresponding yellow
# FN = r'D:\measuring\data\20130901\221059_LaserScansYellowRepump_LT1_Yellow\221059_LaserScansYellowRepump_LT1_Yellow.dat'

# no gate + 500 nW red; red scan
# FN = r'D:\measuring\data\20130902\101632_LaserScansYellowRepump_LT1_Red\101632_LaserScansYellowRepump_LT1_Red.dat'
# yellow
# FN = r'D:\measuring\data\20130902\101632_LaserScansYellowRepump_LT1_Yellow\101632_LaserScansYellowRepump_LT1_Yellow.dat'

d = get_data(FN)

fig, ax = plot_grid_data(d)
fig2,ax2 = plot_offset_lines(d)

# ax.set_xlim(54,69)


