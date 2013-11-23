import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm


x = np.arange(100)/10.
y = np.arange(100)/10.
xx,yy = np.meshgrid(x,y)
zz = xx**2 + yy**2

fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(11,8))
im1 = ax1.imshow(zz, cmap=cm.gist_earth)
im2 = ax2.imshow(zz*2, cmap=cm.gist_earth)
im3 = ax3.imshow(zz*3, cmap=cm.gist_earth)
im4 = ax4.imshow(zz*4, cmap=cm.gist_earth)

cb1 = fig.colorbar(im1, ax=ax1)
cb2 = fig.colorbar(im2, ax=ax2)
cb3 = fig.colorbar(im3, ax=ax3)
cb4 = fig.colorbar(im4, ax=ax4)




