import numpy as np

# name of the group (top-level) in the data that contains all analysis results
ANALYSISGRP = 'analysis'

# where in the stats array of the adwin data is the number of sequence starts?
ADWIN1STATS_SEQSTARTS_IDX = 6

# Color cycle for plotting
# same colors as in matplotlibrc -- handy for barplots, since they don't use the rc colorcycle
COLORS = 'RoyalBlue', 'Crimson', 'LimeGreen', 'DarkOrange', 'BlueViolet'

# settings the correct photon sync time range
TAIL_LENGTH = 50
CH0_START = 92
CH0_STOP = CH0_START + TAIL_LENGTH
CH1_START = 79.6
CH1_STOP = CH1_START + TAIL_LENGTH

# Hardware / sequence settings
ADWINPROCESSCYCLE = 1e-6
ADWINTIMEHISTBIN = 1e-3
PIPULSES = 2
SEQREPS = 1000
PIPULSESEP = 600
CHANNELDELAY = 12.4 # ch0 arrival time minus ch1 arrival time, in ns

# For tail analysis
PHOTONHIST_BINEDGES = np.arange(-0.1,3000.2,0.2)
FIT_TAIL_LENGTH = 300
