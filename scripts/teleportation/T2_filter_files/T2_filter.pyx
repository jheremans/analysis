"""
Cython-based fast data analysis for HH-obtained event data.

authors:
Wolfgang Pfaff
"""

import cython
import time
import numpy as np
from cython.view cimport array as cvarray
cimport numpy as cnp

@cython.boundscheck(False)
@cython.wraparound(False)
def hist2d(cnp.ndarray[cnp.uint8_t, ndim=1, mode='c'] channel not None,
    cnp.ndarray[cnp.uint8_t, ndim=1, mode='c'] special not None,
    cnp.ndarray[cnp.uint64_t, ndim=1, mode='c'] sync_time not None,
    period,
    pts,
    t_start,
    t_stop,
    t_bins,
    hist_channel = 0):
    
    """
    Get a 2D histogram of photon counts, with axes sync time and some second parameter.
    We assume here that we get <period> markers in row for the same sweep parameter, and cycle
    through <pts> points of that parameter.
    The photon events are binned according to <t_start>, <t_stop>, and <t_bins>, where
    we only look at the photons on channel <hist_channel>.

    The function returns the <pts> x <t_bins> histogram, and the values of the bin edges.
    """

    cdef cnp.uint64_t i
    cdef cnp.uint64_t k
    cdef cnp.uint32_t period_cnt = 0
    cdef cnp.uint32_t pt = 0
    cdef cnp.uint32_t bin = 0
    cdef cnp.uint64_t length = channel.shape[0]
    cdef cnp.ndarray[cnp.uint32_t, ndim=2, mode='c'] hist = np.zeros((pts, t_bins), dtype='u4')
    cdef cnp.ndarray[cnp.uint64_t, ndim=1, mode='c'] binedges = np.linspace(t_start, t_stop, t_bins+1).astype('u8')

    for i in range(length):        
        if special[i] == 1 and channel[i] == 1:
            period_cnt += 1
            if period_cnt >= period:
                period_cnt = 0
                pt += 1
                if pt >= pts:
                    pt = 0
            continue

        elif special[i] == 0 and channel[i] == hist_channel:

            # find the correct bin
            if sync_time[i] < binedges[0]:
                continue

            for bin in range(t_bins):
                if binedges[bin+1] > sync_time[i]:
                    hist[pt, bin] += 1
                    break

    return hist, binedges


