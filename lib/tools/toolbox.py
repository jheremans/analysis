# some convenience tools
#

import os
import time
import logging
import numpy as np

try:
    import qt
    datadir = qt.config['datadir']
except:
    datadir = r'd:\measuring\data'

def nearest_idx(array, value):
    '''
    find the index of the value closest to the specified value.
    '''
    return np.abs(array-value).argmin()

def nearest_value(array, value):
    '''
    find the value in the array that is closest to the specified value.
    '''
    return array[nearest_idx(array,value)]

def verify_timestamp(timestamp):
    if len(timestamp) == 6:
        daystamp = time.strftime('%Y%m%d')
        tstamp = timestamp
    elif len(timestamp) == 14:
        daystamp = timestamp[:8]
        tstamp = timestamp[8:]
    else:
        raise Exception("Cannot interpret timestamp '%s'" % timestamp)
        
    return daystamp, tstamp

def is_older(ts0, ts1):
    '''
    returns True if timestamp ts0 is an earlier data than timestamp ts1,
    False otherwise.
    '''
    dstamp0, tstamp0 = verify_timestamp(ts0)
    dstamp1, tstamp1 = verify_timestamp(ts1)

    return (dstamp0+tstamp0) < (dstamp1+tstamp1)

def latest_data(contains='', older_than=None):
    '''
    finds the latest taken data with <contains> in its name.
    returns the full path of the data directory.
    
    if older_than is not None, than the latest data that fits and that 
    is older than the date given by the timestamp older_than is returned.

    If no fitting data is found, an exception is raised.
    '''
    daydirs = os.listdir(datadir)
    if len(daydirs) == 0:
        logging.warning('No data found in datadir')
        return None

    daydirs.sort()
    
    measdirs = []
    i = len(daydirs)-1
    while len(measdirs) == 0 and i >= 0:       
        daydir = daydirs[i]
        all_measdirs = [d for d in os.listdir(os.path.join(datadir, daydir))]
        all_measdirs.sort()
        
        measdirs = []
        for d in all_measdirs:
            
            # this routine verifies that any output directory is a 'valid' directory
            # (i.e, obeys the regular naming convention)
            _timestamp = daydir + d[:6]
            try:
                dstamp,tstamp = verify_timestamp(_timestamp)
            except:
                continue
            timestamp = dstamp+tstamp

            if contains in d:
                if older_than != None:
                    if not is_older(timestamp, older_than):
                        continue
                measdirs.append(d)
        i -= 1
        
    if len(measdirs) == 0:
        raise Exception('No fitting data found.')

    measdirs.sort()
    measdir = measdirs[-1]

    return os.path.join(datadir,daydir,measdir)

# def newer_data(starttimestamp, endtimestamp=None, contains='',):
#     '''
#     finds all matching data in the given range. default end is now.
#     '''
#     if endtimestamp == None:
#         endtimestamp = time.strftime("%Y%m%d%H%M%S")

#     results = []
#     daydirs = os.listdir(datadir)
#     if len(daydirs) == 0:
#         logging.warning('No data found in datadir')
#         return None
#     daydirs.sort()
    
#     if len(endtimestamp) == 6:
#         endday = time.strftime('%Y%m%d')
#         endtime = endtimestamp
#     elif len(endtimestamp) == 14:
#         endday = endtimestamp[:8]
#         endtime = endtimestamp[8:]
#     else:
#         logging.warning("Cannot interpret timestamp '%s'" % endtimestamp)
#         return None

#     if len(starttimestamp) == 6:
#         startday = time.strftime('%Y%m%d')
#         starttime = starttimestamp
#     elif len(starttimestamp) == 14:
#         startday = starttimestamp[:8]
#         starttime = starttimestamp[8:]
#     else:
#         logging.warning("Cannot interpret timestamp '%s'" % starttimestamp)
#         return None

#     # TODO continue here

#     return None

def data_from_time(timestamp):
    '''
    returns the full path of the data specified by its timestamp in the
    form YYYYmmddHHMMSS.
    '''
    daydirs = os.listdir(datadir)
    if len(daydirs) == 0:
        raise Exception('No data in the data directory specified')
    
    daydirs.sort()    
    daystamp, tstamp = verify_timestamp(timestamp)

    if not os.path.isdir(os.path.join(datadir,daystamp)):
        logging.warning("Requested day '%s' not found" % daystamp)
        return None

    measdirs = [ d for d in os.listdir(os.path.join(datadir,daystamp)) \
            if d[:6] == tstamp ]
    
    if len(measdirs) == 0:
        logging.warning("Requested data '%s'/'%s' not found" \
                % (daystamp, tstamp))
        return None
    elif len(measdirs) == 1:
        return os.path.join(datadir,daystamp,measdirs[0])
    else:
        logging.warning('Timestamp is not unique: ', str(measdirs))
        return None

def measurement_filename(directory=os.getcwd(), ext='hdf5'):
    dirname = os.path.split(directory)[1]
    fn = dirname+'.'+ext

    if os.path.exists(os.path.join(directory,fn)):
        return os.path.join(directory,fn)
    else:
        logging.warning("Data path '%s' does not exist" % \
                os.path.join(directory,fn))
        return None

    

