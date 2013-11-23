import numpy as np
import h5py

import events, files, settings

##############################################################################
### sequence starts
##############################################################################

def get_sequence_starts(fp):
    msmt_name = files.get_msmt_name(fp)
    
    f = h5py.File(fp, 'r')
    adwin1stats = f['/{}/adwin_lt1_data/statistics'.format(msmt_name)].value
    f.close()
    
    return adwin1stats[settings.ADWIN1STATS_SEQSTARTS_IDX]

def get_sequence_starts_from_folder(folder):
    filepaths = files.get_all_msmt_filepaths(folder)
    
    starts = 0
    for f in filepaths:
        start += get_sequence_starts(f)

    return starts


##############################################################################
### CR stats
##############################################################################

def _get_CR_hist(fp, adwin, data_name, **kw):
    normalize = kw.pop('normalize', True)

    msmt_name = files.get_msmt_name(fp)
    folder_name = 'adwin_lt1_data' if adwin=='lt1' else 'adwin_lt2_data'
    
    f = h5py.File(fp, 'r')    
    try:
        h = f['{}/{}/{}'.format(msmt_name, folder_name, data_name)].value
        f.close()
        
    except:
        # this typically comes from older, buggy msmts. create some dummy data then.
        print "Cannot find data '{}/{}/{}' in {}".format(msmt_name, folder_name, data_name, fp)
        
        f.close()
        return np.array([0]), np.array([0])
    
    # only return the non-zero part, make sure older, buggy msmts don't break
    nonzeros = np.where(h>0)[0]
    if len(nonzeros) > 1:
        maxbin = nonzeros[-1]
    else:
        maxbin = 0

    div = float(h.sum()) if normalize else 1.
    return h[:maxbin]/div, np.arange(maxbin, dtype='int')

def get_CR_hist_total(fp, adwin, **kw):
    data_name = 'CR_hist_all' if adwin=='lt1' else 'CR_hist'
    return _get_CR_hist(fp, adwin, data_name, **kw)

def get_CR_hist_sequence_timeout(fp, adwin, **kw):
    data_name = 'CR_hist_time_out'
    return _get_CR_hist(fp, adwin, data_name, **kw)

def _make_CR_hist(fp, adwin, data_name, **kw):
    normalize = kw.pop('normalize', True)

    msmt_name = files.get_msmt_name(fp)
    folder_name = 'adwin_lt1_data' if adwin=='lt1' else 'adwin_lt2_data'

    f = h5py.File(fp, 'r')    
    vals = f['{}/{}/{}'.format(msmt_name, folder_name, data_name)].value
    f.close()

    binedges = np.arange(-0.5, vals.max()+0.5)
    h,b = np.histogram(vals, bins=binedges, density=True)

    return h,b

def get_CR_hist_before(fp, adwin, **kw):
    data_name = 'CR_before'
    return _make_CR_hist(fp, adwin, data_name, **kw)

def get_CR_hist_after(fp, adwin, **kw):
    data_name = 'CR_after'
    return _make_CR_hist(fp, adwin, data_name, **kw)





