import h5py
import os

import settings

### data finding and identification
def get_all_msmt_filepaths(folder, suffix='hdf5', pattern=''):
    filepaths = []
    suffixlen = len(suffix)
    
    for root,dirs,files in os.walk(folder):
        for f in files:
            if len(f) > suffixlen and f[-suffixlen:] == suffix and pattern in f:
                filepaths.append(os.path.join(root, f))
    
    return filepaths

def get_msmt_name(fp):
    """
    This assumes that there is only one group, whose name is the msmt name.
    """
    _root, fn = os.path.split(fp)
    f = h5py.File(fp, 'r')
    for k in f.keys():
        if f.get(k, getclass=True) == h5py._hl.group.Group and k in fn:
            f.close()
            return k
    
    raise Exception('Cannot find the name of the measurement.')
    
def get_msmt_fp(folder, ext='hdf5'):
    dirname = os.path.split(folder)[1]
    fn = dirname+'.'+ext
    return os.path.join(folder, fn)

def get_msmt_header(fp):
    _root, fn = os.path.split(fp)
    root, folder = os.path.split(_root)
    daystamp = os.path.split(root)[1]
    return daystamp + '/' + folder


### basic functionality for analysis data
def delete_analysis_data(fp, name, subgroup=None):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    try:
        agrp = f.require_group(settings.ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))
        if name in agrp.keys():
            del agrp[name]
        f.flush()
        f.close()        
    except:
        f.close()
        raise

def clear_analysis_data(fp):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    if settings.ANALYSISGRP in f.keys():
        try:
            del f[settings.ANALYSISGRP]
            f.flush()
            f.close()
        except:
            f.close()
            raise
    else:
        f.close()

def set_analysis_data(fp, name, data, subgroup=None, **kw):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise
    
    try:
        agrp = f.require_group(settings.ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))
        if name in agrp.keys():
            del agrp[name]
        agrp[name] = data
        f.flush()
        
        for k in kw:
            agrp[name].attrs[k] = kw[k]
        
        f.flush()
        f.close()        
    except:
        f.close()
        raise
        
def has_analysis_data(fp, name, subgroup=None):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise
    
    agrp = f.require_group(settings.ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))
    if name in agrp.keys():
        f.close()
        return True
    else:
        f.close()
        return False
    
def get_analysis_data(fp, name, subgroup=None):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise
        
    agrp = f.require_group(settings.ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))

    if name not in agrp.keys():
        return None
    
    dat = agrp[name].value
    attrs = {}
    for (an, av) in agrp[name].attrs.items():
        attrs[an] = av
        
    f.close()
        
    return dat, attrs
