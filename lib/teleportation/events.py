import numpy as np
import h5py

from analysis.lib.fitting import fit
from analysis.lib import tools

import files, stats, settings

# columns: sync no | attempt | ph1 arrival time | ph1 ch | ph2 arrival time | ph2 ch | dt | BSM outcome | Bob outcome |
#          LDE_type (1 = psi-, 0 = psi+) | CR before 1 | CR after 1 | CR before 2 | CR after 2 | PLU mrkr abs. time

T_EV_COLS = 15
T_EV_COL_SYNCNUMBER = 0
T_EV_COL_ATTEMPT = 1
T_EV_COL_PH1_TIME = 2
T_EV_COL_PH1_CHAN = 3
T_EV_COL_PH2_TIME = 4
T_EV_COL_PH2_CHAN = 5
T_EV_COL_DT = 6
T_EV_COL_BSM_RESULT = 7
T_EV_COL_BOB_RESULT = 8
T_EV_COL_PSIMINUS = 9
T_EV_COL_CR1_BEFORE = 10
T_EV_COL_CR1_AFTER = 11
T_EV_COL_CR2_BEFORE = 12
T_EV_COL_CR2_AFTER = 13
T_EV_COL_PLU_ABS_TIME = 14


### identification of events
def get_photons(fp):
    """
    returns two filters (1d-arrays): whether events are ch0-photons/ch1-photons
    """
    f = h5py.File(fp, 'r')
    channel = f['/HH_channel-1'].value
    special = f['/HH_special-1'].value
    f.close()

    is_not_special = special==0
    is_channel_0 = channel==0
    is_channel_1 = channel==1

    is_photon_0 = np.logical_and(is_not_special, is_channel_0)
    is_photon_1 = np.logical_and(is_not_special, is_channel_1)

    return is_photon_0, is_photon_1

def get_markers(fp, chan):
    """
    returns a filter (1d-array): whether events are markers on the given channel
    """
    f = h5py.File(fp, 'r')
    channel = f['/HH_channel-1'].value
    special = f['/HH_special-1'].value
    f.close()
    
    is_special = special==1
    is_channel = channel==chan
    
    return (is_special & is_channel)

def get_coincidences(fp):
    f = h5py.File(fp, 'r')
    sync_time_ns = f['/HH_sync_time-1'].value * 1e-3
    special = f['/HH_special-1'].value
    sync_numbers = f['/HH_sync_number-1'].value
    f.close()

    is_photon = special == 0
    photon_sync_numbers = sync_numbers[is_photon]
    
    # this works nicely on sorted arrays
    multiple_photon_sync_numbers = np.unique(photon_sync_numbers[1:] == \
                                             photon_sync_numbers[:-1])
    dts = np.array([])

    return


### photon histograms
def get_photon_hist(fp, **kw):
    save = kw.pop('save', False)
    fltr = kw.pop('fltr', None)
    force_eval = kw.pop('force_eval', True)
    binedges = kw.pop('binedges', settings.PHOTONHIST_BINEDGES)
    
    if not force_eval and files.has_analysis_data(fp, 'photon_histogram'):
        h, h_attrs = files.get_analysis_data(fp, 'photon_histogram')
        be, be_attrs = files.get_analysis_data(fp, 'photon_histogram_binedges_ns')
        h0 = h[:,0]
        h1 = h[:,1]
        return (h0, be), (h1, be)
    
    f = h5py.File(fp, 'r')
    sync_time_ns = f['/HH_sync_time-1'].value * 1e-3
    f.close()
    
    ph0, ph1 = get_photons(fp)
    if fltr != None:
        _fltr0 = (ph0 & fltr)
        _fltr1 = (ph1 & fltr)
    else:
        _fltr0 = ph0
        _fltr1 = ph1
    
    st0 = sync_time_ns[_fltr0]
    st1 = sync_time_ns[_fltr1]
    
    h0, b0 = np.histogram(st0, bins=binedges)
    h1, b1 = np.histogram(st1, bins=binedges)
    
    if save:
        files.set_analysis_data(fp, 'photon_histogram', vstack((h0,h1)).transpose(),
                          columns=('channel_0', 'channel_1'))
        files.set_analysis_data(fp, 'photon_histogram_binedges_ns', b0)
        files.delete_analysis_data(fp, 'photon_histogram_event_filter')
        if fltr != None:
            files.set_analysis_data(fp, 'photon_histogram_event_filter', fltr)
        
    return (h0, b0), (h1, b1)


def get_photon_hists_from_folder(folder, **kw):
    '''
    return the cumulative photon histogram from all data contained in a folder
    (all sub-levels are searched).
    '''
    filepaths = files.get_all_msmt_filepaths(folder)
    for i,f in enumerate(filepaths):
        if i == 0:
            (h0,b0),(h1,b1) = get_photon_hist(f, **kw)
        else:
            (_h0,_b0),(_h1,_b1) = get_photon_hist(f, **kw)
            h0 += _h0
            h1 += _h1
    return (h0, b0), (h1, b1)

##############################################################################
### filtering
##############################################################################

def filter_synctime(sync_times, t0, t1):
    """
    return a filter that's true for all sync times that are in between t0 and t1
    (t0 is the lower bound, t1 is the upper bound)
    """
    fltr = ((sync_times >= t0) & (sync_times <= t1))
    return fltr

def filter_on_same_sync_number(source_sync_numbers, target_sync_numbers):
    """
    returns a filter for target_sync_numbers that's true for all sync numbers that are also
    in source_sync_numbers.
    """
    return np.in1d(target_sync_numbers, source_sync_numbers)

def filter_marker(fp, chan):
    """
    Note: at the moment this filter includes the marker events on which we filter.
    """
    is_mrkr = get_markers(fp, chan)
        
    f = h5py.File(fp, 'r')
    sync_numbers = f['/HH_sync_number-1'].value
    f.close()
    
    marker_sync_numbers = sync_numbers[is_mrkr]
    
    return filter_on_same_sync_number(marker_sync_numbers, sync_numbers)


##############################################################################
### Tail analysis
##############################################################################

def _fit_tail(h, b):
    binwidth = b[1] - b[0]
    _b = b-b[0]
    
    A = fit.Parameter(h[0], 'A')
    T = fit.Parameter(11, 'T')
    o = fit.Parameter(0, 'o')
    p0 = [A, T, o]
    fitfunc_str = 'A exp(-t/T) + o'
    
    def ff(x):
        return o() + A() * np.exp(-x/T())
    fitres = fit.fit1d(_b, h, None, p0=p0, fitfunc=ff, fitfunc_str=fitfunc_str,
                    do_print=True, ret=True)
    return fitres

def fit_tail(fp, **kw):
    (h0, b0), (h1, b1) = get_photon_hist(fp, **kw)
        
    
    _x0 = (b0[:-1]+b0[1:])/2.
    _x1 = (b1[:-1]+b1[1:])/2.
    binwidth = _x0[1] - _x0[0]
    tail0 = (_x0 >= settings.CH0_START) & (_x0 <= settings.CH0_START + settings.FIT_TAIL_LENGTH)
    tail1 = (_x1 >= settings.CH1_START) & (_x1 <= settings.CH1_START + settings.FIT_TAIL_LENGTH)
    x0, y0 = _x0[tail0], h0[tail0]
    x1, y1 = _x1[tail1], h1[tail1]
    
    fitres0 = _fit_tail(y0, x0)
    fitres1 = _fit_tail(y1, x1)
    
    starts = stats.get_sequence_starts(fp) * settings.SEQREPS
    tcpsh0 = fitres0['params_dict']['A'] * fitres0['params_dict']['T'] / binwidth / starts
    tcpsh1 = fitres1['params_dict']['A'] * fitres1['params_dict']['T'] / binwidth / starts
    
    _x0 = np.linspace(x0[0], x0[-1], 201)
    ret0 = { 'params' : fitres0['params_dict'], 
            'errors' : fitres0['error_dict'], 
            'data_x' : x0,
            'data_y' : y0, 
            'fit_x' : _x0,
            'fit_y' : fitres0['fitfunc'](_x0-_x0[0]),
            'info' : fitres0['fitfunc_str'] + '\n' + fit.str_fit_params(fitres0),
            'tcpsh' : tcpsh0,
            }
    
    _x1 = np.linspace(x1[0], x1[-1], 201)
    ret1 = { 'params' : fitres1['params_dict'], 
            'errors' : fitres1['error_dict'], 
            'data_x' : x1,
            'data_y' : y1,
            'fit_x' : _x1,
            'fit_y' : fitres1['fitfunc'](_x1-_x1[0]),
            'info' : fitres1['fitfunc_str'] + '\n' + fit.str_fit_params(fitres1),
            'tcpsh' : tcpsh1,
            }
    
    keys0 = []
    params0 = np.empty((0,2))
    for k in ret0['params']:
        keys0.append(k)
        params0 = np.vstack((params0, np.array([ret0['params'][k], ret0['errors'][k]])))
    
    files.set_analysis_data(fp, 'ch0_fitparams', params0, subgroup='tail',
                      info=ret0['info'], columns=('fitvalue', 'error'),
                      rows=keys0)
    files.set_analysis_data(fp, 'ch0_fit_data', np.vstack((ret0['data_x'], ret0['data_y'])).transpose(),
                      subgroup='tail', columns=('time bin center (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch0_fit_curve', np.vstack((ret0['fit_x'], ret0['fit_y'])).transpose(),
                      subgroup='tail', columns=('time (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch0_tcpsh', ret0['tcpsh'], subgroup='tail',
                      starts = starts)
    
    keys1 = []
    params1 = np.empty((0,2))
    for k in ret1['params']:
        keys1.append(k)
        params1 = np.vstack((params1, np.array([ret1['params'][k], ret1['errors'][k]])))
    
    files.set_analysis_data(fp, 'ch1_fitparams', params1, subgroup='tail',
                      info=ret1['info'], columns=('fitvalue', 'error'),
                      rows=keys1)
    files.set_analysis_data(fp, 'ch1_fit_data', np.vstack((ret1['data_x'], ret1['data_y'])).transpose(),
                      subgroup='tail', columns=('time bin center (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch1_fit_curve', np.vstack((ret1['fit_x'], ret1['fit_y'])).transpose(),
                      subgroup='tail', columns=('time (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch1_tcpsh', ret1['tcpsh'], subgroup='tail',
                      starts = starts)
    
    return ret0, ret1

##############################################################################
### Teleportation events
##############################################################################

def get_teleportation_events(folder, force_eval=False, verbose=True):
    fp = files.get_msmt_fp(folder)
    fname = files.get_msmt_name(fp)

    if files.has_analysis_data(fp, 'teleportation_events') and not force_eval:
        tev, _a = files.get_analysis_data(fp, 'teleportation_events')
        return tev

    f = h5py.File(fp, 'r')
    sync_times = f['/HH_sync_time-1'].value * 1e-3 # we prefer ns over ps
    sync_numbers = f['/HH_sync_number-1'].value
    channels = f['/HH_channel-1'].value
    abs_times = f['/HH_time-1'].value
    ad1_reps = f['{}/adwin_lt1_data/completed_reps'.format(fname)].value
    ad1_ssro1 = f['{}/adwin_lt1_data/SSRO1_results'.format(fname)].value
    ad1_ssro2 = f['{}/adwin_lt1_data/SSRO2_results'.format(fname)].value
    ad2_reps = f['{}/adwin_lt2_data/completed_reps'.format(fname)].value
    ad2_ssro = f['{}/adwin_lt2_data/SSRO_lt2_data'.format(fname)].value
    ad1_CR_before = f['{}/adwin_lt1_data/CR_before'.format(fname)].value
    ad2_CR_before = f['{}/adwin_lt2_data/CR_before'.format(fname)].value
    ad1_CR_after = f['{}/adwin_lt1_data/CR_after'.format(fname)].value
    ad2_CR_after = f['{}/adwin_lt2_data/CR_after'.format(fname)].value
    
    seq_reps = int(f['{}'.format(fname)].attrs['LDE_attempts_before_CR'])

    f.close()

    # binary encoding: 0 = 00, 1 = 01, 2 = 10, 3 = 11
    BSM_outcomes = (2**1 * ad1_ssro2) + (2**0 * ad1_ssro1)

    ### First the HH events
    PLU_mrkr = 2
    psiminus_mrkr = 3
    psiplus_mrkr = 4
    # this conversion might be neccesary due to a wrong conversion in T2_tools.pyx
    # of the channel bits from a uint_32 to a uint_8. 
    is_PLU_mrkr = get_markers(fp, 2**(PLU_mrkr-1))
    has_same_sync_number_as_PLU_mrkr = filter_marker(fp, 2**(PLU_mrkr-1))
    is_psiplus_mrkr = get_markers(fp, 2**(psiplus_mrkr-1))
    is_psiminus_mrkr = get_markers(fp, 2**(psiminus_mrkr-1))
    has_same_sync_number_as_psiminus_mrkr = filter_marker(fp, 2**(psiminus_mrkr-1))

    is_ph_ch0, is_ph_ch1 = get_photons(fp)
    is_ph = is_ph_ch0 | is_ph_ch1

    is_ph_with_PLU_mrkr = is_ph & filter_marker(fp, 2**(PLU_mrkr-1))
    is_ph0_with_PLU_mrkr = is_ph_ch0 & filter_marker(fp, 2**(PLU_mrkr-1))
    is_ph1_with_PLU_mrkr = is_ph_ch1 & filter_marker(fp, 2**(PLU_mrkr-1))
    is_ph0_with_psiminus_mrkr = is_ph_ch0 & filter_marker(fp, 2**(psiminus_mrkr-1))
    is_ph1_with_psiminus_mrkr = is_ph_ch1 & filter_marker(fp, 2**(psiminus_mrkr-1))

    if settings.VERBOSE:

        print
        print fp
        print
        print 'HH data'
        print '-------'
        print 'Number of PLU markers received:', len(np.where(is_PLU_mrkr)[0])
        print 'Number of Psi+ markers received:', len(np.where(is_psiplus_mrkr)[0])
        print 'Number of Psi- markers received:', len(np.where(is_psiminus_mrkr)[0])
        print 'Number of photons with the same sync as a PLU marker:', len(np.where(is_ph_with_PLU_mrkr)[0])
        print '    ch0:', len(np.where(is_ph0_with_PLU_mrkr)[0])
        print '    ch1:', len(np.where(is_ph1_with_PLU_mrkr)[0])

    is_ph0_with_PLU_mrkr_in_first_window = is_ph0_with_PLU_mrkr & \
        filter_synctime(sync_times, settings.CH0_START,settings.CH0_STOP)
    is_ph1_with_PLU_mrkr_in_first_window = is_ph1_with_PLU_mrkr & \
        filter_synctime(sync_times, settings.CH1_START, settings.CH1_STOP)
    is_ph_with_PLU_mrkr_in_first_window = is_ph0_with_PLU_mrkr_in_first_window | \
        is_ph1_with_PLU_mrkr_in_first_window

    if settings.VERBOSE:

        print
        print 'Number of PLU-marked photons in first window:', \
            len(np.where(is_ph_with_PLU_mrkr_in_first_window)[0])

    is_ph0_with_PLU_mrkr_in_second_window = is_ph0_with_PLU_mrkr & filter_synctime(sync_times, 
        settings.CH0_START+settings.PIPULSESEP, 
        settings.CH0_STOP+settings.PIPULSESEP)
    is_ph1_with_PLU_mrkr_in_second_window = is_ph1_with_PLU_mrkr & filter_synctime(sync_times, 
        settings.CH1_START+settings.PIPULSESEP, 
        settings.CH1_STOP+settings.PIPULSESEP)
    is_ph_with_PLU_mrkr_in_second_window = is_ph0_with_PLU_mrkr_in_second_window | \
        is_ph1_with_PLU_mrkr_in_second_window

    if settings.VERBOSE:

        print 'Number of PLU-marked photons in second window:', \
            len(np.where(is_ph_with_PLU_mrkr_in_second_window)[0])

    LDE_ph0_sync_times = sync_times[is_ph0_with_PLU_mrkr]
    LDE_ph0_sync_numbers = sync_numbers[is_ph0_with_PLU_mrkr]
    LDE_ph1_sync_times = sync_times[is_ph1_with_PLU_mrkr]
    LDE_ph1_sync_numbers = sync_numbers[is_ph1_with_PLU_mrkr]
    PLU_mrkr_sync_numbers = sync_numbers[is_PLU_mrkr]
    PLU_mrkr_abs_times = abs_times[is_PLU_mrkr]
    psi_minus_mrkr_sync_numbers = sync_numbers[is_psiminus_mrkr]
    psi_plus_mrkr_sync_numbers = sync_numbers[is_psiplus_mrkr]

    attempts = np.zeros(PLU_mrkr_sync_numbers.size)
    of = 0
    for i, sn in enumerate(PLU_mrkr_sync_numbers):
        a = (sn-of) % seq_reps
        of = sn
        attempts[i] = a

    ### Adwin SSROs
    if settings.VERBOSE:

        print 
        print 'Adwin LT1'
        print '---------'
        print 'Number of events:', ad1_reps,
        if ad1_reps != len(np.where(is_PLU_mrkr)[0]):
            raise Exception(fp + '\nAdwin LT1 events do not match HH events - data set seems faulty :(')
        else:
            print 'OK :)'
            
        print 
        print 'Adwin LT2'
        print '---------'
        print 'Number of events:', ad2_reps,
        if ad2_reps != len(np.where(is_PLU_mrkr)[0]):
            raise Exception(fp+'\nAdwin LT2 events do not match HH events - data set seems faulty :(')
        else:
            print 'OK :)'
        
    # columns: sync no | attempt | ph1 arrival time | ph1 ch | ph2 arrival time | ph2 ch | dt | BSM outcome | Bob outcome |
    #          LDE_type (1 = psi-, 0 = psi+) | CR before 1 | CR after 1 | CR before 2 | CR after 2 | PLU mrkr abs. time
    teleportation_events = np.empty((0,15))

    for i,s in enumerate(PLU_mrkr_sync_numbers):
        attempt = attempts[i]

        _ph0 = LDE_ph0_sync_numbers==s
        _ph1 = LDE_ph1_sync_numbers==s
        
        if len(LDE_ph0_sync_numbers[_ph0]) == 1 and len(LDE_ph1_sync_numbers[_ph1]) == 1:
            psiminus = 1
            stimes = np.array([ LDE_ph0_sync_times[_ph0], LDE_ph1_sync_times[_ph1] ]).reshape(-1)
            chans = np.array([0,1])            

        elif len(LDE_ph0_sync_numbers[_ph0]) == 2 and len(LDE_ph1_sync_numbers[_ph1]) == 0:
            psiminus = 0
            stimes = LDE_ph0_sync_times[_ph0].reshape(-1)
            chans = np.array([0,0])

        elif len(LDE_ph0_sync_numbers[_ph0]) == 0 and len(LDE_ph1_sync_numbers[_ph1]) == 2:
            psiminus = 0
            stimes = LDE_ph1_sync_times[_ph1].reshape(-1)
            chans = np.array([1,1])

        else:
            continue

        
        idx = np.argsort(stimes)
        atimes = stimes[idx]
        atimes[1] -= settings.PIPULSESEP
        chans = chans[idx]
        dt = atimes[1] - atimes[0]
       
        _evnt = np.array([s, 
            attempt,
            atimes[0],
            chans[0],
            atimes[1],
            chans[1],
            dt,
            BSM_outcomes[i],
            ad2_ssro[i],
            psiminus,
            ad1_CR_before[i],
            ad1_CR_after[i],
            ad2_CR_before[i],
            ad2_CR_after[i],
            PLU_mrkr_abs_times[i],
            ])

        teleportation_events = np.vstack((teleportation_events, _evnt))

    if settings.VERBOSE:
        print
        print 'Found {} valid teleportation events.'.format(int(len(teleportation_events)))
        print '===================================='
        print 

    files.set_analysis_data(fp, 'teleportation_events', teleportation_events)

    return teleportation_events

def filter_event_times(teleportation_events, **kw):
    ch0_start = kw.pop('ch0_start', settings.CH0_START)
    ch1_start = kw.pop('ch1_start', settings.CH1_START)
    window_length = kw.pop('window_length', settings.WINDOW_LENGTH)
    dt_max = kw.pop('dt_max', settings.DT_MAX)
    second_window_length = kw.pop('second_window_length', window_length)

    ch0_stop = ch0_start + window_length
    ch1_stop = ch1_start + window_length

    ch0_stop_w2 = ch0_start + second_window_length
    ch1_stop_w2 = ch1_start + second_window_length

    tev = teleportation_events
    win1_ch0_valid = (tev[:,T_EV_COL_PH1_CHAN] == 0) & \
        (tev[:,T_EV_COL_PH1_TIME] >= ch0_start) & \
        (tev[:,T_EV_COL_PH1_TIME] <= ch0_stop)

    win1_ch1_valid = (tev[:,T_EV_COL_PH1_CHAN] == 1) & \
        (tev[:,T_EV_COL_PH1_TIME] >= ch1_start) & \
        (tev[:,T_EV_COL_PH1_TIME] <= ch1_stop)

    win2_ch0_valid = (tev[:,T_EV_COL_PH2_CHAN] == 0) & \
        (tev[:,T_EV_COL_PH2_TIME] >= ch0_start) & \
        (tev[:,T_EV_COL_PH2_TIME] <= ch0_stop_w2)

    win2_ch1_valid = (tev[:,T_EV_COL_PH2_CHAN] == 1) & \
        (tev[:,T_EV_COL_PH2_TIME] >= ch1_start) & \
        (tev[:,T_EV_COL_PH2_TIME] <= ch1_stop_w2)

    wins_valid = (win1_ch0_valid | win1_ch1_valid) & (win2_ch0_valid | win2_ch1_valid)

    dt_valid = np.abs(tev[:,T_EV_COL_DT]) < dt_max

    is_valid = (wins_valid & dt_valid)
    return is_valid


def filter_attempt_number(teleportation_events, **kw):
    max_attempt_number = kw.pop('max_attempt_number', settings.MAX_ATTEMPT_NUMBER)
    return teleportation_events[:,T_EV_COL_ATTEMPT] <= max_attempt_number

def conditional_results(teleportation_events):

    # note that for the BSM outcome the first qubit is the nitrogen
    BSM_outcome_names = ['00', '01', '10', '11']
    BSM_outcomes = teleportation_events[:,T_EV_COL_BSM_RESULT]
    Bob_outcomes = teleportation_events[:,T_EV_COL_BOB_RESULT]
    is_psiminus_event = teleportation_events[:,T_EV_COL_PSIMINUS] == 1
    is_psiplus_event = teleportation_events[:,T_EV_COL_PSIMINUS] == 0

    Bob_outcomes_psiminus = np.empty((0,2))
    Bob_outcomes_psiplus = np.empty((0,2))
    BSM_outcomes_psiminus = np.array([])
    BSM_outcomes_psiplus = np.array([])

    for outcome in range(4):
        ssro_results_psiminus = Bob_outcomes[((BSM_outcomes==outcome) & is_psiminus_event)]   
        BSM_outcomes_psiminus = np.append(BSM_outcomes_psiminus,
            len(ssro_results_psiminus))
        Bob_outcomes_psiminus = np.vstack((Bob_outcomes_psiminus,
            np.array([len(ssro_results_psiminus[ssro_results_psiminus==0]),
                      len(ssro_results_psiminus[ssro_results_psiminus==1])])))

        ssro_results_psiplus = Bob_outcomes[((BSM_outcomes==outcome) & is_psiplus_event)]   
        BSM_outcomes_psiplus = np.append(BSM_outcomes_psiplus,
            len(ssro_results_psiplus))
        Bob_outcomes_psiplus = np.vstack((Bob_outcomes_psiplus,
            np.array([len(ssro_results_psiplus[ssro_results_psiplus==0]),
                      len(ssro_results_psiplus[ssro_results_psiplus==1])])))

    return BSM_outcomes_psiminus, BSM_outcomes_psiplus, \
        Bob_outcomes_psiminus, Bob_outcomes_psiplus
