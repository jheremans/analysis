import os
import numpy as np
import matplotlib.pyplot as plt

from analysis.lib.pq import hht3

setup = 'lt1'

DATFILEBASE = 'LDE'
DATIDX = 0
DATIDXDIGITS = 3
ADWINLT1FOLDER = 'adwin_lt1'
ADWINLT2FOLDER = 'adwin_lt2'
ADWINPARAMSFN = 'lde_params.npz'
ADWINDATAFN = 'lde_data.npz'

def set_setup(s):
    global setup
    global ADWINFOLDER
    
    setup = s
    ADWINFOLDER = ADWINLT2FOLDER if setup == 'lt2' else ADWINLT1FOLDER

    return

def correct_ssro_mismatch(SSROs_adwin_LT1, SSROs_adwin_LT2, noof_PLU_markers,
        noof_SSROs_adwin_LT1, noof_SSROs_adwin_LT2, CR_adwin_LT1=None):
    added_idxs = 0
    bad_idx = list() #index where 
    for i in range(noof_PLU_markers):
        if (SSROs_adwin_LT2[i]==-1):
            added_idxs += 1
            bad_idx.append(i)
            SSROs_adwin_LT1 = np.insert(SSROs_adwin_LT1,i,-1)
            if CR_adwin_LT1 != None:
                CR_adwin_LT1 = np.insert(CR_adwin_LT1,i,-1)
    
    # print 'Detected missed PLU markers at SSRO LT2 array for indexes: ',bad_idx
    # print 'Compensating SSRO mismatch:'
    # print '\tADwin LT2 received PLU markers: ',noof_PLU_markers
    # print '\tADwin LT2 performed SSROs: ',noof_SSROs_adwin_LT2
    # print '\tADwin LT1 performed SSROs: ',noof_SSROs_adwin_LT1
    # print 'Finished adding -1 to SSRO LT1 data. Added %s element.'%(added_idxs)
    
    if CR_adwin_LT1 == None:
        return SSROs_adwin_LT1, added_idxs
    return SSROs_adwin_LT1, CR_adwin_LT1, added_idxs

def prepare(datfolder, rawfolder, *arg, **kw):
    global setup

    #############################
    #### LOAD ADWIN LT1 DATA ####
    #############################

    # load adwin autodata 
    fn = os.path.join(datfolder, 
            ADWINLT1FOLDER+('-%.'+str(DATIDXDIGITS)+'d') % DATIDX, 
            ADWINDATAFN)
    d = np.load(fn)
    
    noof_SSROs_adwin_LT1 = d['get_noof_SSROs']
    noof_SSRO_triggers_LT1 = d['get_noof_SSRO_triggers']
    d.close()


    #############################
    #### LOAD ADWIN LT2 DATA ####
    #############################
    
    sn = os.path.join(datfolder, 
        ADWINLT2FOLDER+('-%.'+str(DATIDXDIGITS)+'d') % DATIDX, 
        ADWINDATAFN)
    s = np.load(sn)

    noof_PLU_markers = s['get_noof_PLU_markers']
    noof_SSROs_adwin_LT2 = s['get_noof_SSROs']
    s.close()

    ############################################
    ### LOADING SSRO DATA###
    ############################################ 

    fn = os.path.join(datfolder,
         DATFILEBASE+('-%.'+str(DATIDXDIGITS)+'d.npz') % DATIDX)
    d = np.load(fn)
    
    SSROs_adwin_LT1 = d['adwin_lt1_SSRO']
    SSROs_adwin_LT2 = d['adwin_lt2_SSRO']

    ############################################
    ### COMPENSATING A MISMATCH IN SSRO DATA ###
    ############################################ 
    #filling up SSROs_adwin_lt1 with -1 where LT2 missed a PLU

    SSROs_adwin_LT1, added_idxs = correct_ssro_mismatch(SSROs_adwin_LT1, 
            SSROs_adwin_LT2, noof_PLU_markers,noof_SSROs_adwin_LT1, 
            noof_SSROs_adwin_LT2)
    
    ######################################
    ######################################
    ######################################    
    
    if setup == 'lt2':
        CRs_adwin = d['adwin_lt2_CR']
    else:
        CRs_adwin = d['adwin_lt1_CR']
    d.close()
    
    # load HH data
    hh_prefiltered, window1, window2 = split(rawfolder, *arg, **kw)
    hhpludat = hht3.filter_counts_on_marker(
            hh_prefiltered, 2, [-1,0])
    
    np.savez(os.path.join(datfolder, 'hh_prefiltered_events'),
            data=hh_prefiltered)
    np.savez(os.path.join(datfolder, 'window1_events'),
            data=window1)
    np.savez(os.path.join(datfolder, 'window2_events'),
            data=window2)
    
    # more preparations
    SSROs_adwin_LT1 = SSROs_adwin_LT1[:noof_SSRO_triggers_LT1+added_idxs]
    SSROs_adwin_LT2 = SSROs_adwin_LT2[:noof_PLU_markers]
    if setup == 'lt1':
        CRs_adwin = CRs_adwin[:noof_SSRO_triggers_LT1]
    
    ######################################
    ######################################
    ######################################    
    
    print 'PLU markers on adwin lt2:', noof_PLU_markers
    
    print '\nLength of ADwin LT1 SSRO array after addition: '\
            , len(SSROs_adwin_LT1)
    print 'Length of ADwin LT2 SSRO array: ', len(SSROs_adwin_LT2)

    ######################################
    ######################################
    ######################################    

    ### first some sanity checks
    invalid_SSROs_adwin = np.where(SSROs_adwin_LT1==-1)[0] \
            if setup =='lt1' else np.where(SSROs_adwin_LT2==-1)[0]
    
    print ''
    print 'adwin detected SSROs: %d' % \
        (noof_SSRO_triggers_LT1 if setup=='lt1' else noof_PLU_markers)
    
    print '  thereof invalid: %d' % len(invalid_SSROs_adwin)
    
    print 'HH SSROs: %d' % len(hhpludat[
        np.logical_and(hhpludat[:,3]==1, hhpludat[:,2]==2)])
    print ''

    ### plot counts
    if(kw.get('do_plot',False)):
        plot_counts_vs_time(window1, datfolder, window=1)
        plot_counts_vs_time(window2, datfolder, window=2)

    ret = SSROs_adwin_LT1 if setup=='lt1' else SSROs_adwin_LT2
    return ret

def prepare_correlations(datfolder, datfile='hh_prefiltered_events.npz'):
    
    # HH data
    d = np.load(os.path.join(datfolder, datfile))
    hhdata = d['data']
    d.close()
    hhpludat = hht3.filter_counts_on_marker(hhdata, 2, [-1,0])
    hhpludat1 = hht3.filter_counts_on_marker(hhdata, 2, [-1])
    hhpludat2 = hht3.filter_counts_on_marker(hhdata, 2, [0])

    ### some sanity checks
    print ''
    print 'HH SSROs: %d' % len(hhpludat[
        np.logical_and(hhpludat[:,3]==1, hhpludat[:,2]==2)])
    print ''
    print 'counts + PLU: %d' % len(hhpludat[hhpludat[:,3]==0])
    print 'counts in window 1 + PLU: %d' % len(hhpludat1[hhpludat1[:,3]==0])
    print 'counts in window 2 + PLU: %d' % len(hhpludat2[hhpludat2[:,3]==0])
    
    return hhpludat

def correlations(hhpludat, SSROs_adwin, w1_start=(0,0), w1_stop=(700,700),
        w2_start=(0,0), w2_stop=(700,700), do_plot = True, savefile='correlations.pdf'):

    # do timefilter first
    photons = hhpludat[:,3] == 0
    mrkrs = np.logical_not(photons)

    # compute actual correlations             
    windows = np.zeros(len(hhpludat[hhpludat[:,3]==1]), dtype=int)
    
    for _i, nsync in np.ndenumerate(hhpludat[hhpludat[:,3]==1,0]):        
        i = _i[0]
        w1 = hhpludat[np.logical_and(photons, hhpludat[:,0] == nsync-1)]
        w2 = hhpludat[np.logical_and(photons, hhpludat[:,0] == nsync)]

        # print w1, w2        
        try:
            w1ch0 = w1[w1[:,2] == 0]
            w1ch0 = w1ch0[w1ch0[:,1] > w1_start[0]]
            w1ch0 = w1ch0[w1ch0[:,1] < w1_stop[0]]
            w1ch0len = len(w1ch0)
        except IndexError:
            w1ch0len = 0

        try:
            w1ch1 = w1[w1[:,2] == 1]
            w1ch1 = w1ch1[w1ch1[:,1] > w1_start[1]]
            w1ch1 = w1ch1[w1ch1[:,1] < w1_stop[1]]
            w1ch1len = len(w1ch1)
        except IndexError:
            w1ch1len = 0

        try:
            w2ch0 = w2[w2[:,2] == 0]
            w2ch0 = w2ch0[w2ch0[:,1] > w2_start[0]]
            w2ch0 = w2ch0[w2ch0[:,1] < w2_stop[0]]
            w2ch0len = len(w2ch0)
        except IndexError:
            w2ch0len = 0

        try:
            w2ch1 = w2[w2[:,2] == 1]
            w2ch1 = w2ch1[w2ch1[:,1] > w2_start[1]]
            w2ch1 = w2ch1[w2ch1[:,1] < w2_stop[1]]
            w2ch1len = len(w2ch1)
        except IndexError:
            w2ch1len = 0

        w1len = w1ch0len + w1ch1len
        w2len = w2ch0len + w2ch1len

        # print w1len, w2len
        
        if w1len == 1 and w2len == 0:
            windows[i] = 1
        elif w2len == 1 and w1len == 0:
            windows[i] = 2
        else:
            pass

    corr = np.zeros(4)
    total = len(windows[windows>0])
    print len(windows)
    print 'total = ',total
    print 'ssro',len(SSROs_adwin)
    for _i, w in np.ndenumerate(windows):
        i = _i[0]
        if w == 0 or SSROs_adwin[i]<0:
            continue
        corr[(w-1)*2**1 + int(SSROs_adwin[i]>0)*2**0] += 1

    # corr = corr/float(total)
    if do_plot:
        fig = plt.figure()
        plt.bar(np.arange(4)-0.4, corr)
        plt.xticks([])
        plt.ylabel('fraction of events')
        plt.text(0, 1, "Window 1\n0 cts", ha='center', va='bottom')
        plt.text(1, 1, "Window 1\n>0 cts", ha='center', va='bottom')
        plt.text(2, 1, "Window 2\n0 cts", ha='center', va='bottom')
        plt.text(3, 1, "Window 2\n>0 cts", ha='center', va='bottom')
        
        y = plt.ylim()[1]*0.9
        plt.text(0, y, corr[0], ha='center', va='top')
        plt.text(1, y, corr[1], ha='center', va='top')
        plt.text(2, y, corr[2], ha='center', va='top')
        plt.text(3, y, corr[3], ha='center', va='top')

        fig.savefig(savefile, format='pdf')

    return corr

def SSRO_correct(corr, F0, F1, F0err = 0.01, F1err = 0.01):
    w1 = corr[:2].copy()
    w2 = corr[2:].copy()
    w1 = w1[::-1]
    w2 = w2[::-1]
    
    w1total = np.sum(w1)
    w2total = np.sum(w2)
    w1 /= float(w1total)
    w2 /= float(w2total)

    norm = 1./(F0+F1-1.)
    w1staterr = np.sqrt(w1*(1.-w1)/float(w1total))
    w2staterr = np.sqrt(w2*(1.-w2)/float(w2total))
   
    w1p0 = w1[1]
    w1p1 = w1[0]
    w2p0 = w2[1]
    w2p1 = w2[0]

    w1ms0err = np.sqrt( 
            (norm**2 * (-w1p1 + F0*(w1p0+w1p1)) * F1err)**2 +\
            (norm**2 * (w1p0 - F1*(w1p0+w1p1)) * F0err)**2 +\
            (norm * (F1 - 1) * w1staterr[1])**2 +\
            (norm * F1 * w1staterr[0])**2
            )
    w1ms1err = w1ms0err
    w1err = np.array([w1ms0err, w1ms1err])

    w2ms0err = np.sqrt( 
            (norm**2 * (-w2p1 + F0*(w2p0+w2p1)) * F1err)**2 +\
            (norm**2 * (w2p0 - F1*(w2p0+w2p1)) * F0err)**2 +\
            (norm * (F1 - 1) * w2staterr[1])**2 +\
            (norm * F1 * w2staterr[0])**2
            )
    w2ms1err = w2ms0err
    w2err = np.array([w2ms0err, w2ms1err])

    corrmat = 1./(F0*F1 - (1.-F1)*(1.-F0)) * np.array([[F1, F1-1.],[F0-1., F0]])
    
    return np.dot(corrmat, w1.reshape(2,1)).reshape(-1), \
            np.dot(corrmat, w2.reshape(2,1)).reshape(-1), \
            w1err, w2err

def plot_counts_vs_time(windowdata, savepath, window=1, bins=700,
        save_prefix='', save_postfix=''):
    
    gated = hht3.filter_counts_on_marker(windowdata, 2, [window-2])
    
    print 'length of window data = ',len(windowdata)

    fig = plt.figure()
    
    t0, t1 = hht3.get_click_times(windowdata)
    h0y, h0x = np.histogram(t0, bins=bins, range=(0,bins-1))
    h1y, h1x = np.histogram(t1, bins=bins, range=(0,bins-1))
    
    g_t0, g_t1 = hht3.get_click_times(gated)
    g_h0y, g_h0x = np.histogram(g_t0, bins=bins, range=(0,bins-1))
    g_h1y, g_h1x = np.histogram(g_t1, bins=bins, range=(0,bins-1))
    
    ax = plt.subplot(121)
    ax.semilogy(g_h0x[:-1], g_h0y, color = 'g', label='ch0 gated')
    ax.semilogy(h0x[:-1], h0y, label='ch0', color = 'b')
    ax.set_xlabel('Bins')
    ax.legend()

    ax = plt.subplot(122)
    ax.semilogy(h1x[:-1], h1y, label='ch1')
    ax.semilogy(g_h1x[:-1], g_h1y, label='ch1 gated')
    ax.set_xlabel('Bins')
    ax.legend()

    plt.suptitle('Clicks after optical pulse, window %d' % window)
    
    fig.savefig(os.path.join(savepath,'%shistogram_window%s%s.pdf' % \
            (save_prefix,window,save_postfix)))

    np.savez(os.path.join(savepath, 'histogram_window%s.npz'%window), 
        window_number = window, 
        window_data = windowdata, 
        ch0 = h0y, 
        ch0_xaxis = h0x[:-1],
        ch0_gated = g_h0y,
        ch0_gated_xaxis = g_h0x[:-1],
        ch1 = h1y, 
        ch1_xaxis = h1x[:-1],
        ch1_gated = g_h1y,
        ch1_gated_xaxis = g_h1x[:-1])

    return True

def filter_PLU(data, mchan=2):
    filtered = hht3.filter_counts_on_marker(data, mchan, delta_syncs=[0])
    
    return filtered

# in case of large data sets, we need to treat data in chunks;
# can combine back after reducing to some reasonable size
def split(rawfolder, file_contains='LDE', print_status=False, *arg, **kw):
    
    windows = kw.pop('windows', [1,2])
    ch0maxtime = kw.pop('ch0maxtime', 700)
    ch1maxtime = kw.pop('ch1maxtime', 700)

    allfiles = [ f for f in os.listdir(rawfolder) if file_contains in f and \
            f[:1] != '.' and f[-4:] == '.npz' ]
    lastidx = len(allfiles)

    # we don't have too much data left over, so can use this approach
    finaldata = np.zeros((0,4), dtype=np.uint)
    windowdata = [ np.zeros((0,4), dtype=np.uint) for w in windows ]
        
    sync_of = 0
    for i,f in enumerate(allfiles):
        if print_status:
            print 'file %d (%d) ...' % (i+1, len(allfiles))
        
        fp = os.path.join(rawfolder, f)
        data = hht3.data_from_file(fp, ch0_lastbin=ch0maxtime, 
                ch1_lastbin=ch1maxtime, *arg, **kw)
        
        if len(data) == 0:
            continue

        data[:,0] += sync_of + 100
        sync_of = data[-1,0] + 1

        data = hht3.filter_counts_on_marker(data, 1, windows)
        if len(data) == 0:
            continue
        
        _windowdata = []
        for w in windows:
            _wd = hht3.filter_counts_on_marker(data, 1, [w])
            _wd = hht3.delete_markers(_wd,1)
            _windowdata.append(_wd)
            
        data = hht3.delete_markers(data, 1)

        finaldata = np.vstack((finaldata, data))
        for i,w in enumerate(windows):
            windowdata[i] = np.vstack((windowdata[i], _windowdata[i]))

        ret = [finaldata]
        ret.extend(windowdata)
        
    return tuple(ret)
load_hhdata = split

def barplot_correlations(datfolder, correlations, F0, F1, F0err, F1err):

    a = correlations[[0,1]]/float(sum(correlations[[0,1]]))
    b = correlations[[2,3]]/float(sum(correlations[[2,3]]))

    corrected = SSRO_correct(correlations, F0, F1, F0err, F1err)

# corrected: first element: 

    print 'corrected first window:',corrected[0]
    print 'corrected second window:',corrected[1]

    def num2str(num, precision): 
        return "%0.*f" % (precision, num)

    alpha = 0.5
    text_offset = 0.02

    figure2 = plt.figure(figsize = (12,8))
    plt.subplot(121)
    plt.title('Uncorrected')
    plt.bar(np.arange(2), a[::-1], color = 'g', align = 'center', 
            alpha = alpha, label = 'first window')
    plt.bar(np.arange(2,4),b[::-1], color = 'r', align = 'center', 
            alpha = alpha,label = 'second window')
    plt.legend(shadow = True, fancybox = True, loc = 3)
    plt.ylabel('Probability')
    plt.xticks([0,1,2,3],['$m_s = 0$','$m_s = -1$','$m_s = 0$','$m_s = -1$'])
    plt.ylim([0,1])
    for k in range(0,4):
        if k < 2:
            plt.text(k,a[::-1][k]+text_offset, 
                    num2str(a[::-1][k],2), 
                    ha = 'center', va = 'bottom')
        else:
            plt.text(k,b[::-1][k-2]+text_offset, 
                    num2str(b[::-1][k-2],2), 
                    ha = 'center', va = 'bottom')

    plt.subplot(122)
    plt.title('Corrected for read-out')
    plt.bar(np.arange(2), corrected[0], yerr=corrected[2], 
            color = 'g', ecolor='k',
            align = 'center', lw=1,
            label = 'first window')
    plt.bar(np.arange(2,4),corrected[1], yerr=corrected[3], 
            color = 'r', ecolor='k',
            align = 'center', lw=1,
            label = 'second window')
    plt.legend(shadow = True, fancybox = True, loc = 3)
    plt.ylabel('Probability')
    plt.xticks([0,1,2,3],['$m_s = 0$','$m_s = -1$','$m_s = 0$','$m_s = -1$'])
    plt.ylim([0,1])
    for k in range(0,4):
        if k < 2:
            plt.text(k,corrected[0][k]+text_offset, 
                    num2str(corrected[0][k],2), 
                    ha = 'center', va = 'bottom')
        else:
            plt.text(k,corrected[1][k-2]+text_offset, 
                    num2str(corrected[1][k-2],2), 
                    ha = 'center', va = 'bottom')


    if not os.path.exists(os.path.join(datfolder,'figures')):
        os.makedirs(os.path.join(datfolder,'figures'))

    figure2.savefig(os.path.join(datfolder,'figures','correlation_diagram.png'))
    
    raw_pops = dict()
    raw_pops['w1ms0'] = a[::-1][0]
    raw_pops['w1ms1'] = a[::-1][1]
    raw_pops['w2ms0'] = b[::-1][0]
    raw_pops['w2ms1'] = b[::-1][1]

    corr_pops = dict()
    corr_pops['w1ms0'] = corrected[0][0]
    corr_pops['w1ms1'] = corrected[0][1]
    corr_pops['w2ms0'] = corrected[1][0]
    corr_pops['w2ms1'] = corrected[1][1]
    corr_pops['w1ms0_err'] = corrected[2][0]
    corr_pops['w1ms1_err'] = corrected[2][1]
    corr_pops['w2ms0_err'] = corrected[3][0]
    corr_pops['w2ms1_err'] = corrected[3][1]

    F = dict()
    F['F0'] = F0
    F['F1'] = F1
    F['F0_err'] = F0err
    F['F1_err'] = F1err

    np.savez(os.path.join(datfolder, 'correlations.npz'),
            raw = raw_pops, corrected = corrected, 
            corr_pops = corr_pops, F = F)

    return 
