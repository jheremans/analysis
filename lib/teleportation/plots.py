import os

from matplotlib import pyplot as plt
import numpy as np

from analysis.lib.math import readout_correction as roc

import files, events, settings, stats

##############################################################################
### Photon histograms in time
##############################################################################

def _plot_photon_hist(ax, h, b, log=True, **kw):
    label = kw.pop('label', '')

    _h = h.astype(float)
    _h[_h<=1e-1] = 1e-1
    _h = np.append(_h, _h[-1])
           
    ax.plot(b, _h, drawstyle='steps-post', label=label)
    if log:
        ax.set_yscale('log')
    ax.set_xlabel('time (ns)')
    ax.set_ylabel('events')
    ax.set_ylim(bottom=0.1)
    ax.set_xlim(min(b), max(b))

def plot_photon_hist(fp, **kw):    
    ret = kw.pop('ret', 'subplots')

    (h0, b0), (h1, b1) = events.get_photon_hist(fp, **kw)
   
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0)
    _plot_photon_hist(ax1, h1, b1)

    ax0.set_title('photons channel 0')
    ax1.set_title('photons channel 1')
    
    fig.suptitle(files.get_msmt_header(fp) + ' -- Photon histogram')
    
    if ret == 'subplots':
        return fig, (ax0, ax1)
    
def plot_photon_hist_filter_comparison(fp, fltr, **kw):
    ret = kw.pop('ret', 'subplots')
    
    (h0, b0), (h1, b1) = events.get_photon_hist(fp, **kw)
    (h0f, b0f), (h1f, b1f) = events.get_photon_hist(fp, fltr=fltr, **kw)
    
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0, label='unfiltered')
    _plot_photon_hist(ax0, h0f, b0f, label='filtered')

    _plot_photon_hist(ax1, h1, b1, label='unfiltered')
    _plot_photon_hist(ax1, h1f, b1f, label='filtered')
    
    ax0.set_title('photons channel 0')
    ax1.set_title('photons channel 1')

    fig.suptitle(files.get_msmt_header(fp) + ' -- Photon histogram, comparison w/ filter')
        
    if ret=='subplots':
        return fig, (ax0, ax1)


def plot_tail(fp, **kw):
    ret = kw.pop('ret', 'subplots')
    
    (h0, b0), (h1, b1) = events.get_photon_hist(fp, **kw)
    
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0)
    _plot_photon_hist(ax1, h1, b1)
    
    ax0.set_xlim(settings.CH0_START-5, settings.CH0_START + 100)
    ax1.set_xlim(settings.CH1_START-5, settings.CH1_START + 100)
    ax0.set_title('tail channel 0')
    ax1.set_title('tail channel 1')
    
    fitres0, fitres1 = events.fit_tail(fp)
    
    ax0.plot(fitres0['fit_x'], fitres0['fit_y'], lw=2)
    ax1.plot(fitres1['fit_x'], fitres1['fit_y'], lw=2)
    
    ax0.text(settings.CH0_START, 1, fitres0['info'], va='bottom', ha='left', color='k', size='x-small',
             bbox = {'facecolor' : 'white', 'alpha' : 0.5})
    ax1.text(settings.CH1_START, 1, fitres1['info'], va='bottom', ha='left', color='k', size='x-small',
             bbox = {'facecolor' : 'white', 'alpha' : 0.5})
    
    starts = stats.get_sequence_starts(fp)
    ax0.text(settings.CH0_START+5, max(h0)/5., 
             'tailcounts / shot: {:.2e} ({} x {} starts)'.format(fitres0['tcpsh'], starts, settings.SEQREPS))
    ax1.text(settings.CH1_START+5, max(h1)/5., 
             'tailcounts / shot: {:.2e} ({} x {} starts)'.format(fitres1['tcpsh'], starts, settings.SEQREPS))

    fig.suptitle(files.get_msmt_header(fp) + ' -- Tail')

    if ret=='subplots':
        return fig, (ax0, ax1)

def plot_PLU_filter(folder):
    fp = files.get_msmt_fp(folder)

    # get the PLU marked photons first
    is_ph_ch0, is_ph_ch1 = events.get_photons(fp)
    is_ph = is_ph_ch0 | is_ph_ch1
    is_ph_with_PLU_mrkr = is_ph & events.filter_marker(fp, 2)

    # first window
    fig, (ax0, ax1) = plot_photon_hist_filter_comparison(
        fp, save=False, fltr=is_ph_with_PLU_mrkr, log=True,
        binedges=settings.PHOTONHIST_BINEDGES)

    ax0.axvline(settings.CH0_START, color='k', linestyle='--')
    ax0.text(settings.CH0_START+1, 1, '{} ns'.format(settings.CH0_START), color='k')
    ax0.legend()
    ax0.set_title('PLU (hardware) filter first pi-pulse ch0')

    ax1.axvline(settings.CH1_START, color='k', linestyle='--')
    ax1.text(settings.CH1_START+1, 1, '{} ns'.format(settings.CH1_START), color='k')
    ax1.legend()
    ax1.set_title('PLU (hardware) filter first pi-pulse ch1')

    fig.savefig(os.path.join(folder, 'PLU_photons_window1.png'))

    # second window
    fig, (ax0, ax1) = plot_photon_hist_filter_comparison(
        fp, save=False, fltr=is_ph_with_PLU_mrkr, log=True,
        binedges=settings.PHOTONHIST_BINEDGES + settings.PIPULSESEP)

    ax0.axvline(settings.CH0_START+settings.PIPULSESEP, color='k', linestyle='--')
    ax0.text(settings.CH0_START+settings.PIPULSESEP+1, 
        1, '{} ns'.format(settings.CH0_START+settings.PIPULSESEP), color='k')
    ax0.legend()
    ax0.set_title('PLU (hardware) filter second pi-pulse ch0')

    ax1.axvline(settings.CH1_START+settings.PIPULSESEP, color='k', linestyle='--')
    ax1.text(settings.CH1_START+settings.PIPULSESEP+1, 
        1, '{} ns'.format(settings.CH1_START+settings.PIPULSESEP), color='k')
    ax1.legend()
    ax1.set_title('PLU (hardware) filter second pi-pulse ch1')
     
    fig.savefig(os.path.join(folder, 'PLU_photons_window2.png'))


##############################################################################
### CR stat plotting
##############################################################################

def _plot_CR_hist(ax, heights, **kw):
    ax.bar(np.arange(len(heights)), heights, align='center', 
           width=1, color=settings.COLORS[0], lw=0)
    ax.set_xlabel('CR counts')
    ax.set_ylabel('Events')
    ax.set_title(kw.get('title', ''))
    ax.set_xlim(-0.5, len(heights)-0.5 if len(heights) > 0 else 0.5)
   
### plotting of histograms
def plot_CR_hist_total(fp):
    h1, b1 = stats.get_CR_hist_total(fp, 'lt1')
    h2, b2 = stats.get_CR_hist_total(fp, 'lt2')
    
    fig, (ax1, ax2) = files.subplots(1,2, figsize=(12,5))
    _plot_CR_hist(ax1, h1, title='LT1')
    _plot_CR_hist(ax2, h2, title='LT2')
    
    fig.suptitle(files.get_msmt_header(fp) + ' -- CR histogram total')
    folder, _fn = os.path.split(fp)
    fig.savefig(os.path.join(folder, 'CR_histogram_total.png'))
    
def plot_CR_hist_sequence_timeout(fp):
    h1, b1 = stats.get_CR_hist_sequence_timeout(fp, 'lt1')
    h2, b2 = stats.get_CR_hist_sequence_timeout(fp, 'lt2')
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))
    _plot_CR_hist(ax1, h1, title='LT1')
    _plot_CR_hist(ax2, h2, title='LT2')
    
    fig.suptitle(files.get_msmt_header(fp) + ' -- CR histogram sequence timeout')
    folder, _fn = os.path.split(fp)
    fig.savefig(os.path.join(folder, 'CR_histogram_sequence_timeout.png'))

def plot_CR_hist_before_after(fp):
    h1b, b1b = stats.get_CR_hist_before(fp, 'lt1')
    h1a, b1a = stats.get_CR_hist_after(fp, 'lt1')

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))
    ax1.bar(np.arange(len(h1b)), h1b, align='center',
        width=1, color=settings.COLORS[0], lw=0, alpha=0.7,
        label='before')
    ax1.bar(np.arange(len(h1a)), h1a, align='center',
        width=1, color=settings.COLORS[1], lw=0, alpha=0.7,
        label='after')
    ax1.set_xlabel('CR counts')
    ax1.set_ylabel('Events')
    ax1.set_title('LT1')
    ax1.set_xlim(left=-0.5)
    ax1.legend(loc='best')

    h2b, b2b = stats.get_CR_hist_before(fp, 'lt2')
    h2a, b2a = stats.get_CR_hist_after(fp, 'lt2')

    ax2.bar(np.arange(len(h2b)), h2b, align='center',
        width=1, color=settings.COLORS[0], lw=0, alpha=0.7,
        label='before')
    ax2.bar(np.arange(len(h2a)), h2a, align='center',
        width=1, color=settings.COLORS[1], lw=0, alpha=0.7,
        label='after')
    ax2.set_xlabel('CR counts')
    ax2.set_title('LT2')
    ax2.set_xlim(left=-0.5)
    ax2.legend(loc='best')

    fig.suptitle(files.get_msmt_header(fp) + ' -- CR histogram before/after')
    folder, _fn = os.path.split(fp)
    fig.savefig(os.path.join(folder, 'CR_histogram_before-after.png'))


##############################################################################
### Readout correlations
##############################################################################

def plot_raw_readout_correlations(BSM_outcomes_psiminus, BSM_outcomes_psiplus, 
    Bob_outcomes_psiminus, Bob_outcomes_psiplus, **kw):
    
    BSM_outcome_names = ['00', '01', '10', '11']

    fig, ((ax_bsm_psiminus,ax_bob_psiminus), (ax_bsm_psiplus, ax_bob_psiplus)) = \
        plt.subplots(2,2, figsize=(10,8), sharex=True)

    for outcome in range(4):
        ax_bsm_psiminus.bar(outcome, BSM_outcomes_psiminus[outcome], width=0.7, 
                            color=settings.COLORS[0], ecolor='k')    
        # ax_bsm_psiminus.text(outcome+0.35, 1.01,
        #          '{}'.format(int(Bob_outcomes_psiminus[outcome,:].sum())),
        #          ha='center', va='bottom',
        #          color=settings.COLORS[0])
        # ax_bsm_psiminus.text(outcome+0.35, BSM_probs_psiminus[outcome] + 0.05,
        #          "{:.1f}%".format(BSM_probs_psiminus[outcome]*100),
        #          va='bottom', ha='center')
        
        ax_bob_psiminus.bar(outcome, Bob_outcomes_psiminus[outcome,0],
                            width=0.3, color=settings.COLORS[1],
                            label = 'ms = 0' if outcome == 0 else '')
        ax_bob_psiminus.bar(outcome+0.35, Bob_outcomes_psiminus[outcome,1],
                            width=0.3, color=settings.COLORS[2],
                            label = 'ms = 1' if outcome == 0 else '')
        # ax_bob_psiminus.text(outcome+0.35, 1.01,
        #          '{}'.format(int(Bob_outcomes_psiminus[outcome,0])),
        #          ha='center', va='bottom',
        #          color=settings.COLORS[1])
        
        ax_bsm_psiplus.bar(outcome, BSM_outcomes_psiplus[outcome], width=0.7, 
                            color=settings.COLORS[0], ecolor='k')    
        # ax_bsm_psiplus.text(outcome+0.35, 1.01,
        #          '{}'.format(int(Bob_outcomes_psiplus[outcome,:].sum())),
        #          ha='center', va='bottom',
        #          color=settings.COLORS[0])
        # ax_bsm_psiplus.text(outcome+0.35, BSM_probs_psiplus[outcome] + 0.05,
        #          "{:.1f}%".format(BSM_probs_psiplus[outcome]*100),
        #          va='bottom', ha='center')
        
        ax_bob_psiplus.bar(outcome, Bob_outcomes_psiplus[outcome,0],
                            width=0.3, color=settings.COLORS[1])
        ax_bob_psiplus.bar(outcome+0.35, Bob_outcomes_psiplus[outcome,1],
                            width=0.3, color=settings.COLORS[2])
        
        # ax_bob_psiplus.text(outcome+0.35, 1.01,
        #          '{}'.format(int(Bob_outcomes_psiplus[outcome,0])),
        #          ha='center', va='bottom',
        #          color=settings.COLORS[1])

    ax_bsm_psiminus.set_xlim(-0.1, 3.8)
    ax_bsm_psiminus.set_xticks(np.arange(4)+0.35)
    ax_bsm_psiminus.set_xticklabels(BSM_outcome_names)
    ax_bsm_psiminus.set_title('Psi- BSM outcomes')
    ax_bsm_psiminus.set_ylabel('Events')
    # ax_bsm_psiminus.set_ylim(-0.05, 1.2)
    # ax_bsm_psiminus.axhline(0, c='k', ls=':')
    # ax_bsm_psiminus.axhline(1, c='k', ls=':')

    # ax_bob_psiminus.set_ylim(-0.05, 1.2)
    # ax_bob_psiminus.axhline(0, c='k', ls=':')
    # ax_bob_psiminus.axhline(1, c='k', ls=':')
    ax_bob_psiminus.set_title('Psi- Bob readout')
    ax_bob_psiminus.set_ylabel('Conditional results for Bob')
    ax_bob_psiminus.legend(loc='best')

    ax_bsm_psiplus.set_title('Psi+ BSM outcomes')
    ax_bsm_psiplus.set_ylabel('Events')
    # ax_bsm_psiplus.set_ylim(-0.05, 1.2)
    # ax_bsm_psiplus.axhline(0, c='k', ls=':')
    # ax_bsm_psiplus.axhline(1, c='k', ls=':')

    # ax_bob_psiplus.set_ylim(-0.05, 1.2)
    # ax_bob_psiplus.axhline(0, c='k', ls=':')
    # ax_bob_psiplus.axhline(1, c='k', ls=':')
    ax_bob_psiplus.set_title('Psi+ Bob readout')
    ax_bob_psiplus.set_ylabel('Conditional results for Bob')
       
    ax_bsm_psiplus.set_xlabel('BSM results (N,e)')
    ax_bob_psiplus.set_xlabel('BSM results (N,e)')

    fig.suptitle('BSM / Bob SSRO Correlations')
    # fig.savefig(os.path.join(folder, 'BSM-BobSSRO-correlations.png'))


def plot_single_roc(zero_events, one_events, *ssro_fids, **kw):
    fig, ax = plt.subplots(1,1, figsize=(4.5,4))

    frac0 = float(zero_events)/(zero_events+one_events)
    u_frac0 = np.sqrt(frac0*(1.-frac0)/(zero_events+one_events))

    p0, u_p0 = roc.single_qubit_outcome_with_ROC(zero_events, zero_events+one_events, *ssro_fids)

    ax.bar(range(2), [frac0, p0], 
        color=[settings.COLORS[0], settings.COLORS[1]],
        align='center', yerr=[u_frac0, u_p0], 
        ecolor='k', width=0.8)

    ax.text(0, 1.05, '{:.0f}+/-{:.0f} %'.format(frac0*100., u_frac0*100.),
        ha='center', va='bottom', color=settings.COLORS[0])
    ax.text(1, 1.05, '{:.0f}+/-{:.0f} %'.format(p0*100., u_p0*100.),
        ha='center', va='bottom', color=settings.COLORS[1])  

    ax.set_xlim(-0.5,1.5)
    ax.set_ylim(-0.05, 1.15)
    ax.set_ylabel('P ($m_s =\, 0$)')

    ax.set_xticks([0,1])
    ax.set_xticklabels(['Raw', 'Corrected'])
    ax.axvline(0.5, c='k', ls='--')
    ax.axhline(0, c='k', ls=':')
    ax.axhline(1, c='k', ls=':')

    return fig, ax

##############################################################################
### Teleportation event stats
##############################################################################

def plot_success_attempt_distribution(attempts):
    fig, ax = plt.subplots(1,1, figsize=(15,1.5))

    h,b = np.histogram(attempts, 
        bins=np.arange(-0.5, settings.SEQREPS))

    im = ax.imshow(h.reshape((1, settings.SEQREPS)),
        aspect='auto', cmap='Blues')
    ax.set_xlabel('attempt')
    ax.set_yticks([])

    cbar = fig.colorbar(im, orientation='vertical',
        aspect=3, ticks=[0, h.max()])
    cbar.set_label('occurrences')

    # ax.hist(attempts,
    #     bins=np.arange(-0.5, settings.SEQREPS),
    #     ec='None', color=settings.COLORS[0])

    ax.set_title('Successful LDE attempts')

    return fig, ax






