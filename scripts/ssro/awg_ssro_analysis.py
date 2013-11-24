from analysis.lib.m2 import m2
from analysis.lib.tools import toolbox
from matplotlib import pyplot as plt
import numpy as np

class AWGSSROAnalysis(m2.M2Analysis):

    def plot_f0_vs_time(self):

        dh=self.f['HH_sync_time-1'].value
        dd=self.f['HH_sync_number-1'].value

        k=100
        dds=np.zeros(100)
        ats=np.linspace(1e6,11e6,k)
        for i,at in enumerate(ats):
            dds[i]=len(unique(dd[dh<at]))/5000.
        print 'r'


        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        ax.plot(ats,dds)
        fig.savefig(os.path.join(self.folder, 'f0_vs_time.pdf'))



def awgssro(folder=tb.latest_data('AWGSSRO')):
    print folder
    a = AWGSSROAnalysis(folder)
    a.plot_f0_vs_time()
    a.finish()


if __name__ == '__main__':
    awgssro()