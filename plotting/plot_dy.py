from __future__ import print_function, division
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)

blind = True

from Plotter import Plotter
from samples import get_sample_list

def plot(year):
    base = 'dyProcessor'

    # define the samples
    backgrounds = ['TT','Z']
    data = 'DATA'
    signals = []
    sampleMap = { s: get_sample_list(base,year,s) for s in backgrounds+signals+[data] }

    print(sampleMap)
    
    # load the tfiles
    channels = ['ee','mm']
    tfiles = {}
    for s in sampleMap:
        for sample in sampleMap[s]:
            # TODO: we can keep this in the sample map builder
            tfiles[sample] = ROOT.TFile.Open(f'hists/{base}/{year}/{sample}.root')
    
    # styles for plots
    styleMap = {
        'TT'  : {'label': 't#bar{t}',    'linecolor': ROOT.kGreen+3,   'fillcolor': ROOT.kGreen+3,},
        'Z'   : {'label': 'Drell-Yan',   'linecolor': ROOT.kOrange-2,  'fillcolor': ROOT.kOrange-2,},
        'DATA': {'label': 'Observed',}
    }
    
    # setup plotter
    plotter = Plotter('DY',year)
    plotter.setStyleMap(styleMap)
    for bg in backgrounds: plotter.addSampleToStack(bg)
    for sig in signals: plotter.addSampleToPlot(sig)
    plotter.setDataSample(data)
    
    # plot
    plots = {
        'mz1'       : {'hPath': '{sample}_{chan}_massWindow_zmass',        'xlabel': 'm_{ll} (GeV)',       'binning': range(60,120,1),         'ylabel': 'Events / 1 GeV', 'logx': False, 'logy': False,},
        'met'       : {'hPath': '{sample}_{chan}_massWindow_met',          'xlabel': 'E_{T}^{miss} (GeV)', 'binning': range(0,420,20),         'ylabel': 'Events / 20 GeV','logx': False, 'logy': True, 'ymin':0.1,},
    }
    
    def sumHists(name,*hists):
        hlist = ROOT.TList()
        for h in hists:
            if h: hlist.Add(h.Clone())
        if hlist.IsEmpty():
            print('No histograms for',name)
            return None
        hist = hists[0].Clone(name)
        hist.Reset()
        hist.Merge(hlist)
        return hist
    
    for plot in plots:
        # individual channels
        for chan in channels:
            # load the histograms
            hists = {s:[] for s in sampleMap}
            for s in sampleMap:
                for sample in sampleMap[s]:
                    name = f'h_{plot}_{chan}_{s}_{sample}'
                    hist  = tfiles[sample].Get(plots[plot]['hPath'].format(sample=sample,chan=chan))
                    if hist:
                        hists[s] += [hist.Clone(name)]
            # sum the histograms
            for s in hists:
                hname = f'{plot}_{chan}_{s}'
                hist = sumHists(hname,*hists[s])
                # bin the histogram
                if hist:
                    binning = array('d',plots[plot]['binning'])
                    hist = hist.Rebin(len(binning)-1,hname+'_rebin',binning)
                hists[s] = hist
                if s=='SIG': hists[s].Scale(0.001)
            # scale the hists:
            for s in hists:
                if s=='DATA': continue
                hists[s].Scale(float(plotter.intLumi)/1000)
            # send to plotter
            plotter.plot(hists, f'{year}/{chan}/{plot}', **plots[plot])
    
        # combined
        hists = {s:[] for s in sampleMap}
        for chan in channels:
            # load the histograms
            for s in sampleMap:
                for sample in sampleMap[s]:
                    name = f'h_{plot}_{chan}_{s}_{sample}'
                    hist  = tfiles[sample].Get(plots[plot]['hPath'].format(sample=sample,chan=chan))
                    if hist:
                        hists[s] += [hist.Clone(name)]
        # sum the histograms
        for s in hists:
            hname = f'{plot}_{s}'
            hist = sumHists(hname,*hists[s])
            # bin the histogram
            if hist:
                binning = array('d',plots[plot]['binning'])
                hist = hist.Rebin(len(binning)-1,hname+'_rebin',binning)
            hists[s] = hist
            if s=='SIG': hists[s].Scale(0.001)
        # scale the hists:
        for s in hists:
            if s=='DATA': continue
            hists[s].Scale(float(plotter.intLumi)/1000)
        # send to plotter
        plotter.plot(hists, f'{year}/{plot}', **plots[plot])

years = ['2016','2017','2018']
for year in years:
    plot(year)
