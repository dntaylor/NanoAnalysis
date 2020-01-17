from __future__ import print_function, division
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)

blind = True

from Plotter import Plotter
from samples import get_sample_list

def plot(year):
    base = 'mmProcessor'

    # define the samples
    backgrounds = []
    data = 'DATA'
    signals = []
    sampleMap = { s: get_sample_list(base,year,s) for s in backgrounds+signals+[data] }
    print(sampleMap)

    # load the tfiles
    channels = ['mm']
    tfiles = {}
    for s in sampleMap:
        for sample in sampleMap[s]:
            # TODO: we can keep this in the sample map builder
            tfiles[sample] = ROOT.TFile.Open(f'hists/{base}/{year}/{sample}.root')
    
    # styles for plots
    styleMap = {
        'DATA': {'label': 'Observed',}
    }
    
    # setup plotter
    plotter = Plotter('MuMu',year)
    plotter.setStyleMap(styleMap)
    for bg in backgrounds: plotter.addSampleToStack(bg)
    for sig in signals: plotter.addSampleToPlot(sig)
    plotter.setDataSample(data)
    
    # plot
    plots = {
        'mll'       : {'hPath': '{sample}_{chan}_iso_mass',         'xlabel': 'm_{#mu#mu} (GeV)',       'binning': [x*0.1 for x in range(25,625,1)],          'ylabel': 'Events / 100 MeV', 'logx': False, 'logy': True, 'plotratio': False,},
        'jpsi'      : {'hPath': '{sample}_{chan}_iso_mass',         'xlabel': 'm_{#mu#mu} (GeV)',       'binning': [x*0.01 for x in range(250,450,1)],        'ylabel': 'Events / 10 MeV',  'logx': False, 'logy': True, 'plotratio': False,},
        'upsilon'   : {'hPath': '{sample}_{chan}_iso_mass',         'xlabel': 'm_{#mu#mu} (GeV)',       'binning': [x*0.01 for x in range(600,1400,1)],       'ylabel': 'Events / 10 MeV',  'logx': False, 'logy': False, 'plotratio': False,},
        'pileup'    : {'hPath': '{sample}_{chan}_iso_pileup',       'xlabel': 'Number of reconstructed vertices', 'binning': range(0,120,1), 'ylabel': 'Events',           'logx': False, 'logy': False, 'plotratio': False,},
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
