from __future__ import print_function, division
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)

blind = True

from Plotter import Plotter

year = '2018'

# define the samples
sampleMap = {
    'Z': [
        'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    ],
    'TT': [
        'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    ],
    'DATA': [
        'SingleMuon',
        'DoubleMuon',
        'EGamma',
    ],
}
backgrounds = ['TT','Z']
data = 'DATA'
signals = []

# load the tfiles
channels = ['ee','mm']
tfiles = {}
for s in sampleMap:
    for sample in sampleMap[s]:
        tfiles[sample] = ROOT.TFile.Open(f'hists/dyProcessor/{year}/{sample}.root')

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
        plotter.plot(hists, f'{chan}/{plot}', **plots[plot])

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
    plotter.plot(hists, f'{plot}', **plots[plot])
