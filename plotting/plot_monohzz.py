from __future__ import print_function, division
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)

blind = True

from Plotter import Plotter
from samples import get_sample_list

def plot(year):
    base = 'hzzProcessor'
    
    backgrounds = ['ggZZ','qqZZ','HZZ']
    data = 'DATA'
    signals = []
    
    sampleMap = { s: get_sample_list(base,year,s) for s in backgrounds+signals+[data] }
    print(sampleMap)
    
    # load the tfiles
    channels = ['2m2e','2e2m','4e','4m']
    tfiles = {}
    for s in sampleMap:
        if None in sampleMap[s]:
            sampleMap[s] = [x for x in sampleMap[s] if x is not None]
            print(f'Warning: missing samples for {s}')
        for sample in sampleMap[s]:
            tfiles[sample] = ROOT.TFile.Open(f'hists/hzzProcessor/{year}/{sample}.root')
    
    # styles for plots
    styleMap = {
        'TT'  : {'label': 't#bar{t}',                             'linecolor': ROOT.kGreen+3,   'fillcolor': ROOT.kGreen+3,},
        'TTV' : {'label': 't#bar{t}V',                            'linecolor': ROOT.kGreen+4,   'fillcolor': ROOT.kGreen+4,},
        'Z'   : {'label': 'Z+X',                                  'linecolor': ROOT.kGreen+2,   'fillcolor': ROOT.kGreen+2,},
        'ggZZ': {'label': 'gg#rightarrowZZ, Z#gamma*',            'linecolor': ROOT.kBlue,      'fillcolor': ROOT.kBlue,},
        'qqZZ': {'label': 'q#bar{q}#rightarrowZZ, Z#gamma*',      'linecolor': ROOT.kAzure+6,   'fillcolor': ROOT.kAzure+6,},
        'HWW' : {'label': 'H(125)#rightarrowWW#rightarrow2l2#nu', 'linecolor': ROOT.kRed+2,     'fillcolor': ROOT.kRed+2,},
        'HZZ' : {'label': 'H(125)#rightarrowZZ#rightarrow4l',     'linecolor': ROOT.kRed+1,     'fillcolor': ROOT.kRed+1,},
        'SIG' : {'label': '2HDM+a',                               'linecolor': ROOT.kMagenta+1,},
        'DATA': {'label': 'Observed',}
    }
    
    # setup plotter
    plotter = Plotter('MonoHZZ',year)
    plotter.setStyleMap(styleMap)
    for bg in backgrounds: plotter.addSampleToStack(bg)
    for sig in signals: plotter.addSampleToPlot(sig)
    plotter.setDataSample(data)
    
    # plot
    plots = {
        'm4l'       : {'hPath': '{sample}_{chan}_hzz_mass',         'xlabel': 'm_{4l} (GeV)',       'binning': range(70,500,4),         'ylabel': 'Events / 4 GeV', 'logx': False, 'logy': False,},
        'm4l_zoom'  : {'hPath': '{sample}_{chan}_hzz_mass',         'xlabel': 'm_{4l} (GeV)',       'binning': range(70,170,2),         'ylabel': 'Events / 2 GeV', 'logx': False, 'logy': False,},
        'm4l_full'  : {'hPath': '{sample}_{chan}_massWindow_mass',  'xlabel': 'm_{4l} (GeV)',       'binning': range(113,135,3),        'ylabel': 'Events / 3 GeV', 'logx': False, 'logy': False, 'blind':blind,},
        'mz1'       : {'hPath': '{sample}_{chan}_hzz_z1mass',       'xlabel': 'm_{ll} (GeV)',       'binning': range(40,120,1),         'ylabel': 'Events / 1 GeV', 'logx': False, 'logy': False,},
        'mz2'       : {'hPath': '{sample}_{chan}_hzz_z2mass',       'xlabel': 'm_{ll} (GeV)',       'binning': range(12,120,1),         'ylabel': 'Events / 1 GeV', 'logx': False, 'logy': False,},
        'met'       : {'hPath': '{sample}_{chan}_hzz_met',          'xlabel': 'E_{T}^{miss} (GeV)', 'binning': range(0,420,20),         'ylabel': 'Events / 20 GeV','logx': False, 'logy': True, 'ymin':0.1,},
        'met_limit' : {'hPath': '{sample}_{chan}_massWindow_met',   'xlabel': 'E_{T}^{miss} (GeV)', 'binning': [0,25,50,200,500,1000],  'ylabel': 'Events',         'logx': False, 'logy': True, 'ymin':0.1, 'blind':blind,},
        'pileup'    : {'hPath': '{sample}_{chan}_hzz_pileup',       'xlabel': 'Number of reconstructed vertices', 'binning': range(0,120,1), 'ylabel': 'Events',    'logx': False, 'logy': False,},
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
