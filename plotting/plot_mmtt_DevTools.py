from __future__ import print_function, division
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)

blind = True

from Plotter import Plotter
from samples import get_sample_list

def plot(year):
    base = 'mmttProcessor_DevTools'

    hamasses = [(125,5,), (125,7,), (125,9,), (125,11), (125,13), (125,15), (125,17), (125,19)]

    # define the samples
    backgrounds = []
    data = 'DATA'
    signals = [f'haaH{h}A{a}' for h, a in hamasses]
    sampleMap = { s: get_sample_list(base,year,s) for s in backgrounds+signals+[data] }
    print(sampleMap)

    # load the tfiles
    channels = ['mm']
    cats = ['dimuon','ditau','ditauMD']
    tfiles = {}
    for s in sampleMap:
        for sample in sampleMap[s]:
            # TODO: we can keep this in the sample map builder
            tfiles[sample] = ROOT.TFile.Open(f'hists/{base}/{year}/{sample}.root')
    
    # styles for plots
    styleMap = {
        'DATA': {'label': 'Observed',}
    }

    sigcolors = [ROOT.kRed+2, ROOT.kRed-1, ROOT.kRed-3,
                 ROOT.kGreen+2, ROOT.kGreen-1, ROOT.kGreen-3,
                 ROOT.kBlue+2, ROOT.kBlue-1, ROOT.kBlue-3]
    for i, (h,a) in enumerate(hamasses):
        styleMap[f'haaH{h}A{a}'] = {
            'label': f'ggH #rightarrow aa (m_{{H}} = {h} GeV, m_{{a}} = {a} GeV',
            'linecolor': sigcolors[i],  
        }
    
    # setup plotter
    plotter = Plotter('MMTT',year)
    plotter.setStyleMap(styleMap)
    for bg in backgrounds: plotter.addSampleToStack(bg)
    for sig in signals: plotter.addSampleToPlot(sig)
    plotter.setDataSample(data)
    
    # plot
    plots = {
        'mmMass'    : {'hPath': '{sample}_{chan}_{cat}_mmMass',        'xlabel': 'm_{#mu#mu} (GeV)',                 'binning': [x*0.1 for x in range(25,251,1)], 'ylabel': 'Events / 100 MeV', 'logx': False, 'logy': False, 'plotratio': False, 'blind': blind,},
        'ttMass'    : {'hPath': '{sample}_{chan}_{cat}_ttMass',        'xlabel': 'm_{#tau#tau} (GeV)',               'binning': [x*0.1 for x in range(0,251,5)], 'ylabel': 'Events / 0.5 GeV', 'logx': False, 'logy': False, 'plotratio': False, 'blind': blind,},
        'mmttMass'  : {'hPath': '{sample}_{chan}_{cat}_mmttMass',      'xlabel': 'm_{#mu#mu#tau#tau} (GeV)',         'binning': range(0,250,1),    'ylabel': 'Events / 1 GeV', 'logx': False, 'logy': False, 'plotratio': False, 'blind': blind,},
        'pileup'    : {'hPath': '{sample}_{chan}_{cat}_pileup',        'xlabel': 'Number of reconstructed vertices', 'binning': range(0,120,1),    'ylabel': 'Events',           'logx': False, 'logy': False, 'plotratio': False, 'blind': blind, },
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
    
    chan = 'mm'
    for plot in plots:
        for cat in cats:
            if plot in ['ttMass','mmttMass'] and cat=='dimuon': continue
            # load the histograms
            hists = {s:[] for s in sampleMap}
            for s in sampleMap:
                for sample in sampleMap[s]:
                    name = f'h_{plot}_{chan}_{cat}_{s}_{sample}'
                    hist  = tfiles[sample].Get(plots[plot]['hPath'].format(sample=sample,chan=chan,cat=cat))
                    if hist:
                        hists[s] += [hist.Clone(name)]
            # sum the histograms
            for s in hists:
                hname = f'{plot}_{chan}_{cat}_{s}'
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
            plotter.plot(hists, f'{year}/{chan}/{cat}/{plot}', **plots[plot])
    

years = ['2017']
for year in years:
    plot(year)
