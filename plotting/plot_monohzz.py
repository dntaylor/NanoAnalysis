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
    'ZZ': [
        'ZZTo4L_13TeV_powheg_pythia8_TuneCP5',
    ],
    'ggZZ': [
        'GluGluToContinToZZTo2e2mu_13TeV_TuneCP5_MCFM701_pythia8',
        'GluGluToContinToZZTo2e2tau_13TeV_TuneCP5_MCFM701_pythia8',
        'GluGluToContinToZZTo2mu2tau_13TeV_TuneCP5_MCFM701_pythia8',
        'GluGluToContinToZZTo4e_13TeV_TuneCP5_MCFM701_pythia8',
        'GluGluToContinToZZTo4mu_13TeV_TuneCP5_MCFM701_pythia8',
        'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8',
    ],
    'HZZ': [
        'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
        'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
        'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',
        'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',
        'GluGluToZH_HToZZTo4L_M125_13TeV_JHUGenV723_pythia8',
        #'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8',
        'bbH_HToZZTo4L_M125_13TeV_JHUGenV7011_pythia8',
    ],
    'DATA': [
        'SingleMuon',
        'DoubleMuon',
        'MuonEG',
        'EGamma',
    ],
}
backgrounds = ['Z','TT','ggZZ','ZZ','HZZ']
data = 'DATA'
signals = []

# load the tfiles
channels = ['2m2e','2e2m','4e','4m']
tfiles = {}
for s in sampleMap:
    for sample in sampleMap[s]:
        tfiles[sample] = ROOT.TFile.Open(f'hists/hzzProcessor/{year}/{sample}.root')

# styles for plots
styleMap = {
    'TT'  : {'label': 't#bar{t}',                             'linecolor': ROOT.kGreen+3,   'fillcolor': ROOT.kGreen+3,},
    'TTV' : {'label': 't#bar{t}V',                            'linecolor': ROOT.kGreen+4,   'fillcolor': ROOT.kGreen+4,},
    'Z'   : {'label': 'Z+X',                                  'linecolor': ROOT.kGreen+2,   'fillcolor': ROOT.kGreen+2,},
    'ggZZ': {'label': 'gg#rightarrowZZ, Z#gamma*',            'linecolor': ROOT.kBlue,      'fillcolor': ROOT.kBlue,},
    'ZZ'  : {'label': 'q#bar{q}#rightarrowZZ, Z#gamma*',      'linecolor': ROOT.kAzure+6,   'fillcolor': ROOT.kAzure+6,},
    'HWW' : {'label': 'H(125)#rightarrowWW#rightarrow2l2#nu', 'linecolor': ROOT.kRed+2,     'fillcolor': ROOT.kRed+2,},
    'HZZ' : {'label': 'H(125)#rightarrowZZ#rightarrow4l',     'linecolor': ROOT.kRed+1,     'fillcolor': ROOT.kRed+1,},
    'SIG' : {'label': '2HDM+a',                               'linecolor': ROOT.kMagenta+1,},
    'DATA': {'label': 'Observed',}
}

# setup plotter
plotter = Plotter('MonoHZZ')
plotter.setStyleMap(styleMap)
for bg in backgrounds: plotter.addSampleToStack(bg)
for sig in signals: plotter.addSampleToPlot(sig)
plotter.setDataSample(data)

# plot
plots = {
    'm4l'       : {'hPath': '{sample}_{chan}_hzz_mass',         'xlabel': 'm_{4l} (GeV)',       'binning': range(69,402,3),         'ylabel': 'Events / 3 GeV', 'logx': False, 'logy': False,},
    'm4l_full'  : {'hPath': '{sample}_{chan}_massWindow_mass',  'xlabel': 'm_{4l} (GeV)',       'binning': range(113,135,3),        'ylabel': 'Events / 3 GeV', 'logx': False, 'logy': False, 'blind':blind,},
    'met'       : {'hPath': '{sample}_{chan}_hzz_met',          'xlabel': 'E_{T}^{miss} (GeV)', 'binning': range(0,420,20),         'ylabel': 'Events / 20 GeV','logx': False, 'logy': True, 'ymin':0.1,},
    'met_limit' : {'hPath': '{sample}_{chan}_massWindow_met',   'xlabel': 'E_{T}^{miss} (GeV)', 'binning': [0,25,50,200,500,1000],  'ylabel': 'Events',         'logx': False, 'logy': True, 'ymin':0.1, 'blind':blind,},
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
    # send to plotter
    plotter.plot(hists, f'{plot}', **plots[plot])
