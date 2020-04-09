import os

# this defines sample priorities
# sample names have changed through the years

class Sample:
    '''A class that will return a single sample given a list of equivalent samples
    listed in priority order. If no sample is found, returns None.'''

    HISTPATH = 'hists/{base}/{year}/{sample}.root'

    def __init__(self,base,year):
        self._base = base
        self._year = year

    def __call__(self,*samples):
        for sample in samples:
            spath = self.HISTPATH.format(
                base = self._base,
                year = self._year,
                sample = sample,
            )
            if os.path.exists(spath): return sample
        return None
        

# something to consider implementing:
#  if TTTo2L and TTTo1L1Nu available, return both, else fall back to TTJets



def get_sample_list(base,year,sample):
    '''Returns a list of samples for the given sample string'''

    picker = Sample(base,year)

    if sample == 'DATA':
        samples = []
        if year in ['2016','2017']:
            samples += [picker('DoubleEG')]
            samples += [picker('SingleElectron')]
        elif year in ['2018']:
            samples += [picker('EGamma')]
        samples += [picker('DoubleMuon')]
        samples += [picker('SingleMuon')]
        samples += [picker('MuonEG')]
        return [s for s in samples if s is not None]

    Z = [
        picker(
            'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
            'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
            'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',
            'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        ),
        picker(
            'DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
            'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
            'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8',
            'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        ),
    ]
    if sample == 'Z': return Z

    TT = [
        picker(
            'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8',
            'TTJets_TuneCP5_13TeV-madgraphMLM-pythia8',
            'TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8',
        )
    ]
    if sample == 'TT': return TT

    # WZ
    WZ = [
        picker(
            'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8',
            'WZTo3LNu_TuneCP5_13TeV-powheg-pythia8',
        )
    ]
    if sample == 'WZ': return WZ

    # ZZ
    qqZZ = [
        picker(
            'ZZTo4L_13TeV_powheg_pythia8_TuneCP5',
            'ZZTo4L_13TeV_powheg_pythia8',
        )
    ]

    ggZZ = [
        picker(
            'GluGluToContinToZZTo2e2mu_13TeV_TuneCP5_MCFM701_pythia8',
            'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8',
        ),
        picker(
            'GluGluToContinToZZTo2e2tau_13TeV_TuneCP5_MCFM701_pythia8',
            'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8',
        ),
        picker(
            'GluGluToContinToZZTo2mu2tau_13TeV_TuneCP5_MCFM701_pythia8',
            'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8',
        ),
        picker(
            'GluGluToContinToZZTo4e_13TeV_TuneCP5_MCFM701_pythia8',
            'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8',
        ),
        picker(
            'GluGluToContinToZZTo4mu_13TeV_TuneCP5_MCFM701_pythia8',
            'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8',
        ),
        picker(
            'GluGluToContinToZZTo4tau_13TeV_TuneCP5_MCFM701_pythia8',
            'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8',
        ),
    ]
    if sample == 'qqZZ': return qqZZ
    if sample == 'ggZZ': return ggZZ
    if sample == 'ZZ': return qqZZ + ggZZ

    HZZ = [
        picker(
            'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
            'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8',
        ),
        picker(
            'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
            'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8',
        ),
        picker(
            'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',
            'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV709_pythia8',
        ),
        picker(
            'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',
            'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV709_pythia8',
        ),
        picker(
            'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8',
            'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV709_pythia8',
        ),
        picker(
            'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8',
            'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV709_pythia8',
        ),
        picker(
            'bbH_HToZZTo4L_M125_13TeV_JHUGenV7011_pythia8',
            'bbH_HToZZTo4L_M125_13TeV_JHUGenV702_pythia8',
        ),
        picker(
            'tqH_HToZZTo4L_M125_13TeV_JHUgenV7011_pythia8',
            'tqH_HToZZTo4L_M125_13TeV_JHUgenV702_pythia8',
        ),
    ]
    if sample == 'HZZ': return HZZ

    for h in [125,250,500,750,1000]:
        for a in range(1,61):
            if sample == f'haaH{h}A{a}':
                return [
                    picker(
                        f'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-{h}_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8',
                        f'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8', # h 125 2016 only
                    )
                ]
