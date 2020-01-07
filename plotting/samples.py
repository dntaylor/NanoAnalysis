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
        return samples

    if sample == 'Z':
        return [
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

    if sample == 'TT':
        return [
            picker(
                'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8',
                'TTJets_TuneCP5_13TeV-madgraphMLM-pythia8',
                'TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8',
            )
        ]
