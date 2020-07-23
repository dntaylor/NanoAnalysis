#!/usr/bin/env python
import os
import json
import subprocess
from utilities import load, dump, get_das
import logging


update = True # requery everything
verbose = False

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# NanoAODv6
#nano_tag = 'Nano25Oct2019'
# NanoAODv7
nano_tag = '02Apr2020'

# data
year_tags = {
    '2016': 'Run2016',
    '2017': 'Run2017',
    '2018': 'Run2018',
}

# datasets
datasets = [
    'SingleMuon',
    'DoubleMuon',
    'SingleElectron',
    'DoubleEG',
    'MuonEG',
    'EGamma', # 2017 instead of SingleElectron/DoubleEG
]

def get_data(update=False,verbose=False):

    fname = 'data'

    result = load(fname)

    for year in year_tags:
        if year not in result: result[year] = {}
        for dataset in datasets:
            query = 'dataset dataset=/{}/{}*{}*/NANOAOD'.format(dataset,year_tags[year],nano_tag)
            samples = get_das(query,verbose=verbose)
            if not samples: continue
            if dataset not in result[year]: result[year][dataset] = {}
            sampleMap = result[year][dataset].get('files',{})
            for sample in samples:
                if not update and sample in sampleMap: continue
                query = 'file dataset={}'.format(sample)
                sampleMap[sample] = get_das(query,verbose=verbose)
    
            result[year][dataset] = {'datasets': samples, 'files': sampleMap}
    
    dump(fname,result)

get_data(update,verbose)

# mc
# note: this doesnt get the "new_pmx" or weirdly named "ext" samples... maybe okay?
# for example: RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_ext_102X_mc2017_realistic_v7
#              RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_new_pmx_102X_mc2017_realistic_v7
year_tags = {
    # NanoAODv6
    #'2016': f'RunIISummer16NanoAODv6-PUMoriond17_{nano_tag}_102X_mcRun2_asymptotic_v7',
    #'2017': f'RunIIFall17NanoAODv6-PU2017_12Apr2018_{nano_tag}_102X_mc2017_realistic_v7',
    #'2018': f'RunIIAutumn18NanoAODv6-{nano_tag}_102X_upgrade2018_realistic_v20',
    # NanoAODv7
    '2016': f'RunIISummer16NanoAODv7-PUMoriond17_Nano{nano_tag}_102X_mcRun2_asymptotic_v8',
    '2017': f'RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano{nano_tag}_102X_mc2017_realistic_v8',
    '2018': f'RunIIAutumn18NanoAODv7-Nano{nano_tag}_102X_upgrade2018_realistic_v21',
}

# datasets (note, tune changes between 2016 and 2017/2018, but not always)
datasets = [
    # W
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
    'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    'WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8',
    # DY
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',
    # WW
    'WWTo2L2Nu_13TeV-powheg',
    'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8',
    # TT
    'TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8',
    'TTJets_TuneCUETP8M2T4_13TeV-madgraphMLM-pythia8',
    'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'TTJets_TuneCP5_13TeV-madgraphMLM-pythia8',
    # TTZ
    'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
    'TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8',
    'TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8'
    'TTZToLL_M-1to10_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
    'TTZToLL_M-1to10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8',
    'TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8',
    # WZ
    'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8',
    'WZTo3LNu_TuneCP5_13TeV-powheg-pythia8',
    # ZZ
    'ZZTo4L_13TeV_powheg_pythia8',
    'ZZTo4L_13TeV_powheg_pythia8_ext1',
    'ZZTo4L_13TeV_powheg_pythia8_TuneCP5',
    'ZZTo2L2Nu_13TeV_powheg_pythia8',
    'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo4e_13TeV_TuneCP5_MCFM701_pythia8',
    'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo2e2mu_13TeV_TuneCP5_MCFM701_pythia8',
    'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo2e2tau_13TeV_TuneCP5_MCFM701_pythia8',
    'GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo4mu_13TeV_TuneCP5_MCFM701_pythia8',
    'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo2mu2tau_13TeV_TuneCP5_MCFM701_pythia8',
    'GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo4tau_13TeV_TuneCP5_MCFM701_pythia8',
    # Higgs v702/709 is 2016, v7011 is 2017/2018
    'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8',
    'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
    'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8',
    'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
    'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV709_pythia8',
    'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',
    'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV709_pythia8',
    'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',
    'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV709_pythia8',
    'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8',
    'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV709_pythia',
    'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia',
    'bbH_HToZZTo4L_M125_13TeV_JHUgenV702_pythia8',
    'bbH_HToZZTo4L_M125_13TeV_JHUGenV7011_pythia8',
    'tqH_HToZZTo4L_M125_13TeV_JHUgenV702_pythia8',
    'tqH_HToZZTo4L_M125_13TeV_JHUgenV7011_pythia8',
    # HAA
    'SUSY*HToAA*AToMuMu*AToTauTau*',
]


def get_mc(update=False,verbose=False):

    fname = 'mc'

    result = load(fname)

    for year in year_tags:
        if year not in result: result[year] = {}
        for dataset in datasets:
            query = 'dataset dataset=/{}/{}*/NANOAODSIM'.format(dataset,year_tags[year])
            samples = get_das(query,verbose=verbose)
            if not samples: continue
            thesedatasets = set(s.split('/')[1] for s in samples)
            for thisdataset in thesedatasets:
                if thisdataset not in result[year]: result[year][thisdataset] = {}
                sampleMap = result[year][thisdataset].get('files',{})
                goodsamples = []
                for sample in samples:
                    if not update and sample in sampleMap: continue
                    if 'Validation error' in sample: continue
                    if sample.split('/')[1]!=thisdataset: continue
                    query = 'file dataset={}'.format(sample)
                    sampleMap[sample] = get_das(query,verbose=verbose)
                    goodsamples += [sample]
    
                result[year][thisdataset] = {'datasets': goodsamples, 'files': sampleMap}
    
    dump(fname,result)

get_mc(update,verbose)
