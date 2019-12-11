import os
import json
import subprocess


reload = True # requery everything

def run_command(command):
    return subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

def get_das(query):
    command = 'dasgoclient -query "{}"'.format(query)
    return [line.decode('utf-8') for line in run_command(command).split()]

def load(fname):
    jname = '{}.json'.format(fname)
    if not os.path.exists(jname): return {}
    with open(jname) as f:
        return json.load(f)

def dump(fname,content):
    jname = '{}.json'.format(fname)
    with open(jname,'w') as f:
        json.dump(content,f, indent=4, sort_keys=True)

# NanoAODv6
nano_tag = 'Nano25Oct2019'

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

def get_data(reload=False):

    fname = 'data'

    result = load(fname)

    for year in year_tags:
        if year not in result: result[year] = {}
        for dataset in datasets:
            query = 'dataset dataset=/{}/{}*{}*/NANOAOD'.format(dataset,year_tags[year],nano_tag)
            samples = get_das(query)
            if not samples: continue
            if dataset not in result[year]: result[year][dataset] = {}
            sampleMap = result[year][dataset].get('files',{})
            for sample in samples:
                if not reload and sample in sampleMap: continue
                query = 'file dataset={}'.format(sample)
                sampleMap[sample] = get_das(query)
    
            result[year][dataset] = {'datasets': samples, 'files': sampleMap}
    
    dump(fname,result)

get_data(reload)

# mc
year_tags = {
    '2016': 'RunIISummer16NanoAODv6',
    '2017': 'RunIIFall17NanoAODv6',
    '2018': 'RunIIAutumn18NanoAODv6',
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
    # Higgs v6/709 is 2016, v7011 is 2017/2018
    'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
    'GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8',
    'GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUGenV7011_pythia8',
    'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8',
    'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
    'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8', # no 2016?
    'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8', # no 2016?
    'bbH_HToZZTo4L_M125_13TeV_JHUGenV7011_pythia8', # no 2016?
    'GluGluToZH_HToZZTo4L_M125_13TeV_JHUGenV723_pythia8', # 2016/2017/2018
]


def get_mc(reload=False):

    fname = 'mc'

    result = load(fname)

    for year in year_tags:
        if year not in result: result[year] = {}
        for dataset in datasets:
            query = 'dataset dataset=/{}/{}*{}*/NANOAODSIM'.format(dataset,year_tags[year],nano_tag)
            samples = get_das(query)
            if not samples: continue
            if dataset not in result[year]: result[year][dataset] = {}
            sampleMap = result[year][dataset].get('files',{})
            for sample in samples:
                if not reload and sample in sampleMap: continue
                query = 'file dataset={}'.format(sample)
                sampleMap[sample] = get_das(query)
    
            result[year][dataset] = {'datasets': samples, 'files': sampleMap}
    
    dump(fname,result)

get_mc(reload)
