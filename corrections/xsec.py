
# Add STANDARD! cross section values
# (i.e., the ones in the Twiki, the xsec DB, or from LHCHXSWG)
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
# https://cms-gen-dev.cern.ch/xsdb/
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
# if you wish to use a nonstandard cross section for a given sample
# then override it in your processor
# units are PBs
PB = 1.0
FB = 1.0e-3


# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
def get_higgs(prod,decay,lep=False):
    if prod=='ggH': # N3LO
        xsec = 4.852E+01
    elif prod=='VBF':
        xsec = 3.779E+00
    elif prod=='WH':
        xsec = 1.369E+00
    elif prod=='WpH':
        xsec = 8.380E-01
    elif prod=='WmH':
        xsec = 5.313E-01
    elif prod=='WplH':
        xsec = 9.404E-02
    elif prod=='WmlH':
        xsec = 5.967E-02
    elif prod=='ZH':
        xsec = 8.824E-01
    elif prod=='ZllH':
        xsec = 2.977E-02
    elif prod=='ZvvH':
        xsec = 1.773E-01
    elif prod=='ttH':
        xsec = 5.065E-01
    elif prod=='bbH':
        xsec = 4.863E-01
    elif prod=='tH':
        xsec = 7.426E-02
    elif prod=='qqtHb':
        xsec = 2.875E-03
    elif prod=='gbtHW':
        xsec = 1.517E-02
    else:
        raise ValueError(f'No such Higgs production recognized: {prod}')

    if decay=='bb':
        br = 5.809E-01
    elif decay=='tautau':
        br = 6.256E-02
    elif decay=='mumu':
        br = 2.171E-04
    elif decay=='cc':
        br = 2.884E-02
    elif decay=='gg':
        br = 8.180E-02
    elif decay=='gammagamma':
        br = 2.270E-03
    elif decay=='Zgamma':
        br = 1.541E-03
    elif decay=='WW':
        br = 2.152E-01
    elif decay=='ZZ':
        br = 2.641E-02
    else:
        raise ValueError(f'No such Higgs decay recognized: {decay}')
    
    if lep:
        if decay=='Zgamma':
            br += (3.3658e-2*3)
        if decay=='WW':
            br *= (10.86e-2*3)**2
        if decay=='ZZ':
            br *= (3.3658e-2*3)**2

    return xsec * br


xsec = {

    # Data should always be 1
    'SingleMuon': 1.0,
    'DoubleMuon': 1.0,
    'SingleElectron': 1.0,
    'DoubleEG': 1.0,
    'MuonEG': 1.0,
    'EGamma': 1.0,
    'Tau': 1.0,

    # W
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'            : 61526.7,
    'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'             : 61526.7,
    'WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8'                 : 61526.7,
    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'                  : 61526.7,
    # DY
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'   : 18610,
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'    : 18610,
    'DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8'        : 18610,
    'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8'         : 18610,
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'       : 6077.22,
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'        : 6077.22,
    'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8'            : 6077.22,
    'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'             : 6077.22,
    # TT
    'TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8'              : 831.76,
    'TTJets_TuneCUETP8M2T4_13TeV-madgraphMLM-pythia8'               : 831.76,
    'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8'                     : 831.76,
    'TTJets_TuneCP5_13TeV-madgraphMLM-pythia8'                      : 831.76,
    # TTZ
    'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8'          : 0.2529,
    'TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8'     : 0.2529,
    'TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8'               : 0.2529,
    'TTZToLL_M-1to10_TuneCUETP8M1_13TeV-amcatnlo-pythia8'           : 0.05324,
    'TTZToLL_M-1to10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8'      : 0.05324,    
    'TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8'                : 0.05324,
    # WW
    'WWTo2L2Nu_13TeV-powheg'                                        : 12.178,
    'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8'                : 12.178,
    # WZ
    'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8'                    : 4.42965,
    'WZTo3LNu_TuneCP5_13TeV-powheg-pythia8'                         : 4.42965,
    # ZZ
    'ZZTo4L_13TeV_powheg_pythia8'                                   : 1.256,
    'ZZTo4L_13TeV_powheg_pythia8_ext1'                              : 1.256,
    'ZZTo4L_13TeV_powheg_pythia8_TuneCP5'                           : 1.256,
    'ZZTo2L2Nu_13TeV_powheg_pythia8'                                : 0.564,
    'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8'                  : 0.001586,
    'GluGluToContinToZZTo4e_13TeV_TuneCP5_MCFM701_pythia8'          : 0.001586,
    'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8'               : 0.003194,
    'GluGluToContinToZZTo2e2mu_13TeV_TuneCP5_MCFM701_pythia8'       : 0.003194,
    'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8'              : 0.003194,
    'GluGluToContinToZZTo2e2tau_13TeV_TuneCP5_MCFM701_pythia8'      : 0.003194,
    #'GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8'               : 14.93, # confused by this number, missing Z->ll, Z->vv ?
    'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8'                 : 0.001586,
    'GluGluToContinToZZTo4mu_13TeV_TuneCP5_MCFM701_pythia8'         : 0.001586,
    'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8'             : 0.003194,
    'GluGluToContinToZZTo2mu2tau_13TeV_TuneCP5_MCFM701_pythia8'     : 0.003194,
    #'GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8'              : 14.93, # confused by this number, missing Z->ll, Z->vv ?
    'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8'                : 0.001586,
    'GluGluToContinToZZTo4tau_13TeV_TuneCP5_MCFM701_pythia8'        : 0.001586,
    # Higgs 
    # values are N3LO 
    'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8'             : get_higgs('ggH','ZZ',True),
    'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8'            : get_higgs('ggH','ZZ',True),
    'GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUGenV7011_pythia8'   : get_higgs('ggH','ZZ',True),
    'GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8'      : get_higgs('ggH','ZZ',True),
    'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8'               : get_higgs('VBF','ZZ',True),
    'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8'              : get_higgs('VBF','ZZ',True),
    'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV709_pythia8'  : get_higgs('WpH','ZZ',True),
    'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8' : get_higgs('WpH','ZZ',True),
    'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV709_pythia8' : get_higgs('WmH','ZZ',True),
    'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8': get_higgs('WmH','ZZ',True),
    'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV709_pythia8' : get_higgs('ZH','ZZ',True),
    'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8': get_higgs('ZH','ZZ',True),
    'GluGluToZH_HToZZTo4L_M125_13TeV_JHUGenV723_pythia8'                : get_higgs('ZH','ZZ',True),
    'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV709_pythia'           : get_higgs('ttH','ZZ',True),
    'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia'          : get_higgs('ttH','ZZ',True),
    'bbH_HToZZTo4L_M125_13TeV_JHUgenV702_pythia8'                       : get_higgs('bbH','ZZ',True),
    'bbH_HToZZTo4L_M125_13TeV_JHUGenV7011_pythia8'                      : get_higgs('bbH','ZZ',True),
    'tqH_HToZZTo4L_M125_13TeV_JHUgenV702_pythia8'                       : get_higgs('tH','ZZ',True),
    'tqH_HToZZTo4L_M125_13TeV_JHUgenV7011_pythia8'                      : get_higgs('tH','ZZ',True),

}

# for SUSY, manually set to Higgs xsec (ggH) with BR=1e-3
for h in [125,250,500,750,1000]:
    if h == 125:
        ggH = 48.52
    elif h == 250:
        ggH = 10.20
    elif h == 500:
        ggH = 10.20
    elif h == 750:
        ggH = 0.4969
    elif h == 1000:
        ggH = 0.1845
    for a in range(4,51):
        xsec[f'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-{h}_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8'] = ggH * 1e-3
