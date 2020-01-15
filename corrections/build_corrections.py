#!/usr/bin/env python
import numpy as np
import uproot
from coffea import hist, lookup_tools
from coffea.util import load, save
from xsec import xsec

def save_corrections(year):
    corrections = {}

    # golden json
    if year == '2016':
        corrections['golden'] = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'
    if year == '2017':
        corrections['golden'] = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
    if year == '2018':
        corrections['golden'] = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'

    
    # cross sections
    corrections['xsec'] = xsec
    # manually add the test samples
    corrections['xsec']['DY'] = 6077.22
    corrections['xsec']['HZZ'] = 43.92 * 2.64e-02 * (3.3658e-2*3)**2
    corrections['xsec']['DoubleMuon'] = 1.
    
    extractor = lookup_tools.extractor()
    # electron
    # POG
    if year == '2016':
        extractor.add_weight_sets([
            'electron_id_Veto_ * data/scalefactors/electron/2016/2016_ElectronWPVeto_Fall17V2.root',
            'electron_id_Loose_ * data/scalefactors/electron/2016/2016LegacyReReco_ElectronLoose_Fall17V2.root',
            'electron_id_Medium_ * data/scalefactors/electron/2016/2016LegacyReReco_ElectronMedium_Fall17V2.root',
            'electron_id_Tight_ * data/scalefactors/electron/2016/2016LegacyReReco_ElectronTight_Fall17V2.root',
            'electron_id_MVA80_ * data/scalefactors/electron/2016/2016LegacyReReco_ElectronMVA80_Fall17V2.root',
            'electron_id_MVA90_ * data/scalefactors/electron/2016/2016LegacyReReco_ElectronMVA90_Fall17V2.root',
            'electron_id_MVA80noiso_ * data/scalefactors/electron/2016/2016LegacyReReco_ElectronMVA80noiso_Fall17V2.root',
            'electron_id_MVA90noiso_ * data/scalefactors/electron/2016/2016LegacyReReco_ElectronMVA90noiso_Fall17V2.root',
        ])
    elif year == '2017':
        extractor.add_weight_sets([
            'electron_id_Veto_ * data/scalefactors/electron/2017/2017_ElectronWPVeto_Fall17V2.root',
            'electron_id_Loose_ * data/scalefactors/electron/2017/2017_ElectronLoose.root',
            'electron_id_Medium_ * data/scalefactors/electron/2017/2017_ElectronMedium.root',
            'electron_id_Tight_ * data/scalefactors/electron/2017/2017_ElectronTight.root',
            'electron_id_MVA80_ * data/scalefactors/electron/2017/2017_ElectronMVA80.root',
            'electron_id_MVA90_ * data/scalefactors/electron/2017/2017_ElectronMVA90.root',
            'electron_id_MVA80noiso_ * data/scalefactors/electron/2017/2017_ElectronMVA80noiso.root',
            'electron_id_MVA90noiso_ * data/scalefactors/electron/2017/2017_ElectronMVA90noiso.root',
        ])
    elif year == '2018':
        extractor.add_weight_sets([
            'electron_id_Veto_ * data/scalefactors/electron/2018/2018_ElectronWPVeto_Fall17V2.root',
            'electron_id_Loose_ * data/scalefactors/electron/2018/2018_ElectronLoose.root',
            'electron_id_Medium_ * data/scalefactors/electron/2018/2018_ElectronMedium.root',
            'electron_id_Tight_ * data/scalefactors/electron/2018/2018_ElectronTight.root',
            'electron_id_MVA80_ * data/scalefactors/electron/2018/2018_ElectronMVA80.root',
            'electron_id_MVA90_ * data/scalefactors/electron/2018/2018_ElectronMVA90.root',
            'electron_id_MVA80noiso_ * data/scalefactors/electron/2018/2018_ElectronMVA80noiso.root',
            'electron_id_MVA90noiso_ * data/scalefactors/electron/2018/2018_ElectronMVA90noiso.root',
        ])

    # HZZ
    extractor.add_weight_sets([
        # electron reco
        f'electron_reco_ * data/scalefactors/electron/{year}/Ele_Reco_{year}.root',
        # electron hzz id
        f'electron_hzz_id_nogap_ * data/scalefactors/electron/{year}/ElectronSF_Legacy_{year}_NoGap.root',
        f'electron_hzz_id_gap_ * data/scalefactors/electron/{year}/ElectronSF_Legacy_{year}_Gap.root',
    ])

    # muon
    # POG
    if year == '2016':
        extractor.add_weight_sets([
            # id
            'muon_id_ * data/scalefactors/muon/2016/EfficienciesStudies_2016_legacy_rereco_rootfiles_RunBCDEF_SF_ID.root',
            'muon_id_2_ * data/scalefactors/muon/2016/EfficienciesStudies_2016_legacy_rereco_rootfiles_RunGH_SF_ID.root',
            # iso
            'muon_iso_ * data/scalefactors/muon/2016/EfficienciesStudies_2016_legacy_rereco_rootfiles_RunBCDEF_SF_ISO.root',
            'muon_iso_2_ * data/scalefactors/muon/2016/EfficienciesStudies_2016_legacy_rereco_rootfiles_RunGH_SF_ISO.root',
            # jpsi
            'muon_id_jpsi_ * data/scalefactors/muon/2016/EfficienciesStudies_2016_legacy_rereco_Jpsi_rootfiles_RunBCDEF_SF_ID.root',
            'muon_id_jpsi_2_ * data/scalefactors/muon/2016/EfficienciesStudies_2016_legacy_rereco_Jpsi_rootfiles_RunGH_SF_ID.root',
        ])
    elif year == '2017':
        extractor.add_weight_sets([
            # id
            'muon_id_ * data/scalefactors/muon/2017/EfficienciesStudies_2017_rootfiles_RunBCDEF_SF_ID.root',
            # iso
            'muon_iso_ * data/scalefactors/muon/2017/EfficienciesStudies_2017_rootfiles_RunBCDEF_SF_ISO.root',
            # jpsi
            'muon_id_jpsi_ * data/scalefactors/muon/2017/EfficienciesStudies_2017_rootfiles_RunBCDEF_SF_ID_JPsi.root',
        ])
    elif year == '2018':
        extractor.add_weight_sets([
            # id
            'muon_id_ * data/scalefactors/muon/2018/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root',
            # iso
            'muon_iso_ * data/scalefactors/muon/2018/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ISO.root',
            # jpsi
            'muon_id_jpsi_ * data/scalefactors/muon/2018/EfficienciesStudies_2018_Jpsi_rootfiles_RunABCD_SF_ID.root',
        ])

    extractor.finalize()
    evaluator = extractor.make_evaluator()
    
    # EGamma POG corrections
    idnums = [
        'Veto',
        'Loose',
        'Medium',
        'Tight',
        'MVA80',
        'MVA90',
        'MVA80noiso',
        'MVA90noiso',
    ]

    for idnum in idnums:
        corrections[f'electron_id_{idnum}'] = evaluator[f'electron_id_{idnum}_EGamma_SF2D']

    # HZZ corrections
    corrections['electron_reco'] = evaluator['electron_reco_EGamma_SF2D']
    corrections['electron_hzz_id_nogap'] = evaluator['electron_hzz_id_nogap_EGamma_SF2D']
    corrections['electron_hzz_id_gap'] = evaluator['electron_hzz_id_gap_EGamma_SF2D']

    # Muon POG corrections
    if year == '2016':
        idnums = [
            'LooseID',
            'MediumID',
            'TightID',
            'HighPtID',
        ]
        iddenom = 'genTracks'
        effvars = 'eta_pt'
        highpt_effvars = 'eta_pair_newTuneP_probe_pt'
        iso_num_denoms = [
            ('LooseRelTkIso', 'HighPtIDandIPCut'),
            ('TightRelIso', 'MediumID'),
            ('TightRelIso', 'TightIDandIPCut'),
            ('LooseRelIso', 'LooseID'),
            ('LooseRelIso', 'MediumID'),
            ('LooseRelIso', 'TightIDandIPCut'),
        ]
        jpsinums = [
            'LooseID',
            'MediumID',
            'TightID',
            'SoftID',
        ]
        jpsidenom = 'genTracks'
        jpsieffvars = 'pt_abseta'
    elif year in ['2017','2018']:
        idnums = [
            'LooseID',
            'MediumID',
            'MediumPromptID',
            'TightID',
            'SoftID',
            'HighPtID',
            'TrkHighPtID',
        ]
        iddenom = 'genTracks'
        if year == '2018':
            iddenom = 'TrackerMuons'
        effvars = 'pt_abseta'
        highpt_effvars = 'pair_newTuneP_probe_pt_abseta'
        iso_num_denoms = [
            ('LooseRelTkIso', 'HighPtIDandIPCut'),
            ('LooseRelTkIso', 'TrkHighPtID'),
            ('TightRelTkIso', 'HighPtIDandIPCut'),
            ('TightRelTkIso', 'TrkHighPtID'),
            ('TightRelIso', 'MediumID'),
            ('TightRelIso', 'TightIDandIPCut'),
            ('LooseRelIso', 'LooseID'),
            ('LooseRelIso', 'MediumID'),
            ('LooseRelIso', 'TightIDandIPCut'),
        ]
        jpsinums = [
            'LooseID',
            'MediumID',
            'TightID',
            'SoftID',
        ]
        jpsidenom = 'genTracks'
        jpsieffvars = 'pt_abseta'

    lumi2016_BCDEF = 19.721 / (16.146 + 19.721)
    lumi2016_GH    = 16.146 / (16.146 + 19.721)

    for idnum in idnums:
        histkey = f'NUM_{idnum}_DEN_{iddenom}_{effvars}'
        if idnum in ['HighPtID', 'TrkHighPtID']:
            histkey = f'NUM_{idnum}_DEN_{iddenom}_{highpt_effvars}'
        corrections[f'muon_id_{idnum}'] = evaluator[f'muon_id_{histkey}']
        if year == '2016':
            corrections[f'muon_id_{idnum}']._values *= lumi2016_BCDEF
            corrections[f'muon_id_{idnum}']._values += evaluator[f'muon_id_2_{histkey}']._values * lumi2016_GH

    for isonum, isodenom in iso_num_denoms:
        histkey = f'NUM_{isonum}_DEN_{isodenom}_{effvars}'
        if isodenom in ['HighPtIDandIPCut', 'TrkHighPtID']:
            histkey = f'NUM_{isonum}_DEN_{isodenom}_{highpt_effvars}'
        corrections[f'muon_iso_{isonum}_{isodenom}'] = evaluator[f'muon_iso_{histkey}']
        if year == '2016':
            corrections[f'muon_iso_{isonum}_{isodenom}']._values *= lumi2016_BCDEF
            corrections[f'muon_iso_{isonum}_{isodenom}']._values += evaluator[f'muon_iso_2_{histkey}']._values * lumi2016_GH

    for jpsinum in jpsinums:
        histkey = f'NUM_{jpsinum}_DEN_{jpsidenom}_{jpsieffvars}'
        corrections[f'muon_id_jpsi_{jpsinum}'] = evaluator[f'muon_id_jpsi_{histkey}']
        if year == '2016':
            corrections[f'muon_id_jpsi_{jpsinum}']._values *= lumi2016_BCDEF
            corrections[f'muon_id_jpsi_{jpsinum}']._values += evaluator[f'muon_id_jpsi_2_{histkey}']._values * lumi2016_GH




    
    # pileup
    # from NanoAOD tools
    # 2016 has a bug
    #with uproot.open(f'data/pileup/dataPileup{year}.root') as f:
    #    norm = lambda x: x/x.sum()
    #    edges = f['pileup'].edges
    #    dataPileup = norm(f['pileup'].values)
    #    dataPileupUp = norm(f['pileup_plus'].values)
    #    dataPileupDown = norm(f['pileup_minus'].values)
    #with uproot.open(f'data/pileup/mcPileup{year}.root') as f:
    #    mcPileup = f['pu_mc'].values
    #def zeropad(a,n):
    #    _a = np.zeros(n)
    #    _a[:len(a)] = a
    #    return _a
    #nmax = max(len(dataPileup),len(mcPileup))
    #dataPileup = zeropad(dataPileup,nmax)
    #mcPileup = zeropad(mcPileup,nmax)
    #mask = (mcPileup>0)
    #pileupRatio = dataPileup.copy()
    #pileupRatioUp = dataPileupUp.copy()
    #pileupRatioDown = dataPileupDown.copy()
    #pileupRatio[mask] /= mcPileup[mask]
    #pileupRatioUp[mask] /= mcPileup[mask]
    #pileupRatioDown[mask] /= mcPileup[mask]
    # from HZZ
    with uproot.open(f'data/pileup/pu_weights_{year}.root') as f:
        edges = f['weights'].edges
        pileupRatio = f['weights'].values
        pileupRatioUp = f['weights_varUp'].values
        pileupRatioDown = f['weights_varDn'].values
    
    corrections['pileupWeight'] = lookup_tools.dense_lookup.dense_lookup(pileupRatio, edges)
    corrections['pileupWeightUp'] = lookup_tools.dense_lookup.dense_lookup(pileupRatioUp, edges)
    corrections['pileupWeightDown'] = lookup_tools.dense_lookup.dense_lookup(pileupRatioDown, edges)


    # rochester correction
    tag = 'roccor.Run2.v3'
    fname = f'data/rochester/{tag}/RoccoR{year}.txt'
    corrections['rochester_data'] = lookup_tools.txt_converters.convert_rochester_file(fname,loaduncs=True)

    
    save(corrections, f'corrections/corrections_{year}.coffea')

for year in ['2016','2017','2018']:
    save_corrections(year)
