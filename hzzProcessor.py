#!/usr/bin/env python
from __future__ import print_function, division
import os
import sys
import logging
import argparse
import numpy as np
import uproot
import uproot_methods
from coffea import hist, processor, lookup_tools
from coffea.util import load, save
from coffea.analysis_objects import JaggedCandidateArray
from awkward import JaggedArray, IndexedArray

ZMASS = 91.1876

logger = logging.getLogger("MonoHZZProcessor")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class HZZProcessor(processor.ProcessorABC):
    # will sync with
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4lRunIILegacy
    # to the best of NanoAOD ability
    def __init__(self,year='2018'):
        self._year = year

        self._corrections = {}

        # cross sections
        # TODO: load
        self._corrections['xsec'] = {}
        # manually add the test samples
        self._corrections['xsec']['DY'] = 6077.22
        self._corrections['xsec']['HZZ'] = 43.92 * 2.64e-02	* (3.3658e-2*3)**2
        self._corrections['xsec']['DoubleMuon'] = 1.

        extractor = lookup_tools.extractor()
        # electron
        extractor.add_weight_sets([
            # electron reco
            f'electron_reco_ * data/scalefactors/electron/Ele_Reco_{self._year}.root',
            # electron hzz id
            f'electron_hzz_id_nogap_ * data/scalefactors/electron/ElectronSF_Legacy_{self._year}_NoGap.root',
            f'electron_hzz_id_gap_ * data/scalefactors/electron/ElectronSF_Legacy_{self._year}_Gap.root',
        ])
        extractor.finalize()
        evaluator = extractor.make_evaluator()

        self._corrections['electron_reco'] = evaluator['electron_reco_EGamma_SF2D']
        self._corrections['electron_hzz_id_nogap'] = evaluator['electron_hzz_id_nogap_EGamma_SF2D']
        self._corrections['electron_hzz_id_gap'] = evaluator['electron_hzz_id_gap_EGamma_SF2D']

        # pileup
        with uproot.open(f'data/pileup/dataPileup{self._year}.root') as f:
            norm = lambda x: x/x.sum()
            edges = f['pileup'].edges
            dataPileup = norm(f['pileup'].values)
            dataPileupUp = norm(f['pileup_plus'].values)
            dataPileupDown = norm(f['pileup_minus'].values)
        with uproot.open(f'data/pileup/mcPileup{self._year}.root') as f:
            mcPileup = f['pu_mc'].values
        mask = (mcPileup>0)
        pileupRatio = dataPileup.copy()
        pileupRatioUp = dataPileupUp.copy()
        pileupRatioDown = dataPileupDown.copy()
        pileupRatio[mask] /= mcPileup[mask]
        pileupRatioUp[mask] /= mcPileup[mask]
        pileupRatioDown[mask] /= mcPileup[mask]

        self._corrections[f'pileupWeight{self._year}'] = lookup_tools.dense_lookup.dense_lookup(pileupRatio, edges)
        self._corrections[f'pileupWeight{self._year}Up'] = lookup_tools.dense_lookup.dense_lookup(pileupRatioUp, edges)
        self._corrections[f'pileupWeight{self._year}Down'] = lookup_tools.dense_lookup.dense_lookup(pileupRatioDown, edges)

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        channel_axis = hist.Cat("channel", "Channel")
        mass_axis = hist.Bin("mass", r"$m_{4\ell}$ [GeV]", 600, 0.25, 300)
        zmass_axis = hist.Bin("mass", r"$m_{2\ell}$ [GeV]", 240, 0, 120)
        pt_axis = hist.Bin("pt", r"$p_{T,\ell}$ [GeV]", 3000, 0.25, 300)
        met_axis = hist.Bin("pt", r"$E_{T}^{miss}$ [GeV]", 3000, 0, 3000)

        hist.Hist.DEFAULT_DTYPE = 'f'  # save some space by keeping float bin counts instead of double
        self._accumulator = processor.dict_accumulator({
            'mass': hist.Hist("Counts", dataset_axis, channel_axis, mass_axis),
            'z1mass': hist.Hist("Counts", dataset_axis, channel_axis, zmass_axis),
            'z2mass': hist.Hist("Counts", dataset_axis, channel_axis, zmass_axis),
            #'met': hist.Hist("Counts", dataset_axis, channel_axis, met_axis),
            #'pt_lead': hist.Hist("Counts", dataset_axis, channel_axis, pt_axis),
            'cutflow': processor.defaultdict_accumulator(int),
            'sumw': processor.defaultdict_accumulator(int),
        })
        

    @property
    def accumulator(self):
        return self._accumulator

    def _add_muon_id(self, muons):
        # note: input muons must pass
        # slimmedMuons and (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))
        ptCut = (muons.pt>5)
        etaCut = (abs(muons.eta)<2.4)
        dxyCut = (abs(muons.dxy)<0.5)
        dzCut = (abs(muons.dz)<1)
        idCut = ((muons.isGlobal) | ((muons.isTracker) & (muons.nStations>0))) # TODO following not in nanoaod: & (muons.muonBestTrackType != 2)
        sipCut = (abs(muons.sip3d)<4)

        hzzLooseNoIso = (ptCut & etaCut & dxyCut & dzCut & idCut & sipCut)

        isoCut = (muons.pfRelIso03_all<0.35)

        hzzLoose = (hzzLooseNoIso & isoCut)

        tightCut = ((muons.isPFcand) | ((muons.highPtId>0) & (muons.pt>200)))

        hzzTight = (hzzLoose & tightCut)

        # TODO: ghost cleaning

        muons['hzzLooseNoIso'] = hzzLooseNoIso
        muons['hzzLoose'] = hzzLoose
        muons['hzzTight'] = hzzTight
        
    def _add_electron_id(self, electrons):
        # note: input electrons must pass
        # slimmedElectrons and (pt > 5)
        ptCut = (electrons.pt>7)
        etaCut = (abs(electrons.eta)<2.5)
        dxyCut = (abs(electrons.dxy)<0.5)
        dzCut = (abs(electrons.dz)<1)
        sipCut = (abs(electrons.sip3d)<4)

        hzzLooseNoIso = (ptCut & etaCut & dxyCut & dzCut & sipCut)

        isoCut = (electrons.pfRelIso03_all<0.35)

        hzzLoose = (hzzLooseNoIso & isoCut)

        # TODO: check
        pt = electrons.pt
        eta = abs(electrons.etaSC)
        mva = electrons.mva
        # these are the 2017 trainings only, because that is what is available
        tightBin00 = ((pt<=10) & (eta<0.8)                  & (mva>0.8955937602))
        tightBin01 = ((pt<=10) & ((eta>=0.8) & (eta<1.479)) & (mva>0.91106464032))
        tightBin02 = ((pt<=10) & (eta>=1.479)               & (mva>0.94067753025))
        tightBin10 = ((pt>10 ) & (eta<0.8)                  & (mva>0.04240620843))
        tightBin11 = ((pt>10 ) & ((eta>=0.8) & (eta<1.479)) & (mva>0.0047338429))
        tightBin12 = ((pt>10 ) & (eta>=1.479)               & (mva>-0.60423293572))
        tightCut = (tightBin00 | tightBin01 | tightBin02 | tightBin10 | tightBin11 | tightBin12)

        hzzTight = (hzzLoose & tightCut)

        electrons['hzzLooseNoIso'] = hzzLooseNoIso
        electrons['hzzLoose'] = hzzLoose
        electrons['hzzTight'] = hzzTight
        
    def _add_photon_id(self, photons, muons, electrons):
        ptCut = (photons.pt>2)
        etaCut = (abs(photons.eta)<2.4)
        isoCut = (photons.relIso03<1.8)

        hzzPreselection = (ptCut & etaCut & isoCut)

        # matching muons is all that is available now
        drOverEtCut = (photons.dROverEt2<0.012)
        drMuonCut = (photons.p4.delta_r(muons.p4[photons.muonIdx])<0.5)

        hzzMatchMuon = (hzzPreselection & drOverEtCut & drMuonCut)
        hzzMatchElectron = (hzzPreselection & (np.zeros_like(hzzPreselection)))

        photons['hzzPreselection'] = hzzPreselection
        photons['hzzMatchMuon'] = hzzMatchMuon
        photons['hzzMatchElectron'] = hzzMatchElectron
        

    def _add_trigger(self,df):
        dataset = df['dataset']

        triggerPaths = {}
        # DoubleMuon
        if self._year=='2016':
            triggerPaths['DoubleMuon'] = [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                "HLT_TripleMu_12_10_5",
            ]
        elif self._year=='2017':
            triggerPaths['DoubleMuon'] = [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                "HLT_TripleMu_10_5_5_DZ",
                "HLT_TripleMu_12_10_5",
            ]
        elif self._year=='2018':
            triggerPaths['DoubleMuon'] = [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "HLT_TripleMu_10_5_5_DZ",
                "HLT_TripleMu_12_10_5",
            ]

        # DoubleEG
        if self._year=='2016':
            triggerPaths['DoubleEG'] = [
                "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ]

        # MuonEG
        if self._year=='2016':
            triggerPaths['MuonEG'] = [
                "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
            ]
        elif self._year=='2017':
            triggerPaths['MuonEG'] = [
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
            ]
        elif self._year=='2018':
            triggerPaths['MuonEG'] = [
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
            ]

        # EGamma
        if self._year=='2017':
            triggerPaths['EGamma'] = [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_DoubleEle33_CaloIdL_MW",
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
                "HLT_Ele35_WPTight_Gsf",
                "HLT_Ele38_WPTight_Gsf",
                "HLT_Ele40_WPTight_Gsf",
            ]
        elif self._year=='2018':
            triggerPaths['EGamma'] = [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_DoubleEle25_CaloIdL_MW",
                "HLT_Ele32_WPTight_Gsf",
            ]


        # SingleMuon
        if self._year=='2016':
            triggerPaths['SingleMuon'] = [
                "HLT_IsoMu20",
                "HLT_IsoTkMu20",
                "HLT_IsoMu22",
                "HLT_IsoTkMu22",
                "HLT_IsoMu24",
                "HLT_IsoTkMu24",
            ]
        elif self._year=='2017':
            triggerPaths['SingleMuon'] = [
                "HLT_IsoMu27",
            ]
        elif self._year=='2018':
            triggerPaths['SingleMuon'] = [
                "HLT_IsoMu24",
            ]

        # SingleElectron
        if self._year=='2016':
            triggerPaths['SingleElectron'] = [
                "HLT_Ele25_eta2p1_WPTight_Gsf",
                "HLT_Ele27_WPTight_Gsf",
                "HLT_Ele27_eta2p1_WPLoose_Gsf",
                "HLT_Ele32_eta2p1_WPTight_Gsf",
            ]


        # Define priority
        # To avoid double counting in data, for a given dataset
        # all lower and current datasets triggers are accepted
        # and all higher datasets triggers are vetoed
        # in MC, all triggers are accepted
        if self._year=='2016':
            triggerPriority = [
                'DoubleMuon',
                'DoubleEG',
                'MuonEG',
                'SingleMuon',
                'SingleElectron',
            ]
        else:
            triggerPriority = [
                'DoubleMuon',
                'EGamma',
                'MuonEG',
                'SingleMuon',
            ]

        triggersToAccept = []
        triggersToVeto = []
        accept = dataset not in triggerPriority # accept MC, reject data until you reach the correct dataset
        for d in triggerPriority:
            if d==dataset: accept = True # start accepting triggers
            for p in triggerPaths[d]:
                if accept:
                    triggersToAccept += [p]
                else:
                    triggersToVeto += [p]

        result = np.zeros_like(df['event'],dtype=bool)
        for p in triggersToAccept:
            result = ((result) | (df[p]))
        for p in triggersToVeto:
            result = ((result) & (~df[p]))

        df['passHLT'] = result
                


    def process(self, df):
        logging.debug('starting process')
        output = self.accumulator.identity()

        dataset = df['dataset']
        self._isData = dataset in ['SingleMuon','DoubleMuon','SingleElectron','DoubleEG','EGamma','MuonEG']


        selection = processor.PackedSelection()

        # TODO: instead of cutflow, use processor.PackedSelection
        output['cutflow']['all events'] += df['event'].size

        logging.debug('adding trigger')
        self._add_trigger(df)

        passHLT = df['passHLT']
        selection.add('trigger',passHLT)
        output['cutflow']['pass trigger'] += passHLT.sum()

        
        # require one good vertex
        logging.debug('checking vertices')
        passGoodVertex = (df['PV_npvsGood']>0)
        output['cutflow']['good vertex'] += passGoodVertex.sum()
        selection.add('goodVertex',passGoodVertex)


        logging.debug('building muons')
        muons = JaggedCandidateArray.candidatesfromcounts(
            df['nMuon'],
            pt=df['Muon_pt'],
            eta=df['Muon_eta'],
            phi=df['Muon_phi'],
            mass=df['Muon_mass'],
            charge=df['Muon_charge'],
            dxy=df['Muon_dxy'],
            dz=df['Muon_dz'],
            isGlobal=df['Muon_isGlobal'],
            isTracker=df['Muon_isTracker'],
            nStations=df['Muon_nStations'],
            sip3d=df['Muon_sip3d'],
            pfRelIso03_all=df['Muon_pfRelIso03_all'],
            isPFcand=df['Muon_isPFcand'],
            highPtId=df['Muon_highPtId'],
            fsrPhotonIdx=df['Muon_fsrPhotonIdx'],
            pdgId=df['Muon_pdgId'],
        )

        logging.debug('building electrons')
        if self._year=='2016':
            mvaDisc = "mvaSummer16IdIso"
        elif self._year=='2017':
            mvaDisc = "mvaFall17V2Iso"
        else:
            mvaDisc = "mvaAutumn18IdIso"
        mvaDisc = "mvaFall17V2Iso" # only one available
        electrons = JaggedCandidateArray.candidatesfromcounts(
            df['nElectron'],
            pt=df['Electron_pt'],
            eta=df['Electron_eta'],
            phi=df['Electron_phi'],
            mass=df['Electron_mass'],
            charge=df['Electron_charge'],
            dxy=df['Electron_dxy'],
            dz=df['Electron_dz'],
            sip3d=df['Electron_sip3d'],
            pfRelIso03_all=df['Electron_pfRelIso03_all'],
            deltaEtaSC=df['Electron_deltaEtaSC'],
            etaSC=df['Electron_eta']+df['Electron_deltaEtaSC'],
            mva=df[f'Electron_{mvaDisc}'],
            pdgId=df['Electron_pdgId'],
        )

        logging.debug('building fsr photons')
        photons = JaggedCandidateArray.candidatesfromcounts(
            df['nFsrPhoton'],
            pt=df['FsrPhoton_pt'],
            eta=df['FsrPhoton_eta'],
            phi=df['FsrPhoton_phi'],
            mass=np.zeros_like(df['FsrPhoton_pt']),
            dROverEt2=df['FsrPhoton_dROverEt2'],
            relIso03=df['FsrPhoton_relIso03'],
            muonIdx=df['FsrPhoton_muonIdx'],
        )


        logging.debug('adding muon id')
        self._add_muon_id(muons)
        logging.debug('adding electron id')
        self._add_electron_id(electrons)
        logging.debug('add fsr photon id')
        self._add_photon_id(photons,muons,electrons)

        # TODO, select loose cands here, tight cands later

        logging.debug('selecting muons')
        hzzTightMuonId = (muons.hzzTight>0)
        muons = muons[hzzTightMuonId]

        logging.debug('selecting electrons')
        hzzTightElectronId = (electrons.hzzTight>0)
        electrons = electrons[hzzTightElectronId]

        # TODO: cross clean electrons within DR<0.05 of tight muon

        passFourLeptons = (muons.counts >=4) | (electrons.counts >= 4) | ((muons.counts >= 2) & (electrons.counts >= 2))
        output['cutflow']['four leptons'] += passFourLeptons.sum()
        selection.add('fourLeptons',passFourLeptons)

        
        # build cands
        # remake zz to have same columns
        # pt eta phi mass charge pdgId
        logging.debug('rebuilding leptons')
        def rebuild(leptons):
            return JaggedCandidateArray.candidatesfromoffsets(
                leptons.offsets,
                pt=leptons.pt.flatten(),
                eta=leptons.eta.flatten(),
                phi=leptons.phi.flatten(),
                mass=leptons.mass.flatten(),
                charge=leptons.charge.flatten(),
                pdgId=leptons.pdgId.flatten(),
                # needed for electron SF
                etaSC=leptons.etaSC.flatten() if hasattr(leptons,'etaSC') else leptons.eta.flatten(),
            )
        newMuons = rebuild(muons)
        newElectrons = rebuild(electrons)


        logging.debug('building 2 leptons')
        ee_cands = newElectrons.choose(2)
        mm_cands = newMuons.choose(2)
        
        logging.debug('building 4 leptons')
        zz_4e = newElectrons.choose(4)
        zz_4m = newMuons.choose(4)
        zz_2e2m = ee_cands.cross(mm_cands)
        # for some reason cross creates nested
        zz_2e2m['3'] = zz_2e2m['1']['1']
        zz_2e2m['2'] = zz_2e2m['1']['0']
        zz_2e2m['1'] = zz_2e2m['0']['1']
        zz_2e2m['0'] = zz_2e2m['0']['0']

        # combine them
        zz = JaggedArray.concatenate([zz_4e,zz_4m,zz_2e2m], axis=1)

        # TODO: reevaluate best combination to match HZZ
        # TODO: include FSR
        def massmetric(cands, i, j):
            z1mass = (cands['%d' % i]['p4'] + cands['%d' % j]['p4']).mass
            k, l = set(range(4)) - {i, j}
            z2mass = (cands['%d' % k]['p4'] + cands['%d' % l]['p4']).mass
            deltam = np.abs(z1mass - ZMASS)
            deltaq = np.abs(cands['%d' % i]['charge'] + cands['%d' % j]['charge'])
            deltaf1 = np.abs(cands['%d' % i]['pdgId'] + cands['%d' % j]['pdgId'])
            deltaf2 = np.abs(cands['%d' % k]['pdgId'] + cands['%d' % l]['pdgId'])
            # inflate deltam to absurd number if charge sum is nonzero or different flavor
            return z1mass, z2mass, deltam + 1e10*deltaq + 1e10*deltaf1 + 1e10*deltaf2


        def bestcombination(zzcands):
            good_charge = sum(zzcands[str(i)]['charge'] for i in range(4)) == 0
            good_event = good_charge.sum() == 1
            # this downselection keeps all events where exactly one candidate satisfies the requirement
            # but does not reduce the number of events, i.e. len(zz) stays the same
            # TODO: dont veto on 5th, select best
            zzcands = zzcands[good_charge*good_event][:,:1]
            if zzcands.counts.sum() == 0:
                # empty array (because a bug in concatenate makes it fail on empty arrays)
                empty = JaggedArray.fromcounts(np.zeros(len(zzcands), dtype='i'), [])
                return empty, empty
            # now we have to check the permutations of leptons for closest mass to Z boson
            # only 4 of these 6 permutations are valid charge pairs, but its easier
            # to compare them all, and assign a large delta mass rather than figure out which
            # are valid beforehand
            z1mass = []
            z2mass = []
            iperm = []
            z11 = []
            z12 = []
            z21 = []
            z22 = []
            for i,j in [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]:
                z1, z2, idx = massmetric(zzcands, i, j)
                z1mass.append(z1)
                z2mass.append(z2)
                iperm.append(idx)
                k, l = set(range(4)) - {i,j}
                i = z1.ones_like(dtype=int)*i
                j = z1.ones_like(dtype=int)*j
                k = z1.ones_like(dtype=int)*k
                l = z1.ones_like(dtype=int)*l
                z11.append(i)
                z12.append(j)
                z21.append(k)
                z22.append(l)


            z1mass = JaggedArray.concatenate(z1mass, axis=1)
            z2mass = JaggedArray.concatenate(z2mass, axis=1)
            iperm = JaggedArray.concatenate(iperm, axis=1)
            z11 = JaggedArray.concatenate(z11, axis=1)
            z12 = JaggedArray.concatenate(z12, axis=1)
            z21 = JaggedArray.concatenate(z21, axis=1)
            z22 = JaggedArray.concatenate(z22, axis=1)
            z1mass = z1mass[iperm.argmin()]
            z2mass = z2mass[iperm.argmin()]
            z11 = z11[iperm.argmin()].astype(int).astype(str)
            z12 = z12[iperm.argmin()].astype(int).astype(str)
            z21 = z21[iperm.argmin()].astype(int).astype(str)
            z22 = z22[iperm.argmin()].astype(int).astype(str)
            return z1mass, z2mass, z11, z12, z21, z22


        logging.debug('selecting best combinations')
        z1, z2, z11, z12, z21, z22 = bestcombination(zz)

        passZCand = ((z1.counts>0) & (z2.counts>0))
        output['cutflow']['z cand'] += passZCand.sum()
        selection.add('zCand',passZCand)

        # im sure there is a better way, but for now just do this
        def get_lepton_values(zl,key):
            val = np.zeros_like(zl.flatten(),dtype=float)
            for i in range(4):
                mask = (str(i)==zl.flatten())
                if key=='pt':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].pt
                elif key=='eta':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].eta
                elif key=='phi':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].phi
                elif key=='mass':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].mass
                else:
                    val[mask] = zz[passZCand][str(i)].flatten()[mask][key]
            return JaggedArray.fromoffsets(zl.offsets,val)

        

        chanSels = {}
        z11pdg = get_lepton_values(z11,'pdgId')
        z12pdg = get_lepton_values(z12,'pdgId')
        z21pdg = get_lepton_values(z21,'pdgId')
        z22pdg = get_lepton_values(z22,'pdgId')
        for chan in ['4e','4m','2e2m','2m2e']:
            if chan=='4e':
                pdgIds = (11,11,11,11)
            if chan=='4m':
                pdgIds = (13,13,13,13)
            if chan=='2e2m':
                pdgIds = (11,11,13,13)
            if chan=='2m2e':
                pdgIds = (13,13,11,11)
            chanSels[chan] = ((abs(z11pdg)==pdgIds[0])
                            & (abs(z12pdg)==pdgIds[1])
                            & (abs(z21pdg)==pdgIds[2])
                            & (abs(z22pdg)==pdgIds[3]))

        # TODO lepton scalefactors
        weights = processor.Weights(df.size)
        if self._isData: 
            output['sumw'][dataset] = 0 # always set to 0 for data
        else:
            output['sumw'][dataset] += df['genWeight'].sum()
            weights.add('genWeight',df['genWeight'])
            weights.add('pileupWeight',
                        self._corrections[f'pileupWeight{self._year}'](df['Pileup_nPU']),
                        self._corrections[f'pileupWeight{self._year}Up'](df['Pileup_nPU']),
                        self._corrections[f'pileupWeight{self._year}Down'](df['Pileup_nPU']),
                        )
            # electron sf
            # gap: 1.442 1.566
            zls = [z11, z12, z21, z22]
            for ei, zl in enumerate(zls):
                ei = str(ei)
                eta = get_lepton_values(zl,'etaSC')
                pt = get_lepton_values(zl,'pt')
                electronRecoSF = self._corrections['electron_reco'](eta,pt)
                electronIdSF = self._corrections['electron_hzz_id_nogap'](eta,pt)
                gapMask = ((abs(eta)>1.442) & (abs(eta)<1.566))
                electronIdSF[gapMask] = self._corrections['electron_hzz_id_gap'](eta,pt)[gapMask]
                electronSF = np.ones_like(electronRecoSF)
                if ei in ['0','1']:
                    chans = ['4e','2e2m']
                if ei in ['2','3']:
                    chans = ['4e','2m2e']
                for chan in chans:
                    chanSel = (chanSels[chan].ones_like().sum()>0) # turns empty arrays into 0's, nonempty int 1's
                    electronSF[chanSel] *= electronRecoSF[chanSel].prod()
                    electronSF[chanSel] *= electronIdSF[chanSel].prod()
                weights.add('electronSF'+ei,electronSF)
            # TODO: muon sf

        logging.debug('filling')
        for chan in ['4e','4m','2e2m','2m2e']:
            # TODO: all selections and scalefactors
            cut = selection.all('trigger','goodVertex','fourLeptons','zCand')
            chanSel = chanSels[chan]
            weight = chanSel.astype(float) * weights.weight()

            output['mass'].fill(
                dataset=dataset,
                channel=chan,
                mass=zz[cut].p4.mass.flatten(),
                weight=weight[cut].flatten(),
            )
            output['z1mass'].fill(
                dataset=dataset,
                channel=chan,
                mass=z1[cut].flatten(),
                weight=weight[cut].flatten(),
            )
            output['z2mass'].fill(
                dataset=dataset,
                channel=chan,
                mass=z2[cut].flatten(),
                weight=weight[cut].flatten(),
            )
            #output['pt_lead'].fill(
            #    dataset=dataset,
            #    channel=channel,
            #    pt=pt_lead.flatten(),
            #)

        return output


    def postprocess(self, accumulator):
        lumis = {
            '2016': 35920,
            '2017': 41530,
            '2018': 59740,
        }
        lumi = lumis[self._year]
        scale = {}
        for dataset, sumw in accumulator['sumw'].items():
            if not sumw: continue
            if dataset in self._corrections['xsec']:
                scale[dataset] = lumi*self._corrections['xsec'][dataset]/sumw
            else:
                print(f'missing cross section for {dataset}')
                scale[dataset] = lumi / sumw
            
        for h in accumulator.values():
            if isinstance(h, hist.Hist):
                h.scale(scale, axis="dataset")
 
        return accumulator


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='HZZ processor')
    parser.add_argument('year', choices=['2016', '2017', '2018'], default='2018', help='Data taking year')
    args = parser.parse_args()

    processor_instance = HZZProcessor(
        year=args.year,
    )

    save(processor_instance, f'hzzProcessor_{args.year}.coffea')
