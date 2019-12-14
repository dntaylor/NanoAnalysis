#!/usr/bin/env python
from __future__ import print_function, division
import os
import sys
import logging
import argparse
import numpy as np
import uproot
import uproot_methods
from coffea import hist, processor
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

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        channel_axis = hist.Cat("channel", "Channel")
        mass_axis = hist.Bin("mass", r"$m_{4\ell}$ [GeV]", 600, 0.25, 300)
        zmass_axis = hist.Bin("mass", r"$m_{2\ell}$ [GeV]", 240, 0, 120)
        pt_axis = hist.Bin("pt", r"$p_{T,\ell}$ [GeV]", 3000, 0.25, 300)
        met_axis = hist.Bin("pt", r"$E_{T}^{miss}$ [GeV]", 3000, 0, 3000)

        self._accumulator = processor.dict_accumulator({
            'mass': hist.Hist("Counts", dataset_axis, channel_axis, mass_axis),
            'z1mass': hist.Hist("Counts", dataset_axis, channel_axis, zmass_axis),
            'z2mass': hist.Hist("Counts", dataset_axis, channel_axis, zmass_axis),
            #'met': hist.Hist("Counts", dataset_axis, channel_axis, met_axis),
            #'pt_lead': hist.Hist("Counts", dataset_axis, channel_axis, pt_axis),
            'cutflow': processor.defaultdict_accumulator(int),
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
        eta = abs(electrons.eta+electrons.deltaEtaSC)
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

        result = np.zeros_like(df['event'])
        for p in triggersToAccept:
            result = ((result) | (df[p]))
        for p in triggersToVeto:
            result = ((result) & (~df[p]))

        df['passHLT'] = result
                


    def process(self, df):
        logging.debug('starting process')
        output = self.accumulator.identity()

        # TODO: instead of cutflow, use processor.PackedSelection
        output['cutflow']['all events'] += df['event'].size

        logging.debug('adding trigger')
        self._add_trigger(df)

        output['cutflow']['pass trigger'] += df['passHLT'].sum()

        
        # require one good vertex
        logging.debug('checking vertices')
        output['cutflow']['good vertex'] += (df['PV_npvsGood']>0).sum()


        dataset = df['dataset']
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
            mva=df['Electron_{}'.format(mvaDisc)],
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

        fourleptons = (muons.counts >=4) | (electrons.counts >= 4) | ((muons.counts >= 2) & (electrons.counts >= 2))
        output['cutflow']['four leptons'] += fourleptons.sum()
        
        # build cands
        logging.debug('building 2 leptons')
        ee_cands = electrons.choose(2)
        mm_cands = muons.choose(2)
        
        logging.debug('building 4 leptons')
        zz_4e = electrons.choose(4)
        zz_4m = muons.choose(4)
        zz_2e2m = ee_cands.cross(mm_cands)
        # for some reason cross creates nested
        zz_2e2m['3'] = zz_2e2m['1']['1']
        zz_2e2m['2'] = zz_2e2m['1']['0']
        zz_2e2m['1'] = zz_2e2m['0']['1']
        zz_2e2m['0'] = zz_2e2m['0']['0']
        
        # TODO: reevaluate best combination to match HZZ
        # and include FSR
        def massmetric(cands, i, j):
            z1mass = (cands['%d' % i]['p4'] + cands['%d' % j]['p4']).mass
            k, l = set(range(4)) - {i, j}
            z2mass = (cands['%d' % k]['p4'] + cands['%d' % l]['p4']).mass
            deltam = np.abs(z1mass - ZMASS)
            deltaq = np.abs(cands['%d' % i]['charge'] + cands['%d' % j]['charge'])
            # inflate deltam to absurd number if charge sum is nonzero
            return z1mass, z2mass, deltam + 1e10*deltaq


        def bestcombination(zzcands):
            good_charge = sum(zzcands[str(i)]['charge'] for i in range(4)) == 0
            good_event = good_charge.sum() == 1
            # this downselection keeps all events where exactly one candidate satisfies the requirement
            # but does not reduce the number of events, i.e. len(zz_4m) stays the same
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
            for i,j in [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]:
                z1, z2, idx = massmetric(zzcands, i, j)
                z1mass.append(z1)
                z2mass.append(z2)
                iperm.append(idx)

            z1mass = JaggedArray.concatenate(z1mass, axis=1)
            z2mass = JaggedArray.concatenate(z2mass, axis=1)
            iperm = JaggedArray.concatenate(iperm, axis=1)
            z1mass = z1mass[iperm.argmin()]
            z2mass = z2mass[iperm.argmin()]
            return z1mass, z2mass


        logging.debug('selecting best combinations')
        z1_4m, z2_4m = bestcombination(zz_4m)
        z1_4e, z2_4e = bestcombination(zz_4e)
        
        # for 2e2m its a bit simpler
        good_charge = (zz_2e2m.i0['charge'] + zz_2e2m.i1['charge'] == 0) & (zz_2e2m.i2['charge'] + zz_2e2m.i3['charge'] == 0)
        good_event = good_charge.sum() == 1
        zz_2e2m = zz_2e2m[good_event*good_charge][:,:1]
        za_2e2m, zb_2e2m, deltam_a = massmetric(zz_2e2m, 0, 1)
        _, _, deltam_b = massmetric(zz_2e2m, 2, 3)
        # this is a good place for awkward.where, but its not available yet
        z_2e2m = JaggedArray.concatenate([za_2e2m, zb_2e2m], axis=1)
        deltam = JaggedArray.concatenate([deltam_a, deltam_b], axis=1)
        z1_2e2m = z_2e2m[deltam.argmin()]
        z2_2e2m = z_2e2m[deltam.argmax()]

        logging.debug('filling')
        output['mass'].fill(
            dataset=dataset,
            channel='4e',
            mass=zz_4e.p4.mass.flatten(),
        )
        output['z1mass'].fill(
            dataset=dataset,
            channel='4e',
            mass=z1_4e.flatten(),
        )
        output['z2mass'].fill(
            dataset=dataset,
            channel='4e',
            mass=z2_4e.flatten(),
        )
        output['mass'].fill(
            dataset=dataset,
            channel='4m',
            mass=zz_4m.p4.mass.flatten(),
        )
        output['z1mass'].fill(
            dataset=dataset,
            channel='4m',
            mass=z1_4m.flatten(),
        )
        output['z2mass'].fill(
            dataset=dataset,
            channel='4m',
            mass=z2_4m.flatten(),
        )
        output['mass'].fill(
            dataset=dataset,
            channel='2e2m',
            mass=zz_2e2m.p4.mass.flatten(),
        )
        output['z1mass'].fill(
            dataset=dataset,
            channel='2e2m',
            mass=z1_2e2m.flatten(),
        )
        output['z2mass'].fill(
            dataset=dataset,
            channel='2e2m',
            mass=z2_2e2m.flatten(),
        )

        #output['pt_lead'].fill(
        #    dataset=dataset,
        #    channel=channel,
        #    pt=pt_lead.flatten(),
        #)
        return output


    def postprocess(self, accumulator):
        return accumulator


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='HZZ processor')
    parser.add_argument('year', choices=['2016', '2017', '2018'], default='2018', help='Data taking year')
    args = parser.parse_args()

    processor_instance = HZZProcessor(
        year=args.year,
    )

    save(processor_instance, 'hzzProcessor_{year}.coffea'.format(year=args.year))
