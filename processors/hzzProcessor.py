#!/usr/bin/env python
from __future__ import print_function, division
import sys
import logging
import numpy as np
import itertools
from coffea import hist, processor, lookup_tools
from coffea.lumi_tools import lumi_tools
from coffea.util import load, save
from coffea.analysis_objects import JaggedCandidateArray
# TODO: migrate to awkward1
from awkward import JaggedArray

ZMASS = 91.1876

logger = logging.getLogger("MonoHZZProcessor")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                    format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')


class HZZProcessor(processor.ProcessorABC):
    # will sync with
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4lRunIILegacy
    # to the best of NanoAOD ability
    def __init__(self, year='2018', corrections={}):
        self._year = year

        self._corrections = corrections
        self._rochester = lookup_tools.rochester_lookup.rochester_lookup(corrections['rochester_data'])

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        channel_axis = hist.Cat("channel", "Channel")
        mass_axis = hist.Bin("mass", r"$m_{4\ell}$ [GeV]", 1200, 0.25, 600)
        zmass_axis = hist.Bin("mass", r"$m_{2\ell}$ [GeV]", 240, 0, 120)
        met_axis = hist.Bin("met", r"$E_{T}^{miss}$ [GeV]", 3000, 0, 3000)
        npvs_axis = hist.Bin("npvs", "Number of Vertices", 120, 0, 120)

        self._selections = ['hzz', 'massWindow']

        hist.Hist.DEFAULT_DTYPE = 'f'  # save some space by keeping float bin counts instead of double
        self._accumulator = processor.dict_accumulator()
        for sel in self._selections:
            self._accumulator[sel + '_mass'] = hist.Hist("Counts", dataset_axis, channel_axis, mass_axis)
            self._accumulator[sel + '_z1mass'] = hist.Hist("Counts", dataset_axis, channel_axis, zmass_axis)
            self._accumulator[sel + '_z2mass'] = hist.Hist("Counts", dataset_axis, channel_axis, zmass_axis)
            self._accumulator[sel + '_met'] = hist.Hist("Counts", dataset_axis, channel_axis, met_axis)
            self._accumulator[sel + '_pileup'] = hist.Hist("Counts", dataset_axis, channel_axis, npvs_axis)

        self._accumulator['cutflow'] = processor.defaultdict_accumulator(int)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(int)

    @property
    def accumulator(self):
        return self._accumulator

    def _add_muon_id(self, muons):
        ptCut = (muons.pt > 5)
        etaCut = (abs(muons.eta) < 2.4)
        dxyCut = (abs(muons.dxy) < 0.5)
        dzCut = (abs(muons.dz) < 1)
        # TODO following not in nanoaod: & (muons.muonBestTrackType != 2)
        idCut = ((muons.isGlobal) | ((muons.isTracker) & (muons.nStations > 0)))
        sipCut = (abs(muons.sip3d) < 4)

        hzzLooseNoIso = (ptCut & etaCut & dxyCut & dzCut & idCut & sipCut)

        isoCut = (muons.pfRelIso03_all < 0.35)

        hzzLoose = (hzzLooseNoIso & isoCut)

        tightCut = ((muons.isPFcand) | ((muons.highPtId > 0) & (muons.pt > 200)))

        hzzTight = (hzzLoose & tightCut)

        # TODO: ghost cleaning

        muons['hzzLooseNoIso'] = hzzLooseNoIso
        muons['hzzLoose'] = hzzLoose
        muons['hzzTight'] = hzzTight

    def _add_electron_id(self, electrons):
        ptCut = (electrons.pt > 7)
        etaCut = (abs(electrons.eta) < 2.5)
        dxyCut = (abs(electrons.dxy) < 0.5)
        dzCut = (abs(electrons.dz) < 1)
        sipCut = (abs(electrons.sip3d) < 4)

        hzzLooseNoIso = (ptCut & etaCut & dxyCut & dzCut & sipCut)

        isoCut = (electrons.pfRelIso03_all < 0.35)

        hzzLoose = (hzzLooseNoIso & isoCut)

        # TODO: check
        pt = electrons.pt
        eta = abs(electrons.eta + electrons.deltaEtaSC)
        mva = electrons.mvaFall17V2Iso
        # these are the 2017 trainings only, because that is what is available
        tightBin00 = ((pt <= 10) & (eta < 0.8) & (mva > 0.8955937602))
        tightBin01 = ((pt <= 10) & ((eta >= 0.8) & (eta < 1.479)) & (mva > 0.91106464032))
        tightBin02 = ((pt <= 10) & (eta >= 1.479) & (mva > 0.94067753025))
        tightBin10 = ((pt > 10) & (eta < 0.8) & (mva > 0.04240620843))
        tightBin11 = ((pt > 10) & ((eta >= 0.8) & (eta < 1.479)) & (mva > 0.0047338429))
        tightBin12 = ((pt > 10) & (eta >= 1.479) & (mva > -0.60423293572))
        tightCut = (tightBin00 | tightBin01 | tightBin02 | tightBin10 | tightBin11 | tightBin12)

        hzzTight = (hzzLoose & tightCut)

        electrons['hzzLooseNoIso'] = hzzLooseNoIso
        electrons['hzzLoose'] = hzzLoose
        electrons['hzzTight'] = hzzTight

    def _add_photon_id(self, photons, muons, electrons):
        ptCut = (photons.pt > 2)
        etaCut = (abs(photons.eta) < 2.4)
        isoCut = (photons.relIso03 < 1.8)

        hzzPreselection = (ptCut & etaCut & isoCut)

        # matching muons is all that is available now
        drOverEtCut = (photons.dROverEt2 < 0.012)
        drMuonCut = (photons.delta_r(photons.matched_muon) < 0.5)

        hzzMatchMuon = (hzzPreselection & drOverEtCut & drMuonCut)
        hzzMatchElectron = (hzzPreselection & (np.zeros_like(hzzPreselection)))

        photons['hzzPreselection'] = hzzPreselection
        photons['hzzMatchMuon'] = hzzMatchMuon
        photons['hzzMatchElectron'] = hzzMatchElectron

    def _add_trigger(self, events):
        dataset = events.metadata['dataset']

        triggerPaths = {}
        # DoubleMuon
        if self._year == '2016':
            triggerPaths['DoubleMuon'] = [
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                "Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                "TripleMu_12_10_5",
            ]
        elif self._year == '2017':
            triggerPaths['DoubleMuon'] = [
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                "TripleMu_10_5_5_DZ",
                "TripleMu_12_10_5",
            ]
        elif self._year == '2018':
            triggerPaths['DoubleMuon'] = [
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "TripleMu_10_5_5_DZ",
                "TripleMu_12_10_5",
            ]

        # DoubleEG
        if self._year == '2016':
            triggerPaths['DoubleEG'] = [
                "Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "DoubleEle33_CaloIdL_GsfTrkIdVL_MW",
                "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ]
        elif self._year == '2017':
            triggerPaths['DoubleEG'] = [
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "DoubleEle33_CaloIdL_MW",
                "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ]

        # MuonEG
        if self._year == '2016':
            triggerPaths['MuonEG'] = [
                "Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                "Mu8_DiEle12_CaloIdL_TrackIdL",
                "DiMu9_Ele9_CaloIdL_TrackIdL",
            ]
        elif self._year == '2017':
            triggerPaths['MuonEG'] = [
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
                "Mu8_DiEle12_CaloIdL_TrackIdL",
                "Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
            ]
        elif self._year == '2018':
            triggerPaths['MuonEG'] = [
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
                "Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
            ]

        # EGamma
        if self._year == '2018':
            triggerPaths['EGamma'] = [
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "DoubleEle25_CaloIdL_MW",
                "Ele32_WPTight_Gsf",
            ]

        # SingleMuon
        if self._year == '2016':
            triggerPaths['SingleMuon'] = [
                "IsoMu20",
                "IsoTkMu20",
                "IsoMu22",
                "IsoTkMu22",
                "IsoMu24",
                "IsoTkMu24",
            ]
        elif self._year == '2017':
            triggerPaths['SingleMuon'] = [
                "IsoMu27",
            ]
        elif self._year == '2018':
            triggerPaths['SingleMuon'] = [
                "IsoMu24",
            ]

        # SingleElectron
        if self._year == '2016':
            triggerPaths['SingleElectron'] = [
                "Ele25_eta2p1_WPTight_Gsf",
                "Ele27_WPTight_Gsf",
                "Ele27_eta2p1_WPLoose_Gsf",
                "Ele32_eta2p1_WPTight_Gsf",
            ]
        elif self._year == '2017':
            triggerPaths['SingleElectron'] = [
                "Ele35_WPTight_Gsf",
                "Ele38_WPTight_Gsf",
                "Ele40_WPTight_Gsf",
            ]

        # Define priority
        # To avoid double counting in data, for a given dataset
        # all current datasets triggers are accepted
        # and all higher datasets triggers are vetoed
        # no lower datasets triggers are looked at
        # in MC, all triggers are accepted
        if self._year == '2016' or self._year == '2017':
            triggerPriority = [
                'DoubleEG',
                'DoubleMuon',
                'MuonEG',
                'SingleElectron',
                'SingleMuon',
            ]
        else:
            triggerPriority = [
                'EGamma',
                'DoubleMuon',
                'MuonEG',
                'SingleMuon',
            ]

        triggersToAccept = []
        triggersToVeto = []
        accept = not self._isData
        for d in triggerPriority:
            if d == dataset:
                # start accepting triggers
                accept = True
            for p in triggerPaths[d]:
                if accept:
                    triggersToAccept += [p]
                else:
                    triggersToVeto += [p]
            if d == dataset:
                # don't need to look at rest of trigger paths in data
                break

        # TODO: no guarantee the trigger is in every dataset?
        # for now, check, but should find out
        result = np.zeros_like(events.run, dtype=bool)
        for p in triggersToAccept:
            if p not in events.HLT.columns:
                continue
            result = ((result) | (events.HLT[p]))
        for p in triggersToVeto:
            if p not in events.HLT.columns:
                continue
            result = ((result) & (~events.HLT[p]))

        events['passHLT'] = result

    def process(self, events):
        logging.debug('starting process')
        output = self.accumulator.identity()

        dataset = events.metadata['dataset']
        self._isData = dataset in ['SingleMuon', 'DoubleMuon', 'SingleElectron', 'DoubleEG', 'EGamma', 'MuonEG']

        selection = processor.PackedSelection()

        output['cutflow']['all events'] += events.size

        logging.debug('applying lumi mask')
        if self._isData:
            lumiMask = lumi_tools.LumiMask(self._corrections['golden'])
            events['passLumiMask'] = lumiMask(np.array(events.run), np.array(events.luminosityBlock))
        else:
            events['passLumiMask'] = np.ones_like(events.run, dtype=bool)
        passLumiMask = events.passLumiMask
        selection.add('lumiMask', passLumiMask)

        logging.debug('adding trigger')
        self._add_trigger(events)

        passHLT = events.passHLT
        selection.add('trigger', passHLT)
        output['cutflow']['pass trigger'] += passHLT.sum()
        # if no trigger: fast return
        if passHLT.sum() == 0:
            return output

        # require one good vertex
        logging.debug('checking vertices')
        passGoodVertex = (events.PV.npvsGood > 0)
        output['cutflow']['good vertex'] += passGoodVertex.sum()
        selection.add('goodVertex', passGoodVertex)

        # run rochester
        rochester = self._rochester
        _muon_offsets = events.Muon.pt.offsets
        _charge = events.Muon.charge
        _pt = events.Muon.pt
        _eta = events.Muon.eta
        _phi = events.Muon.phi
        if self._isData:
            _k = rochester.kScaleDT(_charge, _pt, _eta, _phi)
            # _kErr = rochester.kScaleDTerror(_charge, _pt, _eta, _phi)
        else:
            # for default if gen present
            _gpt = events.Muon.matched_gen.pt
            # for backup w/o gen
            _nl = events.Muon.nTrackerLayers
            _u = JaggedArray.fromoffsets(_muon_offsets, np.random.rand(*_pt.flatten().shape))
            _hasgen = (_gpt.fillna(-1) > 0)
            _kspread = rochester.kSpreadMC(_charge[_hasgen], _pt[_hasgen], _eta[_hasgen], _phi[_hasgen],
                                           _gpt[_hasgen])
            _ksmear = rochester.kSmearMC(_charge[~_hasgen], _pt[~_hasgen], _eta[~_hasgen], _phi[~_hasgen],
                                         _nl[~_hasgen], _u[~_hasgen])
            _k = np.ones_like(_pt.flatten())
            _k[_hasgen.flatten()] = _kspread.flatten()
            _k[~_hasgen.flatten()] = _ksmear.flatten()
            _k = JaggedArray.fromoffsets(_muon_offsets, _k)
            # _kErrspread = rochester.kSpreadMCerror(_charge[_hasgen], _pt[_hasgen], _eta[_hasgen], _phi[_hasgen],
            #                                        _gpt[_hasgen])
            # _kErrsmear  = rochester.kSmearMCerror(_charge[~_hasgen], _pt[~_hasgen], _eta[~_hasgen], _phi[~_hasgen],
            #                                       _nl[~_hasgen], _u[~_hasgen])
            # _kErr = np.ones_like(_pt.flatten())
            # _kErr[_hasgen.flatten()] = _kErrspread.flatten()
            # _kErr[~_hasgen.flatten()] = _kErrsmear.flatten()
            # _kErr = JaggedArray.fromoffsets(_muon_offsets, _kErr)

        mask = _pt.flatten() < 200
        rochester_pt = _pt.flatten()
        rochester_pt[mask] = (_k * _pt).flatten()[mask]
        events.Muon['pt'] = JaggedArray.fromoffsets(_muon_offsets, rochester_pt)

        logging.debug('adding muon id')
        self._add_muon_id(events.Muon)
        logging.debug('adding electron id')
        self._add_electron_id(events.Electron)
        logging.debug('add fsr photon id')
        self._add_photon_id(events.FsrPhoton, events.Muon, events.Electron)

        # TODO, select loose cands here, tight cands later

        logging.debug('selecting muons')
        hzzTightMuonId = (events.Muon.hzzTight > 0)
        muons = events.Muon[hzzTightMuonId]

        logging.debug('selecting electrons')
        hzzTightElectronId = (events.Electron.hzzTight > 0)
        electrons = events.Electron[hzzTightElectronId]

        # TODO: cross clean electrons within DR<0.05 of tight muon

        passFourLeptons = ((muons.counts >= 4) | (electrons.counts >= 4)
                           | ((muons.counts >= 2) & (electrons.counts >= 2)))
        output['cutflow']['four leptons'] += passFourLeptons.sum()
        selection.add('fourLeptons', passFourLeptons)

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
                etaSC=(leptons.eta + leptons.deltaEtaSC).flatten() if hasattr(leptons, 'deltaEtaSC') else leptons.eta.flatten(),
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
        zz = JaggedArray.concatenate([zz_4e, zz_4m, zz_2e2m], axis=1)

        def best_zz(zz):
            # zz is a jagged array of all possible combinations of 4 leptons (2e2m, 4e, 4m)
            # we now need to decide which combination and iteration are the best
            # first, exclude combinations that fail charge
            zz = zz[(sum(zz[str(i)]['charge'] for i in range(4)) == 0)]
            # next, the combinations
            for i, j in itertools.combinations(range(4), 2):
                zzi = zz[str(i)]
                zzj = zz[str(j)]
                # all combinations of DR > 0.02
                mask = (zzi.p4.delta_r(zzj.p4) > 0.02)
                # all OS combinations m(ll) > 4 not including fsr, regardless of flavor
                opsign = (zzi['charge']+zzj['charge'] == 0)
                mij = (zzi['p4']+zzj['p4']).mass
                mask = mask & ((opsign & (mij > 4)) | ~opsign)
                zz = zz[mask]
            # pt1 > 20, pt2 > 10
            zz = zz[sum(zz[str(i)]['p4'].pt > 20 for i in range(4)) >= 1]
            zz = zz[sum(zz[str(i)]['p4'].pt > 10 for i in range(4)) >= 2]
            # m(4l)>70 including fsr TODO
            zz = zz[(zz['p4'].mass > 70)]
            # Now check the combinatorics
            # z1 = ij, z2 = kl not all will be valid charge combinations, so check them
            # i/j perms are the best combination sorting
            # i is abs(z1-ZMASS)
            # j is sum z2 lepton pt
            # TODO: HZZ actually uses a P_sig/P_bkg to choose best instead of sum z2 lepton pt
            iperm = []
            jperm = []

            z1mass = []
            z2mass = []
            z11 = []
            z12 = []
            z21 = []
            z22 = []
            zzindex = []
            for i, j in itertools.combinations(range(4), 2):
                k, l = set(range(4)) - {i, j}
                zzi = zz[str(i)]
                zzj = zz[str(j)]
                zzk = zz[str(k)]
                zzl = zz[str(l)]
                z1 = (zzi['p4'] + zzj['p4']).mass
                z2 = (zzk['p4'] + zzl['p4']).mass
                deltam = np.abs(z1 - ZMASS)
                idx = deltam
                jdx = -1*(zzk['p4'].pt + zzl['p4'].pt)
                zi = z1.ones_like(dtype=int)*i
                zj = z1.ones_like(dtype=int)*j
                zk = z1.ones_like(dtype=int)*k
                zl = z1.ones_like(dtype=int)*l
                four_flav = (abs(zzi['pdgId']) == abs(zzk['pdgId']))
                ik_os = (zzi['charge'] + zzk['charge']) == 0
                za_ik = (zzi['p4'] + zzk['p4']).mass
                za_il = (zzi['p4'] + zzl['p4']).mass
                zb_jk = (zzj['p4'] + zzk['p4']).mass
                zb_jl = (zzj['p4'] + zzl['p4']).mass
                za = za_ik.flatten()
                zb = zb_jl.flatten()
                za[~ik_os.flatten()] = za_il[~ik_os].flatten()
                zb[~ik_os.flatten()] = zb_jk[~ik_os].flatten()
                zap = np.where(za > zb, za, zb)
                zbp = np.where(za > zb, zb, za)
                za = JaggedArray.fromoffsets(ik_os.offsets, zap)
                zb = JaggedArray.fromoffsets(ik_os.offsets, zbp)

                # verify charge and flavor
                permmask = (zzi['charge'] + zzj['charge']) == 0
                permmask = permmask & (zzk['charge'] + zzl['charge']) == 0
                permmask = permmask & (zzi['pdgId'] + zzj['pdgId']) == 0
                permmask = permmask & (zzk['pdgId'] + zzl['pdgId']) == 0
                # z1 is closest to nominal z mass including fsr TODO
                permmask = permmask & (np.abs(z1 - ZMASS) < np.abs(z2 - ZMASS))
                # 40 < z1 < 120, 12 < z2 < 120 including fsr TODO
                permmask = permmask & ((z1 > 40) & (z1 < 120) & (z2 > 12) & (z2 < 120))
                # smart cut: za and zb are mass sorted alternative pairings (4e/4m).
                # require !(abs(za-mZ) < abs(z1-mZ) && zb<12) including fsr TODO
                permmask = permmask & ((four_flav & ~((abs(za - ZMASS) < abs(z1 - ZMASS)) & (zb < 12))) | ~four_flav)

                iperm.append(idx[permmask])
                jperm.append(jdx[permmask])

                z1mass.append(z1[permmask])
                z2mass.append(z2[permmask])
                z11.append(zi[permmask])
                z12.append(zj[permmask])
                z21.append(zk[permmask])
                z22.append(zl[permmask])
                zzindex.append(zz.localindex[permmask])

            # a bit hacky, but make the second way smaller in scale
            # number chosen such they 14 TeV is smaller than ~10 MeV
            iperm = JaggedArray.concatenate(iperm, axis=1)
            jperm = JaggedArray.concatenate(jperm, axis=1)
            ijperm = iperm+1e-6*jperm

            z1mass = JaggedArray.concatenate(z1mass, axis=1)
            z2mass = JaggedArray.concatenate(z2mass, axis=1)
            z11 = JaggedArray.concatenate(z11, axis=1)
            z12 = JaggedArray.concatenate(z12, axis=1)
            z21 = JaggedArray.concatenate(z21, axis=1)
            z22 = JaggedArray.concatenate(z22, axis=1)
            zzindex = JaggedArray.concatenate(zzindex, axis=1)

            z1mass = z1mass[ijperm.argmin()]
            z2mass = z2mass[ijperm.argmin()]
            z11 = z11[ijperm.argmin()].astype(int).astype(str)
            z12 = z12[ijperm.argmin()].astype(int).astype(str)
            z21 = z21[ijperm.argmin()].astype(int).astype(str)
            z22 = z22[ijperm.argmin()].astype(int).astype(str)
            zzindex = zzindex[ijperm.argmin()]

            return zz[zzindex], z1mass, z2mass, z11, z12, z21, z22

        logging.debug('selecting best combinations')
        zz, z1, z2, z11, z12, z21, z22 = best_zz(zz)

        passZCand = ((z1.counts > 0) & (z2.counts > 0) & ((z1 > 40) & (z1 < 120) & (z2 > 12) & (z2 < 120)).counts > 0)
        output['cutflow']['z cand'] += passZCand.sum()
        selection.add('zCand', passZCand)

        passMassWindow = (passZCand & zz[((zz.p4.mass > 115) & (zz.p4.mass < 135))].counts > 0)
        output['cutflow']['mass window'] += passMassWindow.sum()
        selection.add('massWindow', passMassWindow)

        # im sure there is a better way, but for now just do this
        def get_lepton_values(zl, key):
            val = np.zeros_like(zl.flatten(), dtype=float)
            if len(val) == 0:
                return JaggedArray.fromoffsets(zl.offsets, val)
            for i in range(4):
                mask = (str(i) == zl.flatten())
                if key == 'pt':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].pt
                elif key == 'eta':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].eta
                elif key == 'phi':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].phi
                elif key == 'mass':
                    val[mask] = zz[passZCand][str(i)].flatten()[mask]['p4'].mass
                else:
                    val[mask] = zz[passZCand][str(i)].flatten()[mask][key]
            return JaggedArray.fromoffsets(zl.offsets, val)

        chanSels = {}
        z11pdg = get_lepton_values(z11, 'pdgId')
        z12pdg = get_lepton_values(z12, 'pdgId')
        z21pdg = get_lepton_values(z21, 'pdgId')
        z22pdg = get_lepton_values(z22, 'pdgId')
        for chan in ['4e', '4m', '2e2m', '2m2e']:
            if chan == '4e':
                pdgIds = (11, 11, 11, 11)
            if chan == '4m':
                pdgIds = (13, 13, 13, 13)
            if chan == '2e2m':
                pdgIds = (11, 11, 13, 13)
            if chan == '2m2e':
                pdgIds = (13, 13, 11, 11)
            chanSels[chan] = ((abs(z11pdg) == pdgIds[0])
                              & (abs(z12pdg) == pdgIds[1])
                              & (abs(z21pdg) == pdgIds[2])
                              & (abs(z22pdg) == pdgIds[3]))

        # TODO lepton scalefactors
        weights = processor.Weights(events.run.size)
        if self._isData:
            output['sumw'][dataset] = 0  # always set to 0 for data
        else:
            output['sumw'][dataset] += events.genWeight.sum()
            weights.add('genWeight', events.genWeight)
            weights.add('pileupWeight',
                        self._corrections['pileupWeight'](events.Pileup.nPU),
                        self._corrections['pileupWeightUp'](events.Pileup.nPU),
                        self._corrections['pileupWeightDown'](events.Pileup.nPU),
                        )
            # electron sf
            # gap: 1.442 1.566
            zls = [z11, z12, z21, z22]
            for ei, zl in enumerate(zls):
                ei = str(ei)
                eta = get_lepton_values(zl, 'etaSC')
                pt = get_lepton_values(zl, 'pt')
                electronRecoSF = self._corrections['electron_reco'](eta, pt)
                electronIdSF = self._corrections['electron_hzz_id_nogap'](eta, pt)
                gapMask = ((abs(eta) > 1.442) & (abs(eta) < 1.566))
                electronIdSF[gapMask] = self._corrections['electron_hzz_id_gap'](eta, pt)[gapMask]
                electronSF = np.ones_like(electronRecoSF.prod())
                if ei in ['0', '1']:
                    chans = ['4e', '2e2m']
                if ei in ['2', '3']:
                    chans = ['4e', '2m2e']
                for chan in chans:
                    # turns empty arrays into 0's, nonempty int 1's
                    chanSel = (chanSels[chan].ones_like().sum() > 0)
                    electronSF[chanSel] *= electronRecoSF[chanSel].prod()
                    electronSF[chanSel] *= electronIdSF[chanSel].prod()
                weights.add('electronSF'+ei, electronSF)
            # TODO: muon sf

        logging.debug('filling')
        for sel in self._selections:
            # TODO: all selections and scalefactors
            if sel == 'hzz':
                cut = selection.all('lumiMask', 'trigger', 'goodVertex', 'fourLeptons',
                                    'zCand')
            elif sel == 'massWindow':
                cut = selection.all('lumiMask', 'trigger', 'goodVertex', 'fourLeptons',
                                    'zCand', 'massWindow')
            for chan in ['4e', '4m', '2e2m', '2m2e']:
                chanSel = chanSels[chan]
                weight = chanSel.astype(float) * weights.weight()

                output[sel+'_mass'].fill(
                    dataset=dataset,
                    channel=chan,
                    mass=zz[cut].p4.mass.flatten(),
                    weight=weight[cut].flatten(),
                )
                output[sel+'_z1mass'].fill(
                    dataset=dataset,
                    channel=chan,
                    mass=z1[cut].flatten(),
                    weight=weight[cut].flatten(),
                )
                output[sel+'_z2mass'].fill(
                    dataset=dataset,
                    channel=chan,
                    mass=z2[cut].flatten(),
                    weight=weight[cut].flatten(),
                )
                output[sel+'_met'].fill(
                    dataset=dataset,
                    channel=chan,
                    met=events.MET.pt[cut],
                    weight=weight[cut].flatten(),
                )
                output[sel+'_pileup'].fill(
                    dataset=dataset,
                    channel=chan,
                    npvs=events.PV.npvs[cut],
                    weight=weight[cut].flatten(),
                )

        return output

    def postprocess(self, accumulator):
        # always scale to 1000 pb for plotting
        lumi = 1000
        scale = {}
        for dataset, sumw in accumulator['sumw'].items():
            if not sumw:
                continue
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

    years = ['2016', '2017', '2018']
    for year in years:
        corrections = load(f'corrections/corrections_{year}.coffea')

        processor_instance = HZZProcessor(
            year=year,
            corrections=corrections,
        )

        save(processor_instance, f'processors/hzzProcessor_{year}.coffea')
