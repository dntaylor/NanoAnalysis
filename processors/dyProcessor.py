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
from coffea.lumi_tools import lumi_tools
from coffea.util import load, save
from coffea.analysis_objects import JaggedCandidateArray
from awkward import JaggedArray, IndexedArray

ZMASS = 91.1876

logger = logging.getLogger("MonoDYProcessor")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class DYProcessor(processor.ProcessorABC):

    def __init__(self,year='2018',corrections={}):
        self._year = year

        self._corrections = corrections

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        channel_axis = hist.Cat("channel", "Channel")
        zmass_axis = hist.Bin("mass", r"$m_{2\ell}$ [GeV]", 240, 0, 120)
        pt_axis = hist.Bin("pt", r"$p_{T,\ell}$ [GeV]", 3000, 0.25, 300)
        met_axis = hist.Bin("met", r"$E_{T}^{miss}$ [GeV]", 3000, 0, 3000)

        self._selections = ['massWindow']

        hist.Hist.DEFAULT_DTYPE = 'f'  # save some space by keeping float bin counts instead of double
        self._accumulator = processor.dict_accumulator()
        for sel in self._selections:
            self._accumulator[sel + '_zmass'] = hist.Hist("Counts", dataset_axis, channel_axis, zmass_axis)
            self._accumulator[sel + '_met'] = hist.Hist("Counts", dataset_axis, channel_axis, met_axis)

        self._accumulator['cutflow'] = processor.defaultdict_accumulator(int)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(int)
        

    @property
    def accumulator(self):
        return self._accumulator

    def _add_muon_id(self, muons):
        # note: input muons must pass
        # slimmedMuons and (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))
        ptCut = (muons.pt>20)
        etaCut = (abs(muons.eta)<2.4)
        dxyCut = (abs(muons.dxy)<0.5)
        dzCut = (abs(muons.dz)<1)
        idCut = (muons.mediumId)

        idNoIsoCut = (ptCut & etaCut & dxyCut & dzCut & idCut)

        isoCut = (muons.pfRelIso04_all<0.15)

        idIsoCut = (idNoIsoCut & isoCut)

        muons['passId'] = idIsoCut
        
    def _add_electron_id(self, electrons):
        # note: input electrons must pass
        # slimmedElectrons and (pt > 5)
        ptCut = (electrons.pt>20)
        etaCut = (abs(electrons.eta)<2.5)
        dxyCut = (abs(electrons.dxy)<0.5)
        dzCut = (abs(electrons.dz)<1)
        idCut = (electrons.mvaFall17V2Iso_WP90)

        loose = (ptCut & etaCut & dxyCut & dzCut & idCut)

        electrons['passId'] = loose
        
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
            ]
        elif self._year=='2017':
            triggerPaths['DoubleMuon'] = [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            ]
        elif self._year=='2018':
            triggerPaths['DoubleMuon'] = [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            ]

        # DoubleEG
        if self._year=='2016':
            triggerPaths['DoubleEG'] = [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            ]
        elif self._year=='2017':
            triggerPaths['DoubleEG'] = [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            ]

        # EGamma
        if self._year=='2018':
            triggerPaths['EGamma'] = [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Ele32_WPTight_Gsf",
            ]


        # SingleMuon
        if self._year=='2016':
            triggerPaths['SingleMuon'] = [
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
                "HLT_Ele27_WPTight_Gsf",
            ]
        elif self._year=='2017':
            triggerPaths['SingleElectron'] = [
                "HLT_Ele35_WPTight_Gsf",
            ]


        # Define priority
        # To avoid double counting in data, for a given dataset
        # all current datasets triggers are accepted
        # and all higher datasets triggers are vetoed
        # no lower datasets triggers are looked at
        # in MC, all triggers are accepted
        if self._year=='2016' or self._year=='2017':
            triggerPriority = [
                'DoubleMuon',
                'DoubleEG',
                'SingleMuon',
                'SingleElectron',
            ]
        else:
            triggerPriority = [
                'DoubleMuon',
                'EGamma',
                'SingleMuon',
            ]

        triggersToAccept = []
        triggersToVeto = []
        accept = not self._isData
        for d in triggerPriority:
            if d==dataset: accept = True # start accepting triggers
            for p in triggerPaths[d]:
                if accept:
                    triggersToAccept += [p]
                else:
                    triggersToVeto += [p]
            if d==dataset: break # don't need to look at rest of trigger paths in data

        # TODO: no guarantee the trigger is in every dataset?
        # for now, check, but should find out
        result = np.zeros_like(df['event'],dtype=bool)
        for p in triggersToAccept:
            if p not in df: continue
            result = ((result) | (df[p]))
        for p in triggersToVeto:
            if p not in df: continue
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

        logging.debug('applying lumi mask')
        if self._isData:
            lumiMask = lumi_tools.LumiMask(self._corrections['golden'])
            df['passLumiMask'] = lumiMask(df['run'],df['luminosityBlock'])
        else:
            df['passLumiMask'] = np.ones_like(df['run'],dtype=bool)
        passLumiMask = df['passLumiMask']
        selection.add('lumiMask',passLumiMask)
            

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
            mediumId=df['Muon_mediumId'],
            pfRelIso04_all=df['Muon_pfRelIso04_all'],
            pdgId=df['Muon_pdgId'],
        )

        logging.debug('building electrons')
        electrons = JaggedCandidateArray.candidatesfromcounts(
            df['nElectron'],
            pt=df['Electron_pt'],
            eta=df['Electron_eta'],
            phi=df['Electron_phi'],
            mass=df['Electron_mass'],
            charge=df['Electron_charge'],
            dxy=df['Electron_dxy'],
            dz=df['Electron_dz'],
            deltaEtaSC=df['Electron_deltaEtaSC'],
            etaSC=df['Electron_eta']+df['Electron_deltaEtaSC'],
            mvaFall17V2Iso_WP90=df['Electron_mvaFall17V2Iso_WP90'],
            pdgId=df['Electron_pdgId'],
        )

        logging.debug('adding muon id')
        self._add_muon_id(muons)
        logging.debug('adding electron id')
        self._add_electron_id(electrons)

        logging.debug('selecting muons')
        muonId = (muons.passId>0)
        muons = muons[muonId]

        logging.debug('selecting electrons')
        electronId = (electrons.passId>0)
        electrons = electrons[electronId]

        passTwoLeptons = (muons.counts >= 2) | (electrons.counts >= 2)
        output['cutflow']['two leptons'] += passTwoLeptons.sum()
        selection.add('twoLeptons',passTwoLeptons)

        
        # build cands
        # remake z to have same columns
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
        
        # combine them
        z_cands = JaggedArray.concatenate([ee_cands,mm_cands], axis=1)

        def bestcombination(zcands):
            good_charge = sum(zcands[str(i)]['charge'] for i in range(2)) == 0
            # this keeps the first z cand in each event
            # should instead sort the best first
            # TODO: select best
            zcands = zcands[good_charge][:,:1]
            return zcands


        logging.debug('selecting best combinations')
        z_cands = bestcombination(z_cands)

        z1 = np.zeros_like(z_cands['p4'].pt.flatten(), dtype='i')
        z2 = np.ones_like(z_cands['p4'].pt.flatten(), dtype='i')
        z1[(z_cands['0']['p4'].pt.flatten()<z_cands['1']['p4'].pt.flatten())] = 1
        z2[(z_cands['0']['p4'].pt.flatten()<z_cands['1']['p4'].pt.flatten())] = 0
        z1 = JaggedArray.fromoffsets(z_cands.offsets, z1)
        z2 = JaggedArray.fromoffsets(z_cands.offsets, z2)

        passZCand = (z_cands.counts>0)
        output['cutflow']['z cand'] += passZCand.sum()
        selection.add('zCand',passZCand)

        passMassWindow = (passZCand & z_cands[((z_cands.p4.mass>60) & (z_cands.p4.mass<120))].counts>0)
        output['cutflow']['mass window'] += passMassWindow.sum()
        selection.add('massWindow',passMassWindow)

        # im sure there is a better way, but for now just do this
        def get_lepton_values(zl,key):
            val = np.zeros_like(zl.flatten(),dtype=float)
            if len(val)==0:
                return JaggedArray.fromoffsets(zl.offsets,val) 
            for i in range(2):
                mask = (i==zl.flatten())
                if key=='pt':
                    val[mask] = z_cands[passZCand][str(i)].flatten()[mask]['p4'].pt
                elif key=='eta':
                    val[mask] = z_cands[passZCand][str(i)].flatten()[mask]['p4'].eta
                elif key=='phi':
                    val[mask] = z_cands[passZCand][str(i)].flatten()[mask]['p4'].phi
                elif key=='mass':
                    val[mask] = z_cands[passZCand][str(i)].flatten()[mask]['p4'].mass
                else:
                    val[mask] = z_cands[passZCand][str(i)].flatten()[mask][key]
            return JaggedArray.fromoffsets(zl.offsets,val)

        
        z1pt = get_lepton_values(z1,'pt')
        z2pt = get_lepton_values(z2,'pt')
        passPt = ((z1pt>30) & (z2pt>20)).counts>0
        output['cutflow']['pt threshold'] += passPt.sum()
        selection.add('ptThreshold',passPt)


        chanSels = {}
        z1pdg = get_lepton_values(z1,'pdgId')
        z2pdg = get_lepton_values(z2,'pdgId')
        for chan in ['ee','mm']:
            if chan=='ee':
                pdgIds = (11,11)
            if chan=='mm':
                pdgIds = (13,13)
            chanSels[chan] = ((abs(z1pdg)==pdgIds[0])
                            & (abs(z2pdg)==pdgIds[1]))

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
            zls = [z1, z2]
            # electron sf
            for ei, zl in enumerate(zls):
                ei = str(ei)
                eta = get_lepton_values(zl,'etaSC')
                pt = get_lepton_values(zl,'pt')
                electronRecoSF = self._corrections['electron_reco'](eta,pt)
                electronIdSF = self._corrections['electron_id_MVA90'](eta,pt)
                electronSF = np.ones_like(electronRecoSF)
                if ei in ['0','1']:
                    chans = ['ee']
                else:
                    chans = []
                for chan in chans:
                    chanSel = (chanSels[chan].ones_like().sum()>0) # turns empty arrays into 0's, nonempty int 1's
                    electronSF[chanSel] *= electronRecoSF[chanSel].prod()
                    electronSF[chanSel] *= electronIdSF[chanSel].prod()
                weights.add('electronSF'+ei,electronSF)

            # muon SF
            for mi, zl in enumerate(zls):
                mi = str(mi)
                eta = get_lepton_values(zl,'eta')
                pt = get_lepton_values(zl,'pt')
                if self._year == '2016':
                    idSF = self._corrections['muon_id_MediumID'](eta,pt)
                    isoSF = self._corrections['muon_iso_TightRelIso_MediumID'](eta,pt)
                else:
                    idSF = self._corrections['muon_id_MediumPromptID'](pt,abs(eta))
                    isoSF = self._corrections['muon_iso_TightRelIso_MediumID'](pt,abs(eta))
                    
                muonSF = np.ones_like(idSF)
                if mi in ['0','1']:
                    chans = ['mm']
                else:
                    chans = []
                for chan in chans:
                    chanSel = (chanSels[chan].ones_like().sum()>0) # turns empty arrays into 0's, nonempty int 1's
                    muonSF[chanSel] *= idSF[chanSel].prod()
                    muonSF[chanSel] *= isoSF[chanSel].prod()
                weights.add('muonSF'+mi,muonSF)

        logging.debug('filling')
        for sel in self._selections:
            if sel=='massWindow':
                cut = selection.all('lumiMask','trigger','goodVertex','twoLeptons','zCand','massWindow','ptThreshold')
            for chan in ['ee','mm']:
                chanSel = chanSels[chan]
                weight = chanSel.astype(float) * weights.weight()

                output[sel+'_zmass'].fill(
                    dataset=dataset,
                    channel=chan,
                    mass=z_cands[cut].p4.mass.flatten(),
                    weight=weight[cut].flatten(),
                )
                output[sel+'_met'].fill(
                    dataset=dataset,
                    channel=chan,
                    met=df['MET_pt'][cut],
                    weight=weight[cut].flatten(),
                )

        return output


    def postprocess(self, accumulator):
        # always scale to 1000 pb for plotting
        lumi = 1000
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

    years = ['2016','2017','2018']
    for year in years:
        corrections = load(f'corrections/corrections_{year}.coffea')

        processor_instance = DYProcessor(
            year=year,
            corrections=corrections,
        )

        save(processor_instance, f'processors/dyProcessor_{year}.coffea')
