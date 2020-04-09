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

logger = logging.getLogger("MMTTProcessor_DevToolsor")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class MMTTProcessor_DevToolsor(processor.ProcessorABC):

    def __init__(self,year='2018',corrections={}):
        self._year = year

        self._corrections = corrections
        self._rochester = lookup_tools.rochester_lookup.rochester_lookup(corrections['rochester_data'])

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        channel_axis = hist.Cat("channel", "Channel")
        mass_axis = hist.Bin("mass", r"$m_{2\ell}$ [GeV]", 600, 0, 60)
        hmass_axis = hist.Bin("mass", r"$m_{4\ell}$ [GeV]", 1200, 0, 1200)
        pt_axis = hist.Bin("pt", r"$p_{T,\ell}$ [GeV]", 3000, 0.25, 300)
        npvs_axis = hist.Bin("npvs", "Number of Vertices", 120, 0, 120)

        self._selections = ['dimuon','ditau','ditauMD']

        hist.Hist.DEFAULT_DTYPE = 'f'  # save some space by keeping float bin counts instead of double
        self._accumulator = processor.dict_accumulator()
        for sel in self._selections:
            self._accumulator[sel + '_mmMass'] = hist.Hist("Counts", dataset_axis, channel_axis, mass_axis)
            self._accumulator[sel + '_ttMass'] = hist.Hist("Counts", dataset_axis, channel_axis, mass_axis)
            self._accumulator[sel + '_mmttMass'] = hist.Hist("Counts", dataset_axis, channel_axis, hmass_axis)
            self._accumulator[sel + '_pileup'] = hist.Hist("Counts", dataset_axis, channel_axis, npvs_axis)

        self._accumulator['cutflow'] = processor.defaultdict_accumulator(int)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(int)
        

    @property
    def accumulator(self):
        return self._accumulator

    def _add_muon_id(self, muons):
        ptCut = (muons.pt>3)
        etaCut = (abs(muons.eta)<2.4)
        dxyCut = (abs(muons.dxy)<0.5)
        dzCut = (abs(muons.dz)<1)
        idCut = (muons.mediumId)

        idNoIsoCut = (ptCut & etaCut & dxyCut & dzCut & idCut)

        isoCut = (muons.pfRelIso04_all<0.15)

        idIsoCut = (idNoIsoCut & isoCut)

        muons['passId'] = idNoIsoCut
        muons['passIdIso'] = idIsoCut

    def _add_ditau_id(self, jets):
        scores = {
            'loose':  0.29463472962379456,
            'medium': 0.8535109758377075,
            'tight':  0.9812706708908081,
            'vtight': 0.9978471994400024,
        }
        scoresMD = {
            'loose':  0.7942289710044861,
            'medium': 0.8037583231925964,
            'tight':  0.8053591251373291,
            'vtight': 0.8073446750640869,
        }
        for wp in scores.keys():
            jets[f'ditau_{wp}'] = (jets.ditau2017v1>scores[wp])
            jets[f'ditauMD_{wp}'] = (jets.ditau2017MDv1>scoresMD[wp])
        
    def _add_trigger(self,df):
        dataset = df['dataset']

        triggerPaths = {}

        # SingleMuon
        if self._year=='2016':
            triggerPaths['SingleMuon'] = [
                "IsoMu24",
                "IsoTkMu24",
            ]
        elif self._year=='2017':
            triggerPaths['SingleMuon'] = [
                "IsoMu27",
            ]
        elif self._year=='2018':
            triggerPaths['SingleMuon'] = [
                "IsoMu24",
            ]

        # Define priority
        # To avoid double counting in data, for a given dataset
        # all current datasets triggers are accepted
        # and all higher datasets triggers are vetoed
        # no lower datasets triggers are looked at
        # in MC, all triggers are accepted
        triggerPriority = [
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
            result = ((result) | (df[p+'Pass'].astype(bool)))
        for p in triggersToVeto:
            result = ((result) & (~df[p+'Pass'].astype(bool)))

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
            df['passLumiMask'] = lumiMask(df['run'],df['lumi'])
        else:
            df['passLumiMask'] = np.ones_like(df['run'],dtype=bool)
        passLumiMask = df['passLumiMask']
        selection.add('lumiMask',passLumiMask)
            

        logging.debug('adding trigger')
        self._add_trigger(df)

        passHLT = df['passHLT']
        selection.add('trigger',passHLT)
        output['cutflow']['pass trigger'] += passHLT.sum()
        # if no trigger: fast return
        if passHLT.sum()==0:
            return output

        
        logging.debug('building muons')
        muons = JaggedCandidateArray.candidatesfromcounts(
            df['muons_count'],
            pt=df['muons_pt'],
            eta=df['muons_eta'],
            phi=df['muons_phi'],
            mass=df['muons_mass'],
            charge=df['muons_charge'],
            dxy=df['muons_dxy'],
            dz=df['muons_dz'],
            mediumId=df['muons_isMediumMuon'],
            pfRelIso04_all=df['muons_relPFIsoDeltaBetaR04'],
        )

        logging.debug('building jets')
        jets = JaggedCandidateArray.candidatesfromcounts(
            df['jets_count'],
            pt=df['jets_pt'],
            eta=df['jets_eta'],
            phi=df['jets_phi'],
            mass=df['jets_mass'],
            ditau2017v1=df['jets_deepDiTau_ditau2017v1'],
            ditau2017MDv1=df['jets_deepDiTau_ditau2017MDv1'],
        )

        logging.debug('adding muon id')
        self._add_muon_id(muons)

        logging.debug('adding dita id')
        self._add_ditau_id(jets)

        logging.debug('selecting muons')
        muonId = (muons.passId>0)
        muons = muons[muonId]

        logging.debug('selecting ditau')
        ditauId = (jets.ditau_loose>0)
        ditauMDId = (jets.ditauMD_loose>0)
        jets = jets[(ditauId | ditauMDId)]

        passTwoMuons = (muons.counts >= 2)
        output['cutflow']['two muons'] += passTwoMuons.sum()
        selection.add('twoMuons',passTwoMuons)

        
        logging.debug('building 2 muons')
        mm_cands = muons.choose(2)
        
        def bestcombination(mmcands):
            good_charge = sum(mmcands[str(i)]['charge'] for i in range(2)) == 0
            # this keeps the first z cand in each event
            # should instead sort the best first
            # TODO: select best
            mmcands = mmcands[good_charge][:,:1]
            return mmcands


        logging.debug('selecting best combinations')
        mm_cands = bestcombination(mm_cands)

        m1 = np.zeros_like(mm_cands['p4'].pt.flatten(), dtype='i')
        m2 = np.ones_like(mm_cands['p4'].pt.flatten(), dtype='i')
        m1[(mm_cands['0']['p4'].pt.flatten()<mm_cands['1']['p4'].pt.flatten())] = 1
        m2[(mm_cands['0']['p4'].pt.flatten()<mm_cands['1']['p4'].pt.flatten())] = 0
        m1 = JaggedArray.fromoffsets(mm_cands.offsets, m1)
        m2 = JaggedArray.fromoffsets(mm_cands.offsets, m2)

        passMMCand = (mm_cands.counts>0)
        output['cutflow']['mm cand'] += passMMCand.sum()
        selection.add('mmCand',passMMCand)

        passMassWindow = (passMMCand & mm_cands[((mm_cands.p4.mass>2.5) & (mm_cands.p4.mass<60))].counts>0)
        output['cutflow']['mass window'] += passMassWindow.sum()
        selection.add('massWindow',passMassWindow)

        # im sure there is a better way, but for now just do this
        def get_lepton_values(ml,key):
            val = np.zeros_like(ml.flatten(),dtype=float)
            if len(val)==0:
                return JaggedArray.fromoffsets(ml.offsets,val) 
            for i in range(2):
                mask = (i==ml.flatten())
                if key=='pt':
                    val[mask] = mm_cands[passMMCand][str(i)].flatten()[mask]['p4'].pt
                elif key=='eta':
                    val[mask] = mm_cands[passMMCand][str(i)].flatten()[mask]['p4'].eta
                elif key=='phi':
                    val[mask] = mm_cands[passMMCand][str(i)].flatten()[mask]['p4'].phi
                elif key=='mass':
                    val[mask] = mm_cands[passMMCand][str(i)].flatten()[mask]['p4'].mass
                else:
                    val[mask] = mm_cands[passMMCand][str(i)].flatten()[mask][key]
            return JaggedArray.fromoffsets(ml.offsets,val)

        
        m1pt = get_lepton_values(m1,'pt')
        m2pt = get_lepton_values(m2,'pt')
        passPt = ((m1pt>26) & (m2pt>3)).counts>0
        output['cutflow']['pt threshold'] += passPt.sum()
        selection.add('ptThreshold',passPt)

        mmtt_allcands = mm_cands.cross(jets)

        def bestcombination_mmtt(mmttcands):
            good_ditau = (mmttcands['1'].ditau_loose>0) & (mmttcands['1'].p4.pt>20) & (abs(mmttcands['1'].p4.eta)<2.5)
            good_ditau = good_ditau & (mmttcands['0'].p4.delta_r(mmttcands['1'].p4)>0.8)
            # this keeps the first cand in each event
            # should instead sort the best first
            # TODO: select best
            mmttcands = mmttcands[good_ditau][:,:1]
            return mmttcands

        def bestcombination_mmttMD(mmttcands):
            good_ditau = (mmttcands['1'].ditauMD_loose>0) & (mmttcands['1'].p4.pt>20) & (abs(mmttcands['1'].p4.eta)<2.5)
            good_ditau = good_ditau & (mmttcands['0'].p4.delta_r(mmttcands['1'].p4)>0.8)
            # this keeps the first cand in each event
            # should instead sort the best first
            # TODO: select best
            mmttcands = mmttcands[good_ditau][:,:1]
            return mmttcands

        mmtt_cands = bestcombination_mmtt(mmtt_allcands)
        mmttMD_cands = bestcombination_mmttMD(mmtt_allcands)

        passMMTTCand = (mmtt_cands.counts>0)
        output['cutflow']['mmtt cand'] += passMMTTCand.sum()
        selection.add('mmttCand',passMMTTCand)

        passMMTTMDCand = (mmttMD_cands.counts>0)
        output['cutflow']['mmttMD cand'] += passMMTTMDCand.sum()
        selection.add('mmttMDCand',passMMTTMDCand)


        weights = processor.Weights(df.size)
        if self._isData: 
            output['sumw'][dataset] = 0 # always set to 0 for data
        else:
            # TODO: this is WRONG for DevTools
            # need to grab the LumiTree and get the sumw from there
            output['sumw'][dataset] += df['genWeight'].sum()
            weights.add('genWeight',df['genWeight'])
            weights.add('pileupWeight',
                        self._corrections['pileupWeight'](df['nTrueVertices']),
                        self._corrections['pileupWeightUp'](df['nTrueVertices']),
                        self._corrections['pileupWeightDown'](df['nTrueVertices']),
                        )
            mls = [m1, m2]

            # muon SF
            for mi, ml in enumerate(mls):
                mi = str(mi)
                eta = get_lepton_values(ml,'eta')
                pt = get_lepton_values(ml,'pt')
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
                    muonSF *= idSF.prod()
                    muonSF *= isoSF.prod()
                weights.add('muonSF'+mi,muonSF)

        logging.debug('filling')
        for sel in self._selections:
            if sel=='dimuon':
                cut = selection.all('lumiMask','trigger','twoMuons','mmCand','massWindow','ptThreshold')
            if sel=='ditau':
                cut = selection.all('lumiMask','trigger','twoMuons','mmCand','massWindow','ptThreshold','mmttCand')
                mmtt = mmtt_cands
            if sel=='ditauMD':
                cut = selection.all('lumiMask','trigger','twoMuons','mmCand','massWindow','ptThreshold','mmttMDCand')
                mmtt = mmttMD_cands
            chan = 'mm'
            weight = weights.weight()

            output[sel+'_mmMass'].fill(
                dataset=dataset,
                channel=chan,
                mass=mm_cands[cut].p4.mass.flatten(),
                weight=weight[cut].flatten(),
            )
            output[sel+'_pileup'].fill(
                dataset=dataset,
                channel=chan,
                npvs=df['vertices_count'][cut],
                weight=weight[cut].flatten(),
            )
            if 'ditau' in sel:
                output[sel+'_ttMass'].fill(
                    dataset=dataset,
                    channel=chan,
                    mass=mmtt[cut]['1'].p4.mass.flatten(),
                    weight=weight[cut].flatten(),
                )
                output[sel+'_mmttMass'].fill(
                    dataset=dataset,
                    channel=chan,
                    mass=mmtt[cut].p4.mass.flatten(),
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

        processor_instance = MMTTProcessor_DevToolsor(
            year=year,
            corrections=corrections,
        )

        save(processor_instance, f'processors/mmttProcessor_DevTools_{year}.coffea')
