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

logger = logging.getLogger("MMTTProcessor")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class MMTTProcessor(processor.ProcessorABC):

    def __init__(self,year='2018',corrections={}):
        self._year = year

        self._corrections = corrections
        self._rochester = lookup_tools.rochester_lookup.rochester_lookup(corrections['rochester_data'])

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        channel_axis = hist.Cat("channel", "Channel")
        amass_axis = hist.Bin("mass", r"$m_{2\ell}$ [GeV]", 6500, 0, 65)
        hmass_axis = hist.Bin("mass", r"$m_{4\ell}$ [GeV]", 1500, 0, 1500)
        pt_axis = hist.Bin("pt", r"$p_{T,\ell}$ [GeV]", 3000, 0.25, 300)
        npvs_axis = hist.Bin("npvs", "Number of Vertices", 120, 0, 120)

        self._selections = ['iso','antiiso']

        hist.Hist.DEFAULT_DTYPE = 'f'  # save some space by keeping float bin counts instead of double
        self._accumulator = processor.dict_accumulator()
        for sel in self._selections:
            self._accumulator['_'.join([sel, 'ammmass'])] = hist.Hist("Counts", dataset_axis, channel_axis, amass_axis)
            self._accumulator['_'.join([sel, 'attmass'])] = hist.Hist("Counts", dataset_axis, channel_axis, amass_axis)
            self._accumulator['_'.join([sel, 'hmass'])] = hist.Hist("Counts", dataset_axis, channel_axis, hmass_axis)
            self._accumulator['_'.join([sel, 'pileup'])] = hist.Hist("Counts", dataset_axis, channel_axis, npvs_axis)

        self._accumulator['cutflow'] = processor.defaultdict_accumulator(int)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(int)
        

    @property
    def accumulator(self):
        return self._accumulator

    def _add_muon_id(self, muons):
        # note: input muons must pass
        # slimmedMuons and (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))
        ptCut = (muons.pt>3)
        etaCut = (abs(muons.eta)<2.4)
        dxyCut = (abs(muons.dxy)<0.5)
        dzCut = (abs(muons.dz)<1)
        idCut = (muons.looseId)

        idNoIsoCut = (ptCut & etaCut & dxyCut & dzCut & idCut)

        isoCut = (muons.pfRelIso04_all<0.25)

        idIsoCut = (idNoIsoCut & isoCut)

        muons['passId'] = idNoIsoCut
        muons['passIso'] = isoCut
        muons['passIdIso'] = idIsoCut

    def _add_electron_id(self, electrons):
        # note: input electrons must pass
        # slimmedElectrons and (pt > 5)
        ptCut = (electrons.pt>5)
        etaCut = (abs(electrons.eta)<2.5)
        dxyCut = (abs(electrons.dxy)<0.5)
        dzCut = (abs(electrons.dz)<1)
        idCut = (electrons.mvaFall17V2noIso_WP90)
        isoCut = (electrons.pfRelIso03_all<0.25)

        loose = (ptCut & etaCut & dxyCut & dzCut & idCut)
        looseIso = (loose & isoCut)

        electrons['passId'] = loose
        electrons['passIso'] = isoCut
        electrons['passIdIso'] = looseIso
        
    def _add_trigger(self,df):
        dataset = df['dataset']

        triggerPaths = {}

        # SingleMuon
        if self._year=='2016':
            triggerPaths['SingleMuon'] = [
                "HLT_IsoMu24",
                "HLT_IsoTkMu24",
            ]
        elif self._year=='2017':
            triggerPaths['SingleMuon'] = [
                #"HLT_IsoMu24", # TODO: partially prescaled, check if lower pt threshold is better
                "HLT_IsoMu27",
            ]
        elif self._year=='2018':
            triggerPaths['SingleMuon'] = [
                "HLT_IsoMu24",
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
        # if no trigger: fast return
        if passHLT.sum()==0:
            return output

        
        # require one good vertex
        logging.debug('checking vertices')
        passGoodVertex = (df['PV_npvsGood']>0)
        output['cutflow']['good vertex'] += passGoodVertex.sum()
        selection.add('goodVertex',passGoodVertex)


        # run rochester
        rochester = self._rochester
        _muon_offsets = JaggedArray.counts2offsets(df['nMuon'])
        _charge = JaggedArray.fromoffsets(_muon_offsets, df['Muon_charge'])
        _pt     = JaggedArray.fromoffsets(_muon_offsets, df['Muon_pt'])
        _eta    = JaggedArray.fromoffsets(_muon_offsets, df['Muon_eta'])
        _phi    = JaggedArray.fromoffsets(_muon_offsets, df['Muon_phi'])
        if self._isData:
            _k = rochester.kScaleDT(_charge,_pt,_eta,_phi)
            #_kErr = rochester.kScaleDTerror(_charge,_pt,_eta,_phi)
        else:
            # for default if gen present
            _gen_offsets = JaggedArray.counts2offsets(df['nGenPart'])
            _gid = JaggedArray.fromoffsets(_muon_offsets, df['Muon_genPartIdx'])
            _gpt = JaggedArray.fromoffsets(_gen_offsets, df['GenPart_pt'])
            # for backup w/o gen
            _nl = JaggedArray.fromoffsets(_muon_offsets, df['Muon_nTrackerLayers'])
            _u  = JaggedArray.fromoffsets(_muon_offsets, np.random.rand(*_pt.flatten().shape))
            _hasgen = (_gid>=0)
            _kspread = rochester.kSpreadMC(_charge[ _hasgen], _pt[ _hasgen], _eta[ _hasgen], _phi[ _hasgen], _gpt[_gid[_hasgen]])
            _ksmear  = rochester.kSmearMC( _charge[~_hasgen], _pt[~_hasgen], _eta[~_hasgen], _phi[~_hasgen], _nl[~_hasgen], _u[~_hasgen])
            _k = np.ones_like(_pt.flatten())
            _k[_hasgen.flatten()] = _kspread.flatten()
            _k[~_hasgen.flatten()] = _ksmear.flatten()
            _k = JaggedArray.fromoffsets(_muon_offsets, _k)
            #_kErrspread = rochester.kSpreadMCerror(_charge[ _hasgen], _pt[ _hasgen], _eta[ _hasgen], _phi[ _hasgen], _gpt[_gid[_hasgen]])
            #_kErrsmear  = rochester.kSmearMCerror( _charge[~_hasgen], _pt[~_hasgen], _eta[~_hasgen], _phi[~_hasgen], _nl[~_hasgen], _u[~_hasgen])
            #_kErr = np.ones_like(_pt.flatten())
            #_kErr[_hasgen.flatten()] = _kErrspread.flatten()
            #_kErr[~_hasgen.flatten()] = _kErrsmear.flatten()
            #_kErr = JaggedArray.fromoffsets(_muon_offsets, _kErr)

        mask = _pt.flatten() < 200
        rochester_pt = _pt.flatten()
        rochester_pt[mask] = (_k * _pt).flatten()[mask]

        logging.debug('building muons')
        muons = JaggedCandidateArray.candidatesfromcounts(
            df['nMuon'],
            pt=rochester_pt,
            eta=df['Muon_eta'],
            phi=df['Muon_phi'],
            mass=df['Muon_mass'],
            charge=df['Muon_charge'],
            dxy=df['Muon_dxy'],
            dz=df['Muon_dz'],
            looseId=df['Muon_looseId'],
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
            mvaFall17V2noIso_WP90=df['Electron_mvaFall17V2noIso_WP90'],
            pfRelIso03_all=df['Muon_pfRelIso03_all'],
            pdgId=df['Electron_pdgId'],
        )

        logging.debug('adding muon id')
        self._add_muon_id(muons)

        logging.debug('adding electron id')
        self._add_electron_id(electrons)

        logging.debug('selecting muons')
        muonId = (muons.passId>0)
        muons = muons[muonId]

        passTwoMuons = (muons.counts >= 2)
        output['cutflow']['two muons'] += passTwoMuons.sum()
        selection.add('twoMuons',passTwoMuons)

        
        logging.debug('building 2 muons')
        mm_cands = muons.choose(2)

        # TODO: figure this out....
        # for now, lets assume it is pt sorted already
        def reorder(cands,*ls):
            for ni in range(len(ls)):
                for oi in range(len(ls)):
                    if oi<=ni: continue
                    mask = (ls[oi].content==ni)
                    cands[str(ni)].content[mask], cands[str(oi)].content[mask] = cands[str(oi)].content[mask], cands[str(ni)].content[mask]
            return cands
        #mm_cands = reorder(mm_cands,l1,l2)
        
        def bestcombination(cands):
            # OS
            good_charge = sum(cands[str(i)]['charge'] for i in range(2)) == 0
            cands = cands[good_charge]

            # loose mass window
            mass_window = ((cands['p4'].mass>2.5) & (cands['p4'].mass<62.5))
            cands = cands[mass_window]

            # isolate l1
            l1_iso = (cands['0'].passIso > 0)
            cands = cands[l1_iso]

            # pt threshold
            l1pt = cands['0']['p4'].pt
            l2pt = cands['1']['p4'].pt
            # TODO: special handling for 2017, remove if we find lower threshold is better
            if self._year == '2017':
                pass_pt = ((l1pt>29) & (l2pt>3))
            else:
                pass_pt = ((l1pt>26) & (l2pt>3))
            cands = cands[pass_pt]

            # match l1 to trigger
            # for muons: TrigObj_id == 13
            # for IsoMu: TrigObj_filterBits && 1 << 1 && 1 << 3 # 0 = TrkIsoVVL, 1 = Iso, 3 = 1mu, 10 = 1mu (Mu50), 11 = 1mu (Mu100)
            # will default mass to muon, since filter anyway
            TrigObj_mass = np.ones_like(df['TrigObj_pt']) * 0.1057
            hltmuons = JaggedCandidateArray.candidatesfromcounts(
                df['nTrigObj'],
                pt=df['TrigObj_pt'],
                eta=df['TrigObj_eta'],
                phi=df['TrigObj_phi'],
                mass=TrigObj_mass,
                id=df['TrigObj_id'],
                filterBits=df['TrigObj_filterBits'],
            )
            hltmuons = hltmuons[(hltmuons.id == 13)]
            hltmuons = hltmuons[(hltmuons.filterBits & (1 << 1) & (1 << 3))] # IsoMu
            l1p4 = cands['0']['p4']
            hltp4 = hltmuons['p4']
            l1_hlt = l1p4.cross(hltp4, nested=True)
            matched_hlt = (l1_hlt['0'].delta_r(l1_hlt['1']) < 0.1).any()
            cands = cands[matched_hlt]

            # this keeps the first cand in each event
            # should instead sort the best first
            # TODO: select best
            cands = cands[:,:1]
            return cands


        logging.debug('selecting best combinations')
        mm_cands = bestcombination(mm_cands)

        passMMCand = (mm_cands.counts>0)
        output['cutflow']['mm cand'] += passMMCand.sum()
        selection.add('mmCand',passMMCand)

        # anti iso not working?
        passMMIso = (mm_cands['1'].passIso>0).counts>0
        output['cutflow']['iso'] += passMMIso.sum()
        selection.add('iso',passMMIso)
        selection.add('antiiso',~passMMIso)

        weights = processor.Weights(df.size)
        if self._isData: 
            output['sumw'][dataset] = 0 # always set to 0 for data
        else:
            output['sumw'][dataset] += df['genWeight'].sum()
            weights.add('genWeight',df['genWeight'])
            weights.add('pileupWeight',
                        self._corrections['pileupWeight'](df['Pileup_nPU']),
                        self._corrections['pileupWeightUp'](df['Pileup_nPU']),
                        self._corrections['pileupWeightDown'](df['Pileup_nPU']),
                        )
            for mi in range(2):
                mi = str(mi)
                eta = mm_cands[passMMCand][mi]['p4'].eta
                pt = mm_cands[passMMCand][mi]['p4'].pt
                if self._year == '2016':
                    idSF = self._corrections['muon_id_MediumID'](eta,pt)
                    isoSF = self._corrections['muon_iso_TightRelIso_MediumID'](eta,pt)
                else:
                    idSF = self._corrections['muon_id_MediumPromptID'](pt,abs(eta))
                    isoSF = self._corrections['muon_iso_TightRelIso_MediumID'](pt,abs(eta))
                    
                muonSF = np.ones_like(df.size)
                chan = 'mm'
                muonSF[passMMCand] *= idSF.prod()
                muonSF[passMMCand] *= isoSF.prod()
                weights.add('muonSF'+mi,muonSF)

        logging.debug('filling')
        for sel in self._selections:
            cutSels = ['lumiMask','trigger','goodVertex','twoLeptons','mmCand']
            if sel=='iso':
                cut = selection.all(*cutSels+['iso'])
            if sel=='antiiso':
                cut = selection.all(*cutSels+['antiiso'])
            chan = 'mm'
            weight = weights.weight()

            output['_'.join([sel,'mass'])].fill(
                dataset=dataset,
                channel=chan,
                mass=mm_cands[cut].p4.mass.flatten(),
                weight=weight[cut].flatten(),
            )
            output['_'.join([sel,'pileup'])].fill(
                dataset=dataset,
                channel=chan,
                npvs=df['PV_npvs'][cut],
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

        processor_instance = MMTTProcessor(
            year=year,
            corrections=corrections,
        )

        save(processor_instance, f'processors/mmttProcessor_{year}.coffea')
