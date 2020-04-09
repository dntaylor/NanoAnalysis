#!/usr/bin/env python
from __future__ import print_function, division
import os
import sys
import logging
import argparse
import numpy as np
import uproot
import uproot_methods
from functools import partial
from coffea import hist, processor, lookup_tools
from coffea.lumi_tools import lumi_tools
from coffea.util import load, save
from coffea.analysis_objects import JaggedCandidateArray
from coffea.nanoaod import NanoEvents, NanoCollection
from awkward import JaggedArray, IndexedArray

ZMASS = 91.1876

logger = logging.getLogger("MMTTGenProcessor")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class MMTTGenProcessor(processor.ProcessorABC):
    # will sync with
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4lRunIILegacy
    # to the best of NanoAOD ability
    def __init__(self,year='2018',corrections={}):
        self._year = year

        self._corrections = corrections

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        channel_axis = hist.Cat("channel", "Channel")
        amass_axis = hist.Bin("mass", r"$m_{2\ell}$ [GeV]", 6500, 0, 65)
        hmass_axis = hist.Bin("mass", r"$m_{4\ell}$ [GeV]", 1500, 0, 1500)
        pt_axis = hist.Bin("pt", r"$p_{T,\ell}$ [GeV]", 3000, 0.25, 300)

        hist.Hist.DEFAULT_DTYPE = 'f'  # save some space by keeping float bin counts instead of double
        self._accumulator = processor.dict_accumulator()

        self._accumulator['ammmass'] = hist.Hist("Counts", dataset_axis, channel_axis, amass_axis)
        self._accumulator['attmass'] = hist.Hist("Counts", dataset_axis, channel_axis, amass_axis)
        self._accumulator['hmass'] = hist.Hist("Counts", dataset_axis, channel_axis, hmass_axis)
        self._accumulator['ammvismass'] = hist.Hist("Counts", dataset_axis, channel_axis, amass_axis)
        self._accumulator['attvismass'] = hist.Hist("Counts", dataset_axis, channel_axis, amass_axis)
        self._accumulator['hvismass'] = hist.Hist("Counts", dataset_axis, channel_axis, hmass_axis)

        self._accumulator['cutflow'] = processor.defaultdict_accumulator(int)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(int)
        

    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        logging.debug('starting process')
        output = self.accumulator.identity()

        # need to set the mass for SM particles
        #masses = {
        #    11: 0.510998e-3,
        #    13: 0.105658,
        #    15: 1.77686,
        #}
        #for pdgId in masses:
        #    events.GenPart.mass[abs(events.GenPart.pdgId)==pdgId] = masses[pdgId]

        dataset = events.metadata['dataset']
        self._isData = dataset in ['SingleMuon','DoubleMuon','SingleElectron','DoubleEG','EGamma','MuonEG']
        if self._isData:
            return output

        selection = processor.PackedSelection()

        output['cutflow']['all events'] += events.event.size

        gen = events.GenPart

        #print(gen.pdgId)
        #print(gen.mass)

        # status flags
        # 0 : isPrompt, 
        # 1 : isDecayedLeptonHadron, 
        # 2 : isTauDecayProduct, 
        # 3 : isPromptTauDecayProduct, 
        # 4 : isDirectTauDecayProduct, 
        # 5 : isDirectPromptTauDecayProduct, 
        # 6 : isDirectHadronDecayProduct, 
        # 7 : isHardProcess, 
        # 8 : fromHardProcess, 
        # 9 : isHardProcessTauDecayProduct, 
        #10 : isDirectHardProcessTauDecayProduct, 
        #11 : fromHardProcessBeforeFSR, 
        #12 : isFirstCopy, 
        #13 : isLastCopy, 
        #14 : isLastCopyBeforeFSR,

        #print(gen.pdgId[0])
        #print(np.vectorize(partial(np.binary_repr,width=16))(gen.statusFlags[0]))
        #print(gen.children.pdgId[0])

        isH = (((gen.pdgId==25) | (gen.pdgId==35)) & gen.hasFlags(['isLastCopy']))
        isA = ((gen.pdgId==36) & gen.hasFlags(['isLastCopy']))
        genh = gen[isH]
        gena = gen[isA]

        print(isH[0])
        print(genh.pdgId[0])
        print(genh.children.pdgId[0])
        print(isA[0])
        print(gena.pdgId[0])
        print(gena.children.pdgId[0])

        genamm = gena[(abs(gena.children.pdgId)==13).all()]
        genatt = gena[(abs(gena.children.pdgId)==15).all()]

        # get gen channel
        def get_stable_children(g):
            children = g.children
            stable = (children.status==1)
            stable_children = children[stable]
            get_stable_children(children[~stable])
        

        #output['hmass'].fill(
        #    dataset=dataset,
        #    channel=chan,
        #    mass=genh.p4.mass.flatten(),
        #)

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

        processor_instance = MMTTGenProcessor(
            year=year,
            corrections=corrections,
        )

        save(processor_instance, f'processors/mmttGenProcessor_{year}.coffea')
