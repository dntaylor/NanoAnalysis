#!/usr/bin/env python
import numpy as np
import uproot
from coffea import hist, lookup_tools
from coffea.util import load, save
from xsec import xsec

def save_corrections(year):
    corrections = {}
    
    # cross sections
    corrections['xsec'] = xsec
    # manually add the test samples
    corrections['xsec']['DY'] = 6077.22
    corrections['xsec']['HZZ'] = 43.92 * 2.64e-02 * (3.3658e-2*3)**2
    corrections['xsec']['DoubleMuon'] = 1.
    
    extractor = lookup_tools.extractor()
    # electron
    extractor.add_weight_sets([
        # electron reco
        f'electron_reco_ * data/scalefactors/electron/Ele_Reco_{year}.root',
        # electron hzz id
        f'electron_hzz_id_nogap_ * data/scalefactors/electron/ElectronSF_Legacy_{year}_NoGap.root',
        f'electron_hzz_id_gap_ * data/scalefactors/electron/ElectronSF_Legacy_{year}_Gap.root',
    ])
    extractor.finalize()
    evaluator = extractor.make_evaluator()
    
    corrections['electron_reco'] = evaluator['electron_reco_EGamma_SF2D']
    corrections['electron_hzz_id_nogap'] = evaluator['electron_hzz_id_nogap_EGamma_SF2D']
    corrections['electron_hzz_id_gap'] = evaluator['electron_hzz_id_gap_EGamma_SF2D']
    
    # pileup
    with uproot.open(f'data/pileup/dataPileup{year}.root') as f:
        norm = lambda x: x/x.sum()
        edges = f['pileup'].edges
        dataPileup = norm(f['pileup'].values)
        dataPileupUp = norm(f['pileup_plus'].values)
        dataPileupDown = norm(f['pileup_minus'].values)
    with uproot.open(f'data/pileup/mcPileup{year}.root') as f:
        mcPileup = f['pu_mc'].values
    def zeropad(a,n):
        _a = np.zeros(n)
        _a[:len(a)] = a
        return _a
    nmax = max(len(dataPileup),len(mcPileup))
    dataPileup = zeropad(dataPileup,nmax)
    mcPileup = zeropad(mcPileup,nmax)
    mask = (mcPileup>0)
    pileupRatio = dataPileup.copy()
    pileupRatioUp = dataPileupUp.copy()
    pileupRatioDown = dataPileupDown.copy()
    pileupRatio[mask] /= mcPileup[mask]
    pileupRatioUp[mask] /= mcPileup[mask]
    pileupRatioDown[mask] /= mcPileup[mask]
    
    corrections[f'pileupWeight{year}'] = lookup_tools.dense_lookup.dense_lookup(pileupRatio, edges)
    corrections[f'pileupWeight{year}Up'] = lookup_tools.dense_lookup.dense_lookup(pileupRatioUp, edges)
    corrections[f'pileupWeight{year}Down'] = lookup_tools.dense_lookup.dense_lookup(pileupRatioDown, edges)
    
    save(corrections, f'corrections_{year}.coffea')

for year in ['2016','2017','2018']:
    save_corrections(year)
