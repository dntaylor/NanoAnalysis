#!/usr/bin/env python
from __future__ import print_function, division
import os
import argparse

import uproot
import numpy as np

from coffea import hist
from coffea.util import load, save

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert coffea histogram file')
    parser.add_argument('fnames', type=str, nargs='+', help='Histogram files')
    args = parser.parse_args()

    for fname in args.fnames:
        print(f'Converting {fname}')
        outname = fname.replace('.coffea','.root')
        hists = load(fname)

        if os.path.exists(outname):
            os.remove(outname)
        fout = uproot.create(outname)
        
        if isinstance(hists,tuple):
            hists = hists[0]
        
        for key,h in hists.items():
            if not isinstance(h,hist.Hist): continue
            for dataset in h.identifiers('dataset'):
                for channel in h.identifiers('channel'):
                    newhist = h.integrate('dataset',dataset).integrate('channel',channel)
                    hname = '{}_{}_{}'.format(dataset,channel,key)
                    fout[hname] = hist.export1d(newhist)
        
        fout.close()
        
