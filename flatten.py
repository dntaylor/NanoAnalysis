# coding: utf-8
from __future__ import print_function, division
import os

import uproot
import numpy as np

from coffea import hist
from coffea.util import load, save

hists = load('hists.coffea')


if os.path.exists("templates.root"):
    os.remove("templates.root")
fout = uproot.create("templates.root")

if isinstance(hists,tuple):
    hists = hists[0]

for key,h in hists.items():
    print(key)
    if not isinstance(h,hist.Hist): continue
    for dataset in h.identifiers('dataset'):
        print(dataset)
        for channel in h.identifiers('channel'):
            print(channel)
            newhist = h.integrate('dataset',dataset).integrate('channel',channel)
            hname = '{}_{}_{}'.format(dataset,channel,key)
            fout[hname] = hist.export1d(newhist)

fout.close()

