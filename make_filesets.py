#!/usr/bin/env python
import os
from collections import Counter
from utilities import load, dump


def maybe_dump(path, fileset):
    # if the file exists: read in the current data, and test for differences
    replace = False
    oldfileset = load(path)
    for dataset in fileset:
        oldfilelist = oldfileset.get(dataset,[])
        filelist = fileset[dataset]
        if Counter(oldfilelist) != Counter(filelist):
            replace = True
    if not replace: return
    dump(path,fileset)

        

fulldatafileset = load('data')
fullmcfileset = load('mc')

years = ['2016','2017','2018']
for year in years:
    outpath = 'filesets/{year}/{sample}'
    fileset = {}

    datafileset = {}
    for s in fulldatafileset[year]:
        thisfileset = {s:[]}
        for d in fulldatafileset[year][s]['datasets']:
            thisfileset[s] += fulldatafileset[year][s]['files'][d]
        maybe_dump(outpath.format(year=year,sample=s),thisfileset)
        datafileset.update(thisfileset)
    maybe_dump(outpath.format(year=year,sample='data'),datafileset)
    fileset.update(datafileset)


    mcfileset = {}
    for s in fullmcfileset[year]:
        thisfileset = {s:[]}
        for d in fullmcfileset[year][s]['datasets']:
            thisfileset[s] += fullmcfileset[year][s]['files'][d]
        maybe_dump(outpath.format(year=year,sample=s),thisfileset)
        mcfileset.update(thisfileset)
    maybe_dump(outpath.format(year=year,sample='mc'),mcfileset)
    fileset.update(mcfileset)

    maybe_dump(outpath.format(year=year,sample='all'),fileset)
