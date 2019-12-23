import json
from utilities import python_mkdir, load, dump

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
        dump(outpath.format(year=year,sample=s),thisfileset)
        datafileset.update(thisfileset)
    dump(outpath.format(year=year,sample='data'),datafileset)
    fileset.update(datafileset)


    mcfileset = {}
    for s in fullmcfileset[year]:
        thisfileset = {s:[]}
        for d in fullmcfileset[year][s]['datasets']:
            thisfileset[s] += fullmcfileset[year][s]['files'][d]
        dump(outpath.format(year=year,sample=s),thisfileset)
        mcfileset.update(thisfileset)
    dump(outpath.format(year=year,sample='mc'),mcfileset)
    fileset.update(mcfileset)

    dump(outpath.format(year=year,sample='all'),fileset)
