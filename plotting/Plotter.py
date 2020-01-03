import os
import sys
import logging
import math
from array import array
from collections import OrderedDict
import tempfile

import ROOT
ROOT.gROOT.SetBatch(True)

import CMS_lumi
import tdrstyle

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPalette(1)

class Plotter(object):
    '''Basic plotter utilities'''

    def __init__(self,analysis,**kwargs):
        '''Initialize the plotter'''
        # plot directory
        self.analysis = analysis
        #self.intLumi = kwargs.get('intLumi',35920) #2016
        #self.intLumi = kwargs.get('intLumi',41530) #2017
        self.intLumi = kwargs.get('intLumi',59740) #2018
        self.outputDirectory = kwargs.pop('outputDirectory','plots/{0}'.format(self.analysis))
        # initialize stuff
        self.j = 0
        self.sampleStack = []
        self.plotSamples = []
        self.dataSample = 'DATA'
        self.styleMap = {}

    def _getLegend(self,xstart,ystart,xend,yend,numcol=1):
        '''Get the legend'''
        legend = ROOT.TLegend(xstart,ystart,xend,yend,'','NDC')
        if numcol>1: legend.SetNColumns(int(numcol))
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillColor(0)
        return legend

    def _setStyle(self,pad,position=11,preliminary=True,period_int=4,simulation=False):
        '''Set style for plots based on the CMS TDR style guidelines.
           https://twiki.cern.ch/twiki/bin/view/CMS/Internal/PubGuidelines#Figures_and_tables
           https://ghm.web.cern.ch/ghm/plots/'''
        # set period (used in CMS_lumi)
        # period : sqrts
        # 1 : 7, 2 : 8, 3 : 7+8, 4 : 13, ... 7 : 7+8+13
        # set position
        # 11: upper left, 33 upper right
        CMS_lumi.cmsText = 'CMS'
        CMS_lumi.writeExtraText = preliminary
        CMS_lumi.extraText = "Simulation Preliminary" if simulation else 'Preliminary'
        CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (self.intLumi/1000.)
        if self.intLumi < 1000:
            CMS_lumi.lumi_13TeV = "%0.1f pb^{-1}" % (self.intLumi)
        CMS_lumi.CMS_lumi(pad,period_int,position)

    def _save(self, canvas, savename):
        '''Save the canvas in multiple formats.'''
        logging.debug('Saving {0}'.format(savename))
        canvas.SetName(savename)
        for type in ['pdf','root','png']:
            name = '{0}/{1}/{2}.{1}'.format(self.outputDirectory, type, savename)
            directory = os.path.dirname(name)
            if directory: os.makedirs(directory, exist_ok=True)
            logging.debug('Writing {0}'.format(name))
            canvas.Print(name)

    def _saveTemp(self, canvas):
        '''Save the canvas in multiple formats.'''
        temp = tempfile.NamedTemporaryFile(suffix=".png",delete=False)
        canvas.Print(temp.name)
        return temp.name

    def addSampleToStack(self,sample):
        self.sampleStack += [sample]

    def addSampleToPlot(self,sample):
        self.plotSamples += [sample]

    def setDataSample(self,sample):
        self.dataSample = sample

    def setStyleMap(self,styleMap):
        self.styleMap = styleMap
    
    def _getStyle(self,sample):
        return self.styleMap.get(sample,{})

    def _styleHist(self,sample,hist):
        style = self._getStyle(sample)
        if 'linecolor' in style: hist.SetLineColor(style['linecolor'])
        if 'fillcolor' in style: hist.SetFillColor(style['fillcolor'])
        if 'label' in style: hist.SetTitle(style['label'])

    def _getHist(self,sample,hists):
        self.j += 1
        self._styleHist(sample,hists[sample])
        return hists[sample]


    def _getStack(self,hists):
        self.j += 1
        name = 'h_stack_{}'.format(self.j)
        stack = ROOT.THStack(name,name)
        for s in self.sampleStack:
            if hists[s]:
                self._styleHist(s,hists[s])
                stack.Add(hists[s])
        return stack

    def _getData(self,hists):
        self._styleHist(self.dataSample,hists[self.dataSample])
        hist = hists[self.dataSample]
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.)
        hist.SetLineColor(ROOT.kBlack)
        hist.SetBinErrorOption(ROOT.TH1.kPoisson)
        return hist

    def _get_stat_err(self, hist):
        '''Create statistical errorbars froma histogram'''
        self.j += 1
        staterr = hist.Clone("{0}_staterr_{1}".format(hist.GetName,self.j))
        staterr.SetFillColor(ROOT.kGray+3)
        staterr.SetLineColor(ROOT.kGray+3)
        staterr.SetLineWidth(0)
        staterr.SetMarkerSize(0)
        staterr.SetFillStyle(3013)
        return staterr

    def _get_ratio_stat_err(self, hist, **kwargs):
        '''Return a statistical error bars for a ratio plot'''
        ratiomin = kwargs.pop('ratiomin',0.5)
        ratiomax = kwargs.pop('ratiomax',1.5)
        self.j += 1
        ratiostaterr = hist.Clone("{0}_ratiostaterr_{1}".format(hist.GetName,self.j))
        #ratiostaterr.Sumw2()
        ratiostaterr.SetStats(0)
        ratiostaterr.SetTitle("")
        #ratiostaterr.GetYaxis().SetTitle("Data / MC")
        ratiostaterr.GetYaxis().SetTitle("Obs / Exp")
        ratiostaterr.SetMaximum(ratiomax)
        ratiostaterr.SetMinimum(ratiomin)
        ratiostaterr.SetMarkerSize(0)
        ratiostaterr.SetFillColor(ROOT.kGray+3)
        ratiostaterr.SetFillStyle(3013)
        ratiostaterr.GetXaxis().SetLabelSize(0.19)
        ratiostaterr.GetXaxis().SetTitleSize(0.21)
        ratiostaterr.GetXaxis().SetTitleOffset(1.0)
        ratiostaterr.GetXaxis().SetLabelOffset(0.03)
        ratiostaterr.GetYaxis().SetLabelSize(0.19)
        ratiostaterr.GetYaxis().SetLabelOffset(0.006)
        ratiostaterr.GetYaxis().SetTitleSize(0.21)
        ratiostaterr.GetYaxis().SetTitleOffset(0.35)
        ratiostaterr.GetYaxis().SetNdivisions(503)

        # bin by bin errors
        for i in range(hist.GetNbinsX()+2):
            ratiostaterr.SetBinContent(i, 1.0)
            if hist.GetBinContent(i)>1e-6:  # not empty
                binerror = hist.GetBinError(i) / hist.GetBinContent(i)
                ratiostaterr.SetBinError(i, binerror)
            else:
                ratiostaterr.SetBinError(i, 999.)

        return ratiostaterr

    def _get_ratio_err(self,num,denom,data=False):
        '''Return the ratio of two histograms, taking errors from numerator only'''
        if data:
            # get ratio between two hists with poisson errors
            num_graph = self._get_poisson_err(num)
            graph = ROOT.TGraphAsymmErrors(num.GetNbinsX())
            npoints = 0
            for bin in range(num.GetNbinsX()):
                nval = num.GetBinContent(bin+1)
                dval = denom.GetBinContent(bin+1)
                ey_low = num_graph.GetErrorYlow(bin)
                ey_high = num_graph.GetErrorYhigh(bin)
                if dval > 0:
                    graph.SetPoint(npoints, num.GetBinCenter(bin+1), nval/dval)
                    graph.SetPointEXlow(npoints, 0)
                    graph.SetPointEXhigh(npoints, 0)
                    graph.SetPointEYlow(npoints, ey_low/dval)
                    graph.SetPointEYhigh(npoints, ey_high/dval)
                else:
                    graph.SetPoint(npoints, num.GetBinCenter(bin+1), 0)
                    graph.SetPointEXlow(npoints, 0)
                    graph.SetPointEXhigh(npoints, 0)
                    graph.SetPointEYlow(npoints, 0)
                    graph.SetPointEYhigh(npoints, 0)
                npoints += 1
            graph.Set(npoints)
            return graph
        else:
            self.j += 1
            newnum = num.Clone('ratio_{0}'.format(self.j))
            for b in range(num.GetNbinsX()):
                nVal = num.GetBinContent(b+1)
                nErr = num.GetBinError(b+1)
                dVal = denom.GetBinContent(b+1)
                #if dVal>1e-6: # why?
                if dVal>0:
                    val = (nVal+dVal)/dVal
                    err = nErr/dVal
                else:
                    val = 0
                    err = 0
                newnum.SetBinContent(b+1,val)
                newnum.SetBinError(b+1,err)
            return newnum

    def _get_poisson_err(self,hist):
        # adapted from rootpy to get asymmetric poisson errors
        graph = ROOT.TGraphAsymmErrors(hist.GetNbinsX())
        chisqr = ROOT.TMath.ChisquareQuantile
        npoints = 0
        for bin in range(hist.GetNbinsX()):
            val = hist.GetBinContent(bin+1)
            err2 = hist.GetBinError(bin+1)**2
            width = hist.GetBinWidth(bin+1)
            varbin = val-err2>0.001
            entries = val*width if varbin else val
            if entries <= 0:
                #continue
                entries = 0
            ey_low = entries - 0.5 * chisqr(0.1586555, 2. * entries)
            ey_high = 0.5 * chisqr(1. - 0.1586555, 2. * (entries + 1)) - entries

            if varbin:
                ey_low = val - (entries-ey_low) / width
                ey_high = (entries+ey_high) / width - val

            ex = width / 2.
            graph.SetPoint(npoints, hist.GetBinCenter(bin+1), val)
            graph.SetPointEXlow(npoints, 0)
            graph.SetPointEXhigh(npoints, 0)
            graph.SetPointEYlow(npoints, ey_low)
            graph.SetPointEYhigh(npoints, ey_high)
            npoints += 1
        graph.Set(npoints)
        return graph

    def plot(self,hists,savename,**kwargs):
        xlabel = kwargs.get('xlabel','Variable')
        ylabel = kwargs.get('ylabel','Events')
        logx = kwargs.get('logx',False)
        logy = kwargs.get('logy',False)
        ymax = kwargs.get('ymax',None)
        ymin = kwargs.get('ymin',None)
        blind = kwargs.get('blind',False)


        plotratio = True

        logging.info('Plotting {0}'.format(savename))

        ROOT.gDirectory.Delete('h_*')

        canvas = ROOT.TCanvas(savename,savename,50,50,600,600)

        # ratio plot
        if plotratio:
            plotpad = ROOT.TPad("plotpad", "top pad", 0.0, 0.21, 1.0, 1.0)
            plotpad.SetBottomMargin(0.04)
            plotpad.SetRightMargin(0.03)
            plotpad.Draw()
            plotpad.SetLogy(logy)
            plotpad.SetLogx(logx)
            ratiopad = ROOT.TPad("ratiopad", "bottom pad", 0.0, 0.0, 1.0, 0.21)
            ratiopad.SetTopMargin(0.06)
            ratiopad.SetRightMargin(0.03)
            ratiopad.SetBottomMargin(0.5)
            ratiopad.SetLeftMargin(0.16)
            ratiopad.SetTickx(1)
            ratiopad.SetTicky(1)
            ratiopad.Draw()
            ratiopad.SetLogx(logx)
            if plotpad != ROOT.TVirtualPad.Pad(): plotpad.cd()
        else:
            canvas.SetRightMargin(0.04)
            canvas.SetLogy(logy)
            canvas.SetLogx(logx)

        legend = self._getLegend(0.6,0.5,0.94,0.9)

        # stack
        stack = self._getStack(hists)
        stack.Draw('hist')
        stack.GetXaxis().SetTitle(xlabel)
        stack.GetYaxis().SetTitle(ylabel)
        stack.SetMaximum(stack.GetMaximum()*5 if logy else stack.GetMaximum()*1.2)
        if ymax!=None: stack.SetMaximum(ymax)
        if ymin!=None: stack.SetMinimum(ymin)
        if logx:
            stack.GetXaxis().SetMoreLogLabels()
            stack.GetXaxis().SetNoExponent()
        if plotratio: stack.GetHistogram().GetXaxis().SetLabelOffset(999)
        # stat err
        self.j += 1
        staterr = self._get_stat_err(stack.GetStack().Last().Clone('h_stack_{0}'.format(self.j)))
        staterr.Draw('e2 same')
        # syst err
        #
        for hist in reversed(stack.GetHists()):
            legend.AddEntry(hist,hist.GetTitle(),'f')

        # data
        if not blind:
            data = self._getData(hists)
            data.Draw('ex0 same') # TH1
            #data.Draw('p0 same') # TGraph
            legend.AddEntry(data,data.GetTitle(),'ep')

        # overlay histograms
        sampleHists = OrderedDict()
        for sample in self.plotSamples:
            sampleHists[sample] = self._getHist(sample,hists)
            sampleHists[sample].Draw('hist same')
            sampleHists[sample].SetLineWidth(3)
            legend.AddEntry(sampleHists[sample],sampleHists[sample].GetTitle(),'l')

        legend.Draw()

        if plotratio:
            self._setStyle(plotpad)

        if plotratio:
            self.j += 1
            stackname = 'h_stack_{0}_ratio'.format(self.j)
            denom = stack.GetStack().Last().Clone(stackname)
            ratiostaterr = self._get_ratio_stat_err(denom,**kwargs)
            ratiostaterr.SetXTitle(xlabel)
            ratios = OrderedDict()
            if not blind:
                histName = 'DATA'
                numname = 'h_{0}_{1}_ratio'.format(histName,self.j)
                num = data.Clone(numname)
                numratio = self._get_ratio_err(num,denom,data=histName==self.dataSample)
                ratios[histName] = numratio
            for histName in self.plotSamples:
                numname = 'h_{0}_{1}_ratio'.format(histName,self.j)
                num = sampleHists[sample].Clone(numname)
                numratio = self._get_ratio_err(num,denom,data=False)
                ratios[histName] = numratio

            # and draw
            if ratiopad != ROOT.TVirtualPad.Pad(): ratiopad.cd()
            ratiostaterr.Draw("e2")
            if logx:
                ratiostaterr.GetXaxis().SetMoreLogLabels()
                ratiostaterr.GetXaxis().SetNoExponent()
            unityargs = [ratiostaterr.GetXaxis().GetXmin(),1,ratiostaterr.GetXaxis().GetXmax(),1]
            ratiounity = ROOT.TLine(*unityargs)
            ratiounity.SetLineStyle(2)
            ratiounity.Draw('same')
            for histName, hist in ratios.items():
                if histName==self.dataSample:
                    #hist.Draw('e0 same')
                    hist.Draw('0P same')
                else:
                    hist.SetLineWidth(3)
                    hist.Draw('hist same')
        else:
            self._setStyle(canvas)

        # ratio

        self._save(canvas,savename)



