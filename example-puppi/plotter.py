#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

from array import array

import ROOT

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);
#ROOT.gStyle.SetPadRightMargin(0.16);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetErrorX(0.5);

def makeCanvas(hists, names, canname, type=1, isLog=False):
    
    directory = "plots";
    colors  = [1,4,2,6,7,3,5];
    markers = [20,24,21,25,22,26,23];
    
    max = -999.;
    for hist in hists:
        if max < hist.GetMaximum(): max = hist.GetMaximum();
    
    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for i in range(len(names)):
        hists[i].SetLineColor(colors[i]);
        hists[i].SetMarkerColor(colors[i]);
        hists[i].SetMarkerStyle(markers[i]);
        hists[i].SetLineWidth(2);
        leg.AddEntry(hists[i], names[i], "flp");
    if type == 0:
        hists[0].SetFillStyle(3003);
        hists[0].SetFillColor(colors[0]);

    banner = ROOT.TLatex(0.18,0.92,("Dawg Pound Working Group, dijets"));
    banner.SetNDC()
    banner.SetTextSize(0.035)

    can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
    hists[0].SetMaximum( 1.2*max );
    hists[0].SetMinimum( 0 );
    hists[0].Draw();
    for i in range(1,len(hists)):
        if type == 0: hists[i].Draw("histsames");
        elif type == 1: hists[i].Draw("sames");
        else: print "problem!"

    leg.Draw();
    if isLog:
        ROOT.gPad.SetLogy();
        hists[0].SetMinimum( 1 );
    banner.Draw();
    can.SaveAs(directory+"/"+canname+".eps");
    can.SaveAs(directory+"/"+canname+".png");



############ --------------------------------------------------------
############ --------------------------------------------------------
############ --------------------------------------------------------

if __name__ == '__main__':
    
    file = ROOT.TFile("Output.root");
    tree = file.Get("Tree");

    # define the histograms
    hists_name = []; hists_vars = []; hists_comp = []; hists_pars = []; hists_axis = []; hists_type = [];
    h_name = "pts";
    h_vars = ["Genpt","Puppiptraw","PFpt","CHSpt"];
    h_comp = None;
    h_pars = [20,0,300];
    h_axis = "; pT (GeV); N";
    h_type = 0;
    hists_name.append(h_name); hists_vars.append(h_vars); hists_comp.append(h_comp); hists_pars.append(h_pars); hists_axis.append(h_axis); hists_type.append(h_type);
    
    h_name = "mass";
    h_vars = ["Genm","Puppimraw","PFm","CHSm"];
    h_comp = None;
    h_pars = [20,0,100];
    h_axis = "; mass (GeV); N";
    h_type = 0;
    hists_name.append(h_name); hists_vars.append(h_vars); hists_comp.append(h_comp); hists_pars.append(h_pars); hists_axis.append(h_axis); hists_type.append(h_type);
    
    h_name = "massRes";
    h_vars = ["Puppimraw","PFm","CHSm"];
    h_comp = "Genm";
    h_pars = [20,-100,100];
    h_axis = "; mass - massHS (GeV); N";
    h_type = 1;
    hists_name.append(h_name); hists_vars.append(h_vars); hists_comp.append(h_comp); hists_pars.append(h_pars); hists_axis.append(h_axis); hists_type.append(h_type);
    
    # make the histograms
    hists = [];
    for b in range(len(hists_vars)):
        tmphists = [];
        for i in range(len(hists_vars[b])):
            print b,i
            print hists_name[b]
            print hists_axis[b]
            print hists_vars[b][i]
            tmphists.append( ROOT.TH1F(hists_name[b]+"_"+hists_vars[b][i],hists_axis[b],hists_pars[b][0],hists_pars[b][1],hists_pars[b][2]) );
        hists.append(tmphists);

    # fill the histograms
    for a in range(tree.GetEntriesFast()):
        
        tree.GetEntry(a);
        
        if not tree.index == 0: continue;
        
        for b in range(len(hists)):
            for i in range(len(hists[b])):
                if hists_comp[b] == None: hists[b][i].Fill( getattr(tree,hists_vars[b][i]) );
                else: hists[b][i].Fill( getattr(tree,hists_vars[b][i]) -  getattr(tree,hists_comp[b]) );

    # plot the histograms
    for b in range(len(hists)):
        makeCanvas( hists[b], hists_vars[b], hists_name[b], hists_type[b] );





