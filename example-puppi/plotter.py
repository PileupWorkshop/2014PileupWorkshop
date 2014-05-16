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
ROOT.gStyle.SetPadRightMargin(0.16);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetErrorX(0.5);

def makeCanvas(hists, names, canname, isLog=False):
    
    directory = "plots";
    colors = [2,4,1,6,7];
    
    max = -999.;
    for hist in hists:
        if max < hist.GetMaximum(): max = hist.GetMaximum();
    
    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for i in range(len(names)):
        hists[i].SetLineColor(colors[i]);
        leg.AddEntry(hists[i], names[i], "l");
    
    can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
    hists[0].SetMaximum( 1.2*max );
    hists[0].SetMinimum( 0 );
    hists[0].Draw();
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    if isLog:
        ROOT.gPad.SetLogy();
        hists[0].SetMinimum( 1 );
    can.SaveAs(directory+"/"+canname+".eps");
    can.SaveAs(directory+"/"+canname+".png");
    
    
    for hist in hists:
        hist.Scale(1./hist.Integral());

if __name__ == '__main__':
    
    file = ROOT.TFile("example01.root");
    tree = file.Get("Events");

    h_GEN_pt = ROOT.TH1F("h_GEN_pt",";pT (GeV);N",50,0,500);
    h_GEN_ma = ROOT.TH1F("h_GEN_ma",";mass (GeV);N",50,0,200);
    h_PFCHS_pt = ROOT.TH1F("h_PFCHS_pt",";pT (GeV);N",50,0,500);
    h_PFCHS_ma = ROOT.TH1F("h_PFCHS_ma",";mass (GeV);N",50,0,200);
    h_PUPPI_pt = ROOT.TH1F("h_PUPPI_pt",";pT (GeV);N",50,0,500);
    h_PUPPI_ma = ROOT.TH1F("h_PUPPI_ma",";mass (GeV);N",50,0,200);
    h_PF_pt = ROOT.TH1F("h_PF_pt",";pT (GeV);N",50,0,500);
    h_PF_ma = ROOT.TH1F("h_PF_ma",";mass (GeV);N",50,0,200);

    for i in range(tree.GetEntriesFast()):
        
        tree.GetEntry(i);

        h_GEN_pt.Fill(tree.Jet_GEN_pt[0]);
        h_GEN_ma.Fill(tree.Jet_GEN_ma[0]);

        h_PF_pt.Fill(tree.Jet_PF_pt[0]);
        h_PF_ma.Fill(tree.Jet_PF_ma[0]);

        h_PFCHS_pt.Fill(tree.Jet_PFCHS_pt[0]);
        h_PFCHS_ma.Fill(tree.Jet_PFCHS_ma[0]);

        h_PUPPI_pt.Fill(tree.Jet_PUPPI_pt[0]);
        h_PUPPI_ma.Fill(tree.Jet_PUPPI_ma[0]);

    names = ["Hard","PF","PFCHS","PUPPI"];
    makeCanvas( [h_GEN_pt,h_PF_pt,h_PFCHS_pt,h_PUPPI_pt], names, "pts" )
    makeCanvas( [h_GEN_ma,h_PF_ma,h_PFCHS_ma,h_PUPPI_ma], names, "mass" )





