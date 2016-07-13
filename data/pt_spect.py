#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import ROOT

fdy = "DY.txt"
fqcd = "QCD.txt"
fwjets = "WJets.txt"

with open(fdy) as fid:
    dataDY = np.loadtxt(fid)
with open(fqcd) as fid:
    dataQCD = np.loadtxt(fid)
with open(fwjets) as fid:
    dataWJets = np.loadtxt(fid)

pt_DY = ROOT.TH1F("pt_DY", "Pt", 40, 0, 100)
pt_QCD = ROOT.TH1F("pt_QCD", "Pt", 40, 0, 100)
pt_WJets = ROOT.TH1F("pt_WJets", "Pt", 40, 0, 100)

for i in range(dataDY.shape[0]):
    pt_DY.Fill(dataDY[i,2])
for i in range(dataQCD.shape[0]):
    pt_QCD.Fill(dataQCD[i,2])
for i in range(dataWJets.shape[0]):
    pt_WJets.Fill(dataWJets[i,2])

pt_DY.Scale(1.0/pt_DY.Integral() * 0.00089 * 3600)
pt_QCD.Scale(1.0/pt_QCD.Integral() * 0.099 * 3600)
pt_WJets.Scale(1.0/pt_WJets.Integral() * 0.0024 * 3600)

pt_all = ROOT.TH1F()
pt_DY.Copy(pt_all)
pt_all.Add(pt_QCD)
pt_all.Add(pt_WJets)

pt_all.SetLineColor(ROOT.kBlack)
pt_all.SetLineWidth(2)

pt_all.Scale(1.0/pt_all.Integral())

fid_DY = ROOT.TFile("../p_eta_dists/DY_all.root")
fid_QCD = ROOT.TFile("../p_eta_dists/QCD_Pt15_all.root")
fid_WJets = ROOT.TFile("../p_eta_dists/WJets_all.root")
peta_DY = fid_DY.Get("peta")
peta_QCD = fid_QCD.Get("peta")
peta_WJets = fid_WJets.Get("peta")

pta_DY = ROOT.TH1F("pta_DY", "Pt", 40, 0, 100)
pta_QCD = ROOT.TH1F("pta_QCD", "Pt", 40, 0, 100)
pta_WJets = ROOT.TH1F("pta_WJets", "Pt", 40, 0, 100)

for i in range(1,200):
    pta_DY.SetBinContent(i, peta_DY.GetBinContent(i,1))
    pta_QCD.SetBinContent(i, peta_QCD.GetBinContent(i,1))
    pta_WJets.SetBinContent(i, peta_WJets.GetBinContent(i,1))

pta_all = ROOT.TH1F()
pta_DY.Copy(pta_all)
pta_all.Add(pta_QCD)
pta_all.Add(pta_WJets)

pta_all.Scale(1.0/pta_all.Integral())

pta_all.SetLineColor(ROOT.kRed)
pta_all.SetLineWidth(2)

pt_all.Divide(pta_all)

c1 = ROOT.TCanvas()
c1.SetLogy()

pt_all.Draw()
#pt_all.Draw("SAME")

raw_input()
