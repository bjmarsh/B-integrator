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

pt_DY = ROOT.TH1F("pt_DY", "Pt", 50, 0, 100)
pt_QCD = ROOT.TH1F("pt_QCD", "Pt", 50, 0, 100)
pt_WJets = ROOT.TH1F("pt_WJets", "Pt", 50, 0, 100)

pti_DY = ROOT.TH1F("pti_DY", "Pt", 50, 0, 100)
pti_QCD = ROOT.TH1F("pti_QCD", "Pt", 50, 0, 100)
pti_WJets = ROOT.TH1F("pti_WJets", "Pt", 50, 0, 100)

for i in range(dataDY.shape[0]):
    pt_DY.Fill(dataDY[i,2])
    pti_DY.Fill(dataDY[i,10]/1000)
for i in range(dataQCD.shape[0]):
    pt_QCD.Fill(dataQCD[i,2])
    pti_QCD.Fill(dataQCD[i,10]/1000)
for i in range(dataWJets.shape[0]):
    pt_WJets.Fill(dataWJets[i,2])
    pti_WJets.Fill(dataWJets[i,10]/1000)

pt_DY.Scale(1.0/pt_DY.Integral() * 0.00089 * 3600)
pt_QCD.Scale(1.0/pt_QCD.Integral() * 0.099 * 3600)
pt_WJets.Scale(1.0/pt_WJets.Integral() * 0.0024 * 3600)
pti_DY.Scale(1.0/pti_DY.Integral() * 0.00089 * 3600)
pti_QCD.Scale(1.0/pti_QCD.Integral() * 0.099 * 3600)
pti_WJets.Scale(1.0/pti_WJets.Integral() * 0.0024 * 3600)

pti_all = ROOT.TH1F()
pti_DY.Copy(pti_all)
pti_all.Add(pti_WJets)
pti_all.Add(pti_QCD)

pt_DY.SetFillColor(ROOT.kRed-7)
pt_WJets.SetFillColor(ROOT.kSpring-5)
pt_QCD.SetFillColor(ROOT.kAzure+7)

pti_all.SetLineColor(ROOT.kBlack)
pti_all.SetLineWidth(2)

s1 = ROOT.THStack("pt_stack","P_{T} (produced) of incident muons")
s1.Add(pt_DY)
s1.Add(pt_WJets)
s1.Add(pt_QCD)

c1 = ROOT.TCanvas()
c1.SetLogy()

s1.Draw()
pti_all.Draw("SAME")

s1.GetYaxis().SetTitle("Hits per hour / 2 GeV")
s1.GetXaxis().SetTitle("P_{T} (GeV)")

leg = ROOT.TLegend(0.60,0.70,0.88,0.88);
leg.AddEntry(pt_QCD,"QCD","f")
leg.AddEntry(pt_WJets,"WJets","f")
leg.AddEntry(pt_DY,"DY","f")
leg.AddEntry(pti_all,"P_{T} upon hit","l")
leg.Draw()

raw_input()
