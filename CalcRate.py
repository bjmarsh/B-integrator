import numpy as np
import ROOT

#rootfile = ROOT.TFile("p_eta_dists/DY_all.root")
rootfile = ROOT.TFile("p_eta_dists/QCD_Pt15_all.root")
#rootfile = ROOT.TFile("p_eta_dists/WJets_all.root")
peta = rootfile.Get("peta")

nbinsx = peta.GetNbinsX()
nbinsy = peta.GetNbinsY()

totint = peta.Integral(0,nbinsx+1,0,nbinsy+1)

plow = 10
etalow = 0.0
etahigh = 0.055
philow = 0.0
phihigh = 0.22

nhits = 10000
eff = 0.0281

## find integral of subregion in p/eta plane
xax = peta.GetXaxis()
yax = peta.GetYaxis()

bx1 = xax.FindBin(plow)
bx2 = nbinsx+1
by1 = yax.FindBin(etalow)
by2 = yax.FindBin(etahigh)

fx1 = (xax.GetBinUpEdge(bx1)-plow)/xax.GetBinWidth(bx1)
fx2 = 1
fy1 = (yax.GetBinUpEdge(by1)-etalow)/yax.GetBinWidth(by1)
fy2 = (etahigh-yax.GetBinLowEdge(by2))/yax.GetBinWidth(by2)
if by1==by2:
    fy2 = (etahigh-etalow)/(yax.GetBinUpEdge(by1)-etalow)

print bx1, fx1, bx2, fx2, by1, fy1, by2, fy2

subint = 0
for i in range(bx1, bx2+1):
    for j in range(by1, by2+1):
        content = peta.GetBinContent(i,j)
        if i==bx1:
            content *= fx1
        if i==bx2:
            content *= fx2
        if j==by1:
            content *= fy1
        if j==by2:
            content *= fy2

        subint += content

print subint, totint

lumi = (1.4e34) * (1e-39)  ## in fb-1/s
hitrate = subint * (phihigh-philow)/(2*np.pi) * eff * lumi  ## in hits/sec


print "Hit rate:", hitrate, "hits/sec"
