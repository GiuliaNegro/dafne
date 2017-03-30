import FWCore.ParameterSet.Config as cms
import json
import math
import os
from ROOT import TH2F, TFile


#Run2016B (5.929/fb)
#Run2016C (2.646/fb)
#Run2016D (4.353/fb)
#Run2016E (4.117/fb)
#Run2016F (3.186/fb)
#Run2016G (7.721/fb)
#Run2016H_v2 (8.636/fb)
#Run2016H_v3 (0.221/fb)

w1 = (5.929+2.646+4.353+4.117+3.186)/(5.929+2.646+4.353+4.117+3.186+7.721+8.636+0.221)
w2 = 1.0 - w1


TRSF1 = TFile(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_TrigBF.root'))
TRSF2 = TFile(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_TrigGH.root'))

TRsf1 = TH2F()
TRSF1.GetObject("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio", TRsf1)
TReffData1 = TH2F()
TRSF1.GetObject("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/abseta_pt_DATA", TReffData1)
TReffMC1 = TH2F()
TRSF1.GetObject("Mu50_OR_TkMu50_PtEtaBins/efficienciesMC/abseta_pt_MC", TReffMC1)

TRsf2 = TH2F()
TRSF2.GetObject("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio", TRsf2)
TReffData2 = TH2F()
TRSF2.GetObject("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/abseta_pt_DATA", TReffData2)
TReffMC2 = TH2F()
TRSF2.GetObject("Mu50_OR_TkMu50_PtEtaBins/efficienciesMC/abseta_pt_MC", TReffMC2)



print "SF1: ", TRsf1.GetNbinsX()+1, " xbins, ", TRsf1.GetNbinsY()+1, " ybins"
# for x in range(1, TRsf1.GetNbinsX()+1):
# 	for y in range(1, TRsf1.GetNbinsY()+1):
# 		sf1 = TRsf1.GetBinContent(x,y)
# 		e1 = TRsf1.GetBinError(x,y)
		# print "binX ", x, ", binY ", y, ", sf1 = ", sf1, ", e1 = ", e1

print "SF2: ", TRsf2.GetNbinsX()+1, " xbins, ", TRsf2.GetNbinsY()+1, " ybins"
# for x in range(1, TRsf2.GetNbinsX()+1):
# 	for y in range(1, TRsf2.GetNbinsY()+1):
# 		sf2 = TRsf2.GetBinContent(x,y)
# 		e2 = TRsf2.GetBinError(x,y)
		# print "binX ", x, ", binY ", y, ", sf2 = ", sf2, ", e2 = ", e2

if (TRsf1.GetNbinsX() != TRsf2.GetNbinsX()) or (TRsf1.GetNbinsY() != TRsf2.GetNbinsY()):
	print "different nbins"

print "SF: ", TRsf1.GetNbinsX()+1, " xbins, ", TRsf1.GetNbinsY()+1, " ybins"
for x in range(1, TRsf1.GetNbinsX()+1):
	for y in range(1, TRsf1.GetNbinsY()+1):
		sf1 = TRsf1.GetBinContent(x,y)
		e1 = TRsf1.GetBinError(x,y)
		sf2 = TRsf2.GetBinContent(x,y)
		e2 = TRsf2.GetBinError(x,y)

		sf = sf1*w1 + sf2*w2
		err = math.sqrt(pow(e1*w1,2) + pow(e2*w2,2)) 
		print "binX ", x, ", binY ", y, ", sf = ", sf, ", err = ", err



print "effData1: ", TReffData1.GetNbinsX()+1, " xbins, ", TReffData1.GetNbinsY()+1, " ybins"
# for x in range(1, TReffData1.GetNbinsX()+1):
# 	for y in range(1, TReffData1.GetNbinsY()+1):
# 		effData1 = TReffData1.GetBinContent(x,y)
# 		eData1 = TReffData1.GetBinError(x,y)
# 		print "binX ", x, ", binY ", y, ", effData1 = ", effData1, ", eData1 = ", eData1

print "effData2: ", TReffData2.GetNbinsX()+1, " xbins, ", TReffData2.GetNbinsY()+1, " ybins"
# for x in range(1, TReffData2.GetNbinsX()+1):
# 	for y in range(1, TReffData2.GetNbinsY()+1):
# 		effData2 = TReffData2.GetBinContent(x,y)
# 		eData2 = TReffData2.GetBinError(x,y)
# 		print "binX ", x, ", binY ", y, ", effData2 = ", effData2, ", eData2 = ", eData2

if (TReffData1.GetNbinsX() != TReffData2.GetNbinsX()) or (TReffData1.GetNbinsY() != TReffData2.GetNbinsY()):
	print "different nbins"

print "effData: ", TReffData1.GetNbinsX()+1, " xbins, ", TReffData1.GetNbinsY()+1, " ybins"
for x in range(1, TReffData1.GetNbinsX()+1):
	for y in range(1, TReffData1.GetNbinsY()+1):
		effData1 = TReffData1.GetBinContent(x,y)
		eData1 = TReffData1.GetBinError(x,y)
		effData2 = TReffData2.GetBinContent(x,y)
		eData2 = TReffData2.GetBinError(x,y)

		effData = effData1*w1 + effData2*w2
		errData = math.sqrt(pow(eData1*w1,2) + pow(eData2*w2,2)) 
		print "binX ", x, ", binY ", y, ", effData = ", effData, ", errData = ", errData




# print "effMC: ", TReffMC1.GetNbinsX()+1, " xbins, ", TReffMC1.GetNbinsY()+1, " ybins"
# for x in range(1, TReffMC1.GetNbinsX()+1):
# 	for y in range(1, TReffMC1.GetNbinsY()+1):
# 		effMC1 = TReffMC1.GetBinContent(x,y)
# 		eMC1 = TReffMC1.GetBinError(x,y)
# 		print "binX ", x, ", binY ", y, ", effMC1 = ", effMC1, ", eMC1 = ", eMC1


# print "ratio effMC/effData: ", TReffMC1.GetNbinsX()+1, " xbins, ", TReffMC1.GetNbinsY()+1, " ybins"
# for x in range(1, TReffMC1.GetNbinsX()+1):
# 	for y in range(1, TReffMC1.GetNbinsY()+1):
# 		effData1 = TReffData1.GetBinContent(x,y)
# 		effMC1 = TReffMC1.GetBinContent(x,y)
# 		print "binX ", x, ", binY ", y, ", ratio effMC/effData = ", effMC1/effData1


# print "ratio effData/effMC: ", TReffMC1.GetNbinsX()+1, " xbins, ", TReffMC1.GetNbinsY()+1, " ybins"
# for x in range(1, TReffMC1.GetNbinsX()+1):
# 	for y in range(1, TReffMC1.GetNbinsY()+1):
# 		effData1 = TReffData1.GetBinContent(x,y)
# 		effMC1 = TReffMC1.GetBinContent(x,y)
# 		print "binX ", x, ", binY ", y, ", ratio effData/effMC = ", effData1/effMC1

