import ROOT
import sys
from array import array
from math import sqrt
from ROOT import *
from DrawingAndComparisonFunctions import *


# if len(sys.argv) == 1:
#         print "Usage: %s <dirEE, dirMuMu, dirEMu>" % sys.argv[0]  
#         sys.exit(1)

#python Analysis/TTbarEstimation.py DoubleEG/ SingleMuon/ MuonEG/

# dirEE = str(sys.argv[1])
# dirMuMu = str(sys.argv[2])
# dirEMu = str(sys.argv[3])

# EEinputFile = = TFile(dirEE+"TTJets/distributions_histos.root")
# MuMuInputFile = = TFile(dirMuMu+"TTJets/distributions_histos.root")
# EMuInputFile = = TFile(dirEMu+"TTJets/distributions_histos.root")
# EMuDataInputFile = = TFile(dirEMu+"distributions_histos.root")


#python Analysis/TTbarEstimation.py

EEinputFile = = TFile("DoubleEG/TTJets/distributions_histos.root")
MuMuInputFile = = TFile("SingleMuon/TTJets/distributions_histos.root")
EMuInputFile = = TFile("MuonEG/TTJets/distributions_histos.root")
EMuDataInputFile = = TFile("MuonEG/distributions_histos.root")


outputDir = "Comparisons-Moriond17/TTbarEstimation"
# outputDir = "gnegro@lxplus.cern.ch:/afs/cern.ch/user/g/gnegro/www/cmsWR/Comparisons-Moriond17/"+dataDir+"data-allMCbkg"

# legendList = [legendMC1, legendMC2, legendMC3, legendMC4, legendMC5, legendMC6]
# colorList = [kYellow, kGreen+1, kAzure+7, kMagenta+1, kBlue, kRed-4]#kOrange+10]#kViolet+7] 
# xsecList = [1, 1, 1, 1, 1, 1]
# sumWeightsList = [1, 1, 1, 1, 1, 1]   

lumi = 35.867 
# CMS_lumi.lumi_13TeV = str(lumi)+" fb^{-1}"


EESF = 6.68110e-01  ##da cambiare quando fatto ratio
MuMuSF = 4.78271e-01  ##da cambiare quando fatto ratio

longEESF = str(EESF)
shortEESF = longEESF.substr(0,4)

longMuMuSF = str(MuMuSF)
shortMuMuSF = longMuMuSF.substr(0,4)


histoMWRsignalRegionName = "mass_multiLeptonMultiJets_signalRegion_histo"
histoMWRflavorSidebandName = "mass_multiLeptonMultiJets_eMuSidebandCR_histo"


histoEMuData = TH1D()
EMuDataInputFile.GetObject(histoMWRflavorSidebandName,histoData)

histoEMu = TH1D()
EMuInputFile.GetObject(histoMWRflavorSidebandName,histoDY)

histoEE = TH1D()
EEinputFile.GetObject(histoMWRsignalRegionName,histoEE)

histoMuMu = TH1D()
MuMuInputFile.GetObject(histoMWRsignalRegionName,histoMuMu)


#comparison EMu-EE MC
canvasEE = TCanvas("canvasEE", "", 0, 0, 600, 600)
histoEMu.DrawNormalized()
histoEE.SetLineColor(kRed)
histoEE.DrawNormalized("same histo")
legendEE = TLegend(0.72, 0.50, 0.98, 0.70)
legendEE.AddEntry( histoEMu, "EMu" );
legendEE.AddEntry( histoEE, "EE" );
legendEE.Draw()
canvasEE.SaveAs(outputDir+"/EMu_EE_MC.png")
canvasEE.SaveAs(outputDir+"/EMu_EE_MC.root")

#comparison EMu-MuMu MC
canvasMuMu = TCanvas("canvasMuMu", "", 0, 0, 600, 600)
histoEMu.DrawNormalized()
histoMuMu.SetLineColor(kRed)
histoMuMu.DrawNormalized("same histo")
legendMuMu = TLegend(0.72, 0.50, 0.98, 0.70)
legendMuMu.AddEntry( histoEMu, "EMu" );
legendMuMu.AddEntry( histoMuMu, "MuMu" );
legendMuMu.Draw()
canvasMuMu.SaveAs(outputDir+"/EMu_MuMu_MC.png")
canvasMuMu.SaveAs(outputDir+"/EMu_MuMu_MC.root")


#ratio EE-EMu MC
canvasRatioEE = TCanvas("canvasRatioEE", "", 0, 0, 600, 600)
histoRatioEE = histoEE.Clone()
histoRatioEE.Sumw2()
histoRatioEE.Divide(histoEMu)
# histoRatioEE.GetYaxis().SetRangeUser(0.3,0.6)
histoRatioEE.GetYaxis().SetTitle("ratio M_{EEJJ}/M_{EMuJJ}")

f_EE = TF1("f_EE","[0]",600,2500)   ## controllare estremi plot
f_EE.SetLineColor(kBlue)

chiSquaredBoxEE = TPaveText(1500.,0.54,2000.,0.58)
chiSquaredBoxEE.SetFillColor(kWhite)

histoRatioEE.Fit("f_EE")
histoRatioEE.Draw()
f_EE.Draw("same")

chiSquaredBoxEE.AddText( TString(chiSquaredNdofString(f_EE)) )
chiSquaredBoxEE.AddText( TString("ratio = " + shortEESF) )
chiSquaredBoxEE.Draw("same")

canvasRatioEE.SaveAs(outputDir+"/ratio_EE-EMu.png")
canvasRatioEE.SaveAs(outputDir+"/ratio_EE-EMu.root")

canvasRatioEE.SetLogx(1)
chiSquaredBoxEE.DrawPave(500.,0.54,1000.,0.58,4,"same")
canvasRatioEE.SaveAs(outputDir+"/ratio_EE-EMu_logx.png")
canvasRatioEE.SaveAs(outputDir+"/ratio_EE-EMu_logx.root")


#ratio MuMu-EMu MC
canvasRatioMuMu = TCanvas("canvasRatioMuMu", "", 0, 0, 600, 600)
histoRatioMuMu = histoMuMu.Clone()
histoRatioMuMu.Sumw2()
histoRatioMuMu.Divide(histoEMu)
# histoRatioMuMu.GetYaxis().SetRangeUser(0.4,0.8)
histoRatioMuMu.GetYaxis().SetTitle("ratio M_{MuMuJJ}/M_{EMuJJ}")

f_MuMu = TF1("f_MuMu","[0]",600,2500)   ## controllare estremi plot
f_MuMu.SetLineColor(kBlue)

chiSquaredBoxMuMu = TPaveText(1500.,0.73,2000.,0.79)
chiSquaredBoxMuMu.SetFillColor(kWhite)

histoRatioMuMu.Fit("f_MuMu")
histoRatioMuMu.Draw()
f_MuMu.Draw("same")

chiSquaredBoxMuMu.AddText( TString(chiSquaredNdofString(f_MuMu)) )
chiSquaredBoxMuMu.AddText( TString("ratio = " + shortMuMuSF) )
chiSquaredBoxMuMu.Draw("same")

canvasRatioMuMu.SaveAs(outputDir+"/ratio_MuMu-EMu.png")
canvasRatioMuMu.SaveAs(outputDir+"/ratio_MuMu-EMu.root")

canvasRatioMuMu.SetLogx(1)
chiSquaredBoxMuMu.DrawPave(500.,0.73,1000.,0.79,4,"same")
canvasRatioMuMu.SaveAs(outputDir+"/ratio_MuMu-EMu_logx.png")
canvasRatioMuMu.SaveAs(outputDir+"/ratio_MuMu-EMu_logx.root")








# gStyle->SetOptStat("nemr");
canvasEEEMu = TCanvas("canvasEEEMu","",600,600)
canvEEEMu.cd()
legEEEMu = TLegend(0.72,0.6,0.98,0.8)
legEEEMu.AddEntry(histoEMu,"EMu")
legEEEMu.AddEntry(histoEE,"EE")
histoEMu.Draw("histo")
histoEE.Draw("Psame")
legEEEMu.Draw()
canvasEEEMu.SaveAs(outputDir+"/EMu_EE_signalRegion.png","recreate")
canvasEEEMu.SaveAs(outputDir+"/EMu_EE_signalRegion.png","recreate")

canvasEEEMuData = TCanvas("canvasEEEMuData","",600,600)
canvasEEEMuData.cd()
legEEEMuData = TLegend(0.72,0.6,0.98,0.8)
legEEEMuData.AddEntry(histoEMuData,"Rescaled EMu Data")
legEEEMuData.AddEntry(histoEE,"EE MC")
histoEE.Draw("histo")
histoEMuData.Scale(EESF)
histoEMuData.SetMarkerStyle(2)
histoEMuData.SetMarkerSize(2)
histoEMuData.Draw("Psame")
legEEEMuData.Draw()
canvasEEEMuData.SaveAs(outputDir+"/rescaled_EMuData_EEMC_signalRegion.png","recreate")
canvasEEEMuData.SaveAs(outputDir+"/rescaled_EMuData_EEMC_signalRegion.root","recreate")


canvasMuMuEMu = TCanvas("canvasMuMuEMu","",600,600)
canvMuMuEMu.cd()
legMuMuEMu = TLegend(0.72,0.6,0.98,0.8)
legMuMuEMu.AddEntry(histoEMu,"EMu")
legMuMuEMu.AddEntry(histoMuMu,"MuMu")
histoEMu.Draw("histo")
histoMuMu.Draw("Psame")
legMuMuEMu.Draw()
canvasMuMuEMu.SaveAs(outputDir+"/EMu_MuMu_signalRegion.png","recreate")
canvasMuMuEMu.SaveAs(outputDir+"/EMu_MuMu_signalRegion.png","recreate")

canvasMuMuEMuData = TCanvas("canvasMuMuEMuData","",600,600)
canvasMuMuEMuData.cd()
legMuMuEMuData = TLegend(0.72,0.6,0.98,0.8)
legMuMuEMuData.AddEntry(histoEMuData,"Rescaled EMu Data")
legMuMuEMuData.AddEntry(histoMuMu,"MuMu MC")
histoMuMu.Draw("histo")
histoEMuData.Scale(1/EESF)  #undo the previous scaling
histoEMuData.Scale(MuMuSF)
histoEMuData.SetMarkerStyle(2)
histoEMuData.SetMarkerSize(2)
histoEMuData.Draw("Psame")
legMuMuEMuData.Draw()
canvasMuMuEMuData.SaveAs(outputDir+"/rescaled_EMuData_MuMuMC_signalRegion.png","recreate")
canvasMuMuEMuData.SaveAs(outputDir+"/rescaled_EMuData_MuMuMC_signalRegion.root","recreate")





def chiSquaredNdofString(fit):
	tempchiSqd = "#chi^{2}  =  "
	chiSqdVal = str(fit.GetChisquare())
	ndof = str(fit.GetNDF())
	tempchiSqd += chiSqdVal.substr(0,4)
	tempchiSqd += " / "
	tempchiSqd += ndof.substr(0,2)
	return tempchiSqd





histoMWRname = "mass_multiLeptonMultiJets_lowMllCR_histo"
histoPtDiLeptonName = "pt_dileptons_lowMllCR_histo"
# histoMDiLeptonName = "mass_dileptons_lowMllCR_histo"



histoData_MWR = TH1D()
dataInputFile.GetObject(histoMWRname,histoData_MWR)

histoData_PtDiLepton = TH1D()
dataInputFile.GetObject(histoPtDiLeptonName,histoData_PtDiLepton)


histoDY_MWR = TH1D()
MC1inputFile.GetObject(histoMWRname,histoDY_MWR)

histoDY_PtDiLepton = TH1D()
dataInputFile.GetObject(histoPtDiLeptonName,histoDY_PtDiLepton)


histoList_MWR = []
histoList_PtDiLepton = []

for file in MCFileList:
	histo_MWR = TH1D()
	file.GetObject(histoMWRname,histo_MWR)
	histoList_MWR.append(histo_MWR)

	histo_PtDiLepton = TH1D()
	file.GetObject(histoPtDiLeptonName,histo_PtDiLepton)
	histoList_PtDiLepton.append(histo_PtDiLepton)


histoOthers_MWR = histoList_MWR[0].Clone("histoOthers_MWR")
histoOthers_MWR.Reset()
for histo_MWR in histoList_MWR:
	histoOthers_MWR.Add(histo_MWR)


histoOthers_PtDiLepton = histoList_PtDiLepton[0].Clone("histoOthers_PtDiLepton")
histoOthers_PtDiLepton.Reset()
for histo_PtDiLepton in histoList_PtDiLepton:
	histoOthers_PtDiLepton.Add(histo_PtDiLepton)



histoOthers_MWR.Scale(-1.)
histoData_MWR.Add(histoOthers_MWR)
histoData_MWR.Divide(histoDY_MWR)
histoData_MWR.Draw()
histoData_MWR.Write()

histoOthers_PtDiLepton.Scale(-1.)
histoData_PtDiLepton.Add(histoOthers_MWR)
histoData_PtDiLepton.Divide(histoDY_MWR)
histoData_PtDiLepton.Draw()
histoData_PtDiLepton.Write()

