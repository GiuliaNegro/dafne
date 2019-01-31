import ROOT
import sys
import tdrstyle,CMS_lumi
from array import array
from math import sqrt
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        print "Usage: %s <data_dir, MC1_dir, MC2_dir, MC3_dir, MC4_dir, MC5_dir, MC6_dir, legendData, legendMC1, legendMC2, legendMC3, legendMC4, legendMC5, legendMC6, isMuMu>" % sys.argv[0]  
        sys.exit(1)

#python Analysis/DYestimation.py DoubleEG/ DYJetsPtBinned/prova4 TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets False -b > DY_EE_estimation

#python Analysis/DYestimation.py SingleMuon/ DYJetsPtBinned/prova4 TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets True -b > DY_MuMu_estimation


dataDir = str(sys.argv[1])
MC1dir = str(sys.argv[2])
MC2dir = str(sys.argv[3])
MC3dir = str(sys.argv[4])
MC4dir = str(sys.argv[5])
MC5dir = str(sys.argv[6])
MC6dir = str(sys.argv[7])
legendData = str(sys.argv[8])
legendMC1 = str(sys.argv[9])
legendMC2 = str(sys.argv[10])
legendMC3 = str(sys.argv[11])
legendMC4 = str(sys.argv[12])
legendMC5 = str(sys.argv[13])
legendMC6 = str(sys.argv[14])
isMuMu = str(sys.argv[15])

inputDir = "Distributions-Moriond17/"+dataDir 
dataInputFile = TFile(inputDir+"distributions_histos.root")
MC1inputFile = TFile(inputDir+MC1dir+"/distributions_histos.root")
MC2inputFile = TFile(inputDir+MC2dir+"/distributions_histos.root")
MC3inputFile = TFile(inputDir+MC3dir+"/distributions_histos.root")
MC4inputFile = TFile(inputDir+MC4dir+"/distributions_histos.root")
MC5inputFile = TFile(inputDir+MC5dir+"/distributions_histos.root")
MC6inputFile = TFile(inputDir+MC6dir+"/distributions_histos.root")
MCFileList = [MC2inputFile, MC3inputFile, MC4inputFile, MC5inputFile, MC6inputFile]

outputDir = "Comparisons-Moriond17/"+dataDir+"DYestimation/prova"
# outputDir = "gnegro@lxplus.cern.ch:/afs/cern.ch/user/g/gnegro/www/cmsWR/Comparisons-Moriond17/"+dataDir+"data-allMCbkg"

# legendList = [legendMC1, legendMC2, legendMC3, legendMC4, legendMC5, legendMC6]
# colorList = [kYellow, kGreen+1, kAzure+7, kMagenta+1, kBlue, kRed-4]#kOrange+10]#kViolet+7] 
# xsecList = [1, 1, 1, 1, 1, 1]
# sumWeightsList = [1, 1, 1, 1, 1, 1]   

lumi = 35.867 
# CMS_lumi.lumi_13TeV = str(lumi)+" fb^{-1}"



fileData = TFile(inputDir+"/Zmass_histos.root")
file1 = TFile(inputDir+MC1dir+"/Zmass_histos.root")
file2 = TFile(inputDir+MC2dir+"/Zmass_histos.root")
file3 = TFile(inputDir+MC3dir+"/Zmass_histos.root")
file4 = TFile(inputDir+MC4dir+"/Zmass_histos.root")
file5 = TFile(inputDir+MC5dir+"/Zmass_histos.root")
file6 = TFile(inputDir+MC6dir+"/Zmass_histos.root")
ZMCFileList = [file2,file3,file4,file5,file6]


histoName = "ZtoEE_mass_histo"
# if isMuMu:
# 	histoName = "ZtoMuMu_mass_histo"


histoData = TH1D()
fileData.GetObject(histoName,histoData)
# histoData.Sumw2()

# c = TCanvas()
# histoData.Draw()
# histoData.Write()
# c.SaveAs("prova.png")

histoDY = TH1D()
file1.GetObject(histoName,histoDY)
# histoDY.Sumw2()
histoDY.Scale(lumi)

# c = TCanvas()
# histoDY.Draw()
# histoDY.Write()
# c.SaveAs("provaDY.png")


histoList = []
for file in ZMCFileList:
	histo = TH1D()
	file.GetObject(histoName,histo)
	# histo.Sumw2()
	histo.Scale(lumi)
	histoList.append(histo)

histoOthers = histoList[0].Clone("histoOthers")
histoOthers.Reset()
for histo in histoList:
	histoOthers.Add(histo)


integralData = histoData.Integral(histoData.FindBin(80), histoData.FindBin(100))
integralDY = histoDY.Integral(histoDY.FindBin(80), histoDY.FindBin(100))
integralOthers = histoOthers.Integral(histoOthers.FindBin(80), histoOthers.FindBin(100))

print "integralData = ", integralData 
print "integralDY = ", integralDY
print "integralOthers = ", integralOthers


# histoDY.Scale(integralData/integralDY)
# integralDYNorm = histoDY.Integral(histoDY.FindBin(80), histoDY.FindBin(100))
# print "integralDYNorm = ", integralDYNorm  # = integralData

# for histo in histoList:
# 	histo.Scale(integralData/integralOthers)
# histoOthersNorm = histoList[0].Clone("histoOthersNorm")
# histoOthersNorm.Reset()
# for histo in histoList:
# 	histoOthersNorm.Add(histo)
# integralOthersNorm = histoOthersNorm.Integral(histoOthers.FindBin(80), histoOthers.FindBin(100))
# print "integralOthersNorm = ", integralOthersNorm  # = integralData


SF = (integralData - integralOthers) / integralDY
eSF = sqrt( (1./integralDY)*(1./integralDY)*integralData + (1./integralDY)*(1./integralDY)*integralOthers + ((integralData - integralOthers)/(integralDY*integralDY))*((integralData - integralOthers)/(integralDY*integralDY))*integralDY)

print "SF = ", SF, ", eSF = ", eSF
print "SF data/DY+others= ", (integralData) / (integralDY + integralOthers)  #circa uguale a sopra




# histoMWRname = "mass_multiLeptonMultiJets_lowMllCR_histo"
# histoPtDiLeptonName = "pt_dileptons_lowMllCR_histo"
# # histoMDiLeptonName = "mass_dileptons_lowMllCR_histo"



# histoData_MWR = TH1D()
# dataInputFile.GetObject(histoMWRname,histoData_MWR)

# histoData_PtDiLepton = TH1D()
# dataInputFile.GetObject(histoPtDiLeptonName,histoData_PtDiLepton)


# histoDY_MWR = TH1D()
# MC1inputFile.GetObject(histoMWRname,histoDY_MWR)

# histoDY_PtDiLepton = TH1D()
# dataInputFile.GetObject(histoPtDiLeptonName,histoDY_PtDiLepton)


# histoList_MWR = []
# histoList_PtDiLepton = []

# for file in MCFileList:
# 	histo_MWR = TH1D()
# 	file.GetObject(histoMWRname,histo_MWR)
# 	histoList_MWR.append(histo_MWR)

# 	histo_PtDiLepton = TH1D()
# 	file.GetObject(histoPtDiLeptonName,histo_PtDiLepton)
# 	histoList_PtDiLepton.append(histo_PtDiLepton)


# histoOthers_MWR = histoList_MWR[0].Clone("histoOthers_MWR")
# histoOthers_MWR.Reset()
# for histo_MWR in histoList_MWR:
# 	histoOthers_MWR.Add(histo_MWR)


# histoOthers_PtDiLepton = histoList_PtDiLepton[0].Clone("histoOthers_PtDiLepton")
# histoOthers_PtDiLepton.Reset()
# for histo_PtDiLepton in histoList_PtDiLepton:
# 	histoOthers_PtDiLepton.Add(histo_PtDiLepton)



# histoOthers_MWR.Scale(-1.)
# histoData_MWR.Add(histoOthers_MWR)
# histoData_MWR.Divide(histoDY_MWR)
# histoData_MWR.Draw()
# histoData_MWR.Write()

# histoOthers_PtDiLepton.Scale(-1.)
# histoData_PtDiLepton.Add(histoOthers_MWR)
# histoData_PtDiLepton.Divide(histoDY_MWR)
# histoData_PtDiLepton.Draw()
# histoData_PtDiLepton.Write()

