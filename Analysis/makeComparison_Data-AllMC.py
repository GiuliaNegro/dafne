import ROOT
import sys
from array import array
from math import sqrt
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        print "Usage: %s <data_dir, MC1_dir, MC2_dir, MC3_dir, MC4_dir, MC5_dir, MC6_dir, legendData, legendMC1, legendMC2, legendMC3, legendMC4, legendMC5, legendMC6, norm, log, suff>" % sys.argv[0]  
        sys.exit(1)

#python Analysis/makeComparison_Data-AllMC.py DoubleEG/ DYJetsPtBinned/ptBinnedScaled TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets False True "" -b > DoubleEG
#python Analysis/makeComparison_Data-AllMC.py DoubleEG/ DYJetsPtBinned/ptBinnedScaled TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets True True "_normArea" -b > DoubleEG_normArea

#python Analysis/makeComparison_Data-AllMC.py SingleMuon/ DYJetsPtBinned/ptBinnedScaled TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets False True "" -b > SingleMuon
#python Analysis/makeComparison_Data-AllMC.py SingleMuon/ DYJetsPtBinned/ptBinnedScaled TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets True True "_normArea" -b > SingleMuon_normArea

#python Analysis/makeComparison_Data-AllMC.py MuonEG/ DYJetsPtBinned/prova4 TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets False True "" -b > MuonEG
#python Analysis/makeComparison_Data-AllMC.py MuonEG/ DYJetsPtBinned/prova4 TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets True True "_normArea" -b > MuonEG_normArea


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
norm = str(sys.argv[15])
log = str(sys.argv[16])
suff = str(sys.argv[17])

inputDir = "Distributions-Moriond17/"+dataDir 
dataInputFile = TFile(inputDir+"distributions_histos.root")
MC1inputFile = TFile(inputDir+MC1dir+"/distributions_histos.root")
MC2inputFile = TFile(inputDir+MC2dir+"/distributions_histos.root")
MC3inputFile = TFile(inputDir+MC3dir+"/distributions_histos.root")
MC4inputFile = TFile(inputDir+MC4dir+"/distributions_histos.root")
MC5inputFile = TFile(inputDir+MC5dir+"/distributions_histos.root")
MC6inputFile = TFile(inputDir+MC6dir+"/distributions_histos.root")
MCFileList = [MC1inputFile, MC2inputFile, MC3inputFile, MC4inputFile, MC5inputFile, MC6inputFile]

outputDir = "Comparisons-Moriond17/"+dataDir+"/DYptBinnedScaled/data-allMCbkg"
# outputDir = "gnegro@lxplus.cern.ch:/afs/cern.ch/user/g/gnegro/www/cmsWR/Comparisons-Moriond17/"+dataDir+"data-allMCbkg"

legendList = [legendMC1, legendMC2, legendMC3, legendMC4, legendMC5, legendMC6]
colorList = [kYellow, kGreen+1, kAzure+7, kMagenta+1, kBlue, kRed-4]#kOrange+10]#kViolet+7] 
xsecList = [1, 1, 1, 1, 1, 1]
sumWeightsList = [1, 1, 1, 1, 1, 1]   

lumi = 35.867 
# CMS_lumi.lumi_13TeV = str(lumi)+" fb^{-1}"
CMS_lumi.lumi_13TeV = "35.9 fb^{-1}"

# (histoName, etaRegion, dataInputFile, MCFileList, output_name, suff, zoomX, xRangeMin, xRangeMax, zoomY, yRangeMin, yRangeMax, leg1_name, legendList, 
	# leftLegends, log, doRebin, doRestrictedIntegral, normArea, colorList, xsecList, sumWeightsList, lumi, title="", xTitle="", yTitle=""):






etaMassNames = [""]#, "_EB-EB", "_EE-EE"]#, "_EB-EE", "_noEB-EB"]
regionNames = ["_lowMllCR"]#,"_lowMlljjCR", "_eMuSidebandCR", "_signalRegion"]
objNames = ["leadingEle", "subLeadingEle", "leadingJet", "subLeadingJet"]#, "leadingMu", "subLeadingMu","leadingLepton", "subLeadingLepton"]
if dataDir == "SingleMuon/":
	objNames = ["leadingMu", "subLeadingMu", "leadingJet", "subLeadingJet"]
if dataDir == "MuonEG/":
	regionNames = ["_eMuSidebandCR"]
	objNames = ["leadingEle", "leadingMu", "leadingJet", "subLeadingJet"]

for regionName in regionNames:

	for objName in objNames:
		if objName == "leadingEle" or objName== "subLeadingEle" or objName == "leadingMu" or objName== "subLeadingMu":
			plotDataAndMCHistosAndRatio("pt_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 300, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T} "+objName, yTitle="# Events")
		else:
			plotDataAndMCHistosAndRatio("pt_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 700, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T} "+objName, yTitle="# Events")

		plotDataAndMCHistosAndRatio("eta_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#eta "+objName, yTitle="# Events")
		plotDataAndMCHistosAndRatio("phi_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#phi "+objName, yTitle="# Events")


	for etaMassName in etaMassNames:
		# if "lowMlljjCR" in regionName:
		# 	plotDataAndMCHistosAndRatio("mass_multiLeptonMultiJets"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 1000, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{lljj}", yTitle="# Events")
		# else:
		plotDataAndMCHistosAndRatio("mass_multiLeptonMultiJets"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 3000, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{lljj}", yTitle="# Events")

		plotDataAndMCHistosAndRatio("mass_dileptons"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 300, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{ll}", yTitle="# Events")
		plotDataAndMCHistosAndRatio("mass_dijets"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 1000, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{ll}", yTitle="# Events")

	
	plotDataAndMCHistosAndRatio("pt_dileptons_"+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 700, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T}(ll)", yTitle="# Events")




# objNames = ["ele", "goodEle", "jet", "goodJet"]#, "muon", "goodMuon", "muonRochCor", "goodMuonRochCor"]
# for objName in objNames:
# 	plotDataAndMCHistosAndRatio("pt_"+objName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 700, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T} "+objName, yTitle="# Events")
# 	plotDataAndMCHistosAndRatio("eta_"+objName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#eta "+objName, yTitle="# Events")
# 	plotDataAndMCHistosAndRatio("phi_"+objName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#phi "+objName, yTitle="# Events")


plotDataAndMCHistosAndRatio("nvtx_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 50, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="nvtx", yTitle="# Events")
plotDataAndMCHistosAndRatio("rho_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 50, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="rho", yTitle="# Events")


# vtxHistoNames = ["vtx_x", "vtx_y", "vtx_z"]
# for histoName in vtxHistoNames:
	# plotDataAndMCHistosAndRatio(histoName+"_histo", "", dataInputFile, MCFileList, outputDir, "", False, 0, 0, False, 0, 0, legendData, legendList, True, False, False, False, norm, colorList)




if dataDir == "DoubleEG/" or dataDir == "SingleMuon/":
	fileData = TFile(inputDir+"/Zmass_histos.root")
	file1 = TFile(inputDir+MC1dir+"/Zmass_histos.root")
	file2 = TFile(inputDir+MC2dir+"/Zmass_histos.root")
	file3 = TFile(inputDir+MC3dir+"/Zmass_histos.root")
	file4 = TFile(inputDir+MC4dir+"/Zmass_histos.root")
	file5 = TFile(inputDir+MC5dir+"/Zmass_histos.root")
	file6 = TFile(inputDir+MC6dir+"/Zmass_histos.root")
	ZMCFileList = [file1,file2,file3,file4,file5,file6]

	etaMassNames = ["", "_EB-EB", "_EE-EE", "_EB-EE", "_noEB-EB"]
	Znames = ["toEE"]#_MLMJelectrons"]#,"toEE_electrons"]
	if dataDir == "SingleMuon/":
		Znames = ["toMuMu"]#_MLMJmuons"]#,"toMuMu_muons"]

	for etaMassName in etaMassNames:
		for Zname in Znames:
			plotDataAndMCHistosAndRatio("Z"+Zname+"_mass", etaMassName+"_histo", fileData, ZMCFileList, outputDir, suff, True, 70, 110, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{ll}", yTitle="# Events")	





# etaNames = ["", "_EB", "_EE"]
# eleNames = ["ele", "goodEle", "leadingEle", "subLeadingEle", "goodLeadingEle", "goodSubLeadingEle"]
# for eleName in eleNames:
# 	for etaName in etaNames:

# 		plotDataAndMCHistosAndRatio("etaSC_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff+"/HEEP", False, 0, 0, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="etaSC_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("isEcalDriven_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff+"/HEEP", False, 0, 0, False, 0, 0, legendData, legendList, True, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="isEcalDriven_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("dEtaIn_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, -0.1, 0.1, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="dEtaIn_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("dPhiIn_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, -0.1, 0.1, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="dPhiIn_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("hOverE_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 0.3, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="hOverE_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("sigmaIetaIeta_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="sigmaIetaIeta_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("e2x5_e5x5_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, True, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="e2x5_e5x5_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("e1x5_e5x5_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, True, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="e1x5_e5x5_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("EmHadDepth1Iso_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 20, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="EmHadDepth1Iso_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("missingHits_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="missingHits_"+eleName+etaName, yTitle="# Events")

# 		plotDataAndMCHistosAndRatio("dxy_"+eleName, etaName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 0.05, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="dxy_"+eleName+etaName, yTitle="# Events")



