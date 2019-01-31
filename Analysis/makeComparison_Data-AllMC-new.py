import ROOT
import sys
from array import array
from math import sqrt
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        print "Usage: %s <data_dir, MC1_dir, MC2_dir, MC3_dir, MC4_dir, MC5_dir, MC6_dir, legendData, legendMC1, legendMC2, legendMC3, legendMC4, legendMC5, legendMC6, norm, log, suff>" % sys.argv[0]  
        sys.exit(1)

#python Analysis/makeComparison_Data-AllMC-new.py DoubleEG DYJetsPtBinned TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets False True "" -b > DoubleEG-new
#python Analysis/makeComparison_Data-AllMC-new.py DoubleEG DYJetsPtBinned TTJets WZ ZZ WW WJets data DYJets TTJets WZ ZZ WW WJets True True "_normArea" -b > DoubleEG_normArea-new



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

# inputdir = "Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger/"
# inputdir = "Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger_L-SL_swapped/"
# inputdir = "Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger/rightPt_L-SL/"
inputdir = "Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger/rightPt_L-SL/plotsForEleZ/"

outDir = "Comparisons-Moriond17/MLMJwithSingleObjects/selectedTrigger/rightPt_L-SL/plotsForEleZ/"+dataDir+"/data-allMCbkg"+suff+"/"
suff = ""
legendList = [legendMC1, legendMC2, legendMC3, legendMC4, legendMC5, legendMC6]
colorList = [kYellow, kGreen+1, kAzure+7, kMagenta+1, kBlue, kRed-4]#kOrange+10]#kViolet+7] 
xsecList = [1, 1, 1, 1, 1, 1]
sumWeightsList = [1, 1, 1, 1, 1, 1]   

lumi = 35.867 
# CMS_lumi.lumi_13TeV = str(lumi)+" fb^{-1}"
CMS_lumi.lumi_13TeV = "35.9 fb^{-1}"



etaMassNames = [""]#, "_EB-EB", "_EE-EE"]#, "_EB-EE", "_noEB-EB"]

objNames = ["leadingEle", "subLeadingEle", "leadingJet", "subLeadingJet"]#,"leadingLepton", "subLeadingLepton"]
if dataDir == "SingleMuon/":
	objNames = ["leadingMu", "subLeadingMu", "leadingJet", "subLeadingJet"]
if dataDir == "MuonEG/":
	regionNames = ["_eMuSidebandCR"]
	objNames = ["leadingEle", "leadingMu", "leadingJet", "subLeadingJet"]

triggerDirectories = ["TnPEE/"]#, "TnPMuMu/", "signalEE/", "signalMuMu/", "eMuSideband/"]
triggerNames = ["_TnPEE"]#, "_TnPMuMu", "_signalEE", "_signalMuMu", "_eMuSideband"]

# triggerDirectories = ["signalEE/"]
# triggerNames = ["_signalEE"]

for triggerDirectory, triggerName in zip(triggerDirectories,triggerNames):
	inputDir = inputdir+triggerDirectory

	regionNames = ["_lowMllCR"]#,"_lowMlljjCR", "_signalRegion"]

	if triggerDirectory == "TnPEE/" or triggerDirectory == "TnPMuMu/": 
		regionNames = ["_TnPCR"]

	if triggerDirectory == "eMuSideband/":
		regionNames = ["_eMuSidebandCR"]


	for regionName in regionNames:
		dataInputFile = TFile(inputDir+dataDir+"/distributions_histos"+regionName+".root")
		MC1inputFile = TFile(inputDir+MC1dir+"/distributions_histos"+regionName+".root")
		MC2inputFile = TFile(inputDir+MC2dir+"/distributions_histos"+regionName+".root")
		MC3inputFile = TFile(inputDir+MC3dir+"/distributions_histos"+regionName+".root")
		MC4inputFile = TFile(inputDir+MC4dir+"/distributions_histos"+regionName+".root")
		MC5inputFile = TFile(inputDir+MC5dir+"/distributions_histos"+regionName+".root")
		MC6inputFile = TFile(inputDir+MC6dir+"/distributions_histos"+regionName+".root")
		MCFileList = [MC1inputFile, MC2inputFile, MC3inputFile, MC4inputFile, MC5inputFile, MC6inputFile]

		outputDir = outDir+regionName

# (histoName, etaRegion, dataInputFile, MCFileList, output_name, suff, zoomX, xRangeMin, xRangeMax, zoomY, yRangeMin, yRangeMax, leg1_name, legendList, 
	# leftLegends, log, doRebin, doRestrictedIntegral, normArea, colorList, xsecList, sumWeightsList, lumi, title="", xTitle="", yTitle=""):

		for objName in objNames:
			# if objName == "leadingEle" or objName== "subLeadingEle" or objName == "leadingMu" or objName== "subLeadingMu":
			plotDataAndMCHistosAndRatio("pt_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 300, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T} "+objName, yTitle="# Events")
			# else:
				# plotDataAndMCHistosAndRatio("pt_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 700, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T} "+objName, yTitle="# Events")

			plotDataAndMCHistosAndRatio("eta_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#eta "+objName, yTitle="# Events")
			plotDataAndMCHistosAndRatio("phi_"+objName+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#phi "+objName, yTitle="# Events")


		for etaMassName in etaMassNames:
			# if "lowMlljjCR" in regionName:
			# 	plotDataAndMCHistosAndRatio("mass_multiLeptonMultiJets"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 1000, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{lljj}", yTitle="# Events")
			# else:
			# plotDataAndMCHistosAndRatio("mass_multiLeptonMultiJets"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 3000, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{lljj}", yTitle="# Events")
			plotDataAndMCHistosAndRatio("mass_multiLeptonMultiJets"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 150, 7000, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{lljj}", yTitle="# Events")

			plotDataAndMCHistosAndRatio("mass_dileptons"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 300, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{ll}", yTitle="# Events")
			plotDataAndMCHistosAndRatio("mass_dijets"+regionName, etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 0, 1000, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{ll}", yTitle="# Events")

		plotDataAndMCHistosAndRatio("pt_dileptons"+regionName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 700, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T}(ll)", yTitle="# Events")



	for triggerName in triggerNames:
		dataInputFile = TFile(inputDir+dataDir+"/distributions_histos"+triggerName+".root")
		MC1inputFile = TFile(inputDir+MC1dir+"/distributions_histos"+triggerName+".root")
		MC2inputFile = TFile(inputDir+MC2dir+"/distributions_histos"+triggerName+".root")
		MC3inputFile = TFile(inputDir+MC3dir+"/distributions_histos"+triggerName+".root")
		MC4inputFile = TFile(inputDir+MC4dir+"/distributions_histos"+triggerName+".root")
		MC5inputFile = TFile(inputDir+MC5dir+"/distributions_histos"+triggerName+".root")
		MC6inputFile = TFile(inputDir+MC6dir+"/distributions_histos"+triggerName+".root")
		MCFileList = [MC1inputFile, MC2inputFile, MC3inputFile, MC4inputFile, MC5inputFile, MC6inputFile]

		outputDir = outDir+triggerName


		plotDataAndMCHistosAndRatio("nvtx_histo"+triggerName, "", dataInputFile, MCFileList, outputDir, suff, True, 0, 50, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="nvtx", yTitle="# Events")
		plotDataAndMCHistosAndRatio("rho_histo"+triggerName, "", dataInputFile, MCFileList, outputDir, suff, True, 0, 50, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="rho", yTitle="# Events")

		# vtxHistoNames = ["vtx_x", "vtx_y", "vtx_z"]
		# for histoName in vtxHistoNames:
			# plotDataAndMCHistosAndRatio(histoName+"_histo"+triggerName, "", dataInputFile, MCFileList, outputDir, "", False, 0, 0, False, 0, 0, legendData, legendList, True, False, False, False, norm, colorList)


		objNames = ["ele", "goodEle", "jet", "goodJet"]#, "muon", "goodMuon", "muonRochCor", "goodMuonRochCor"]
		for objName in objNames:
			# plotDataAndMCHistosAndRatio("pt_"+objName+triggerName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 700, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T} "+objName, yTitle="# Events")
			plotDataAndMCHistosAndRatio("pt_"+objName+triggerName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, True, 0, 300, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="p_{T} "+objName, yTitle="# Events")
			plotDataAndMCHistosAndRatio("eta_"+objName+triggerName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#eta "+objName, yTitle="# Events")
			plotDataAndMCHistosAndRatio("phi_"+objName+triggerName+"_histo", "", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="#phi "+objName, yTitle="# Events")


		if dataDir == "DoubleEG" or dataDir == "SingleMuon":
			etaMassNames = ["", "_EB-EB", "_EE-EE", "_EB-EE", "_noEB-EB"]

			Znames = ["toEE_MLMJelectrons","toEE_electrons","toEE_electronsAndJets"]
			if dataDir == "SingleMuon":
				Znames = ["toMuMu_MLMJmuons","toMuMu_muons","toMuMu_muonsAndJets"]

			for Zname in Znames:
				for etaMassName in etaMassNames:
					plotDataAndMCHistosAndRatio("Z"+Zname+triggerName+"_mass", etaMassName+"_histo", dataInputFile, MCFileList, outputDir, suff, True, 70, 110, False, 0, 0, legendData, legendList, False, log, False, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="m_{ll}", yTitle="# Events")	
				plotDataAndMCHistosAndRatio("Z"+Zname+triggerName+"_pt", "_histo", dataInputFile, MCFileList, outputDir, suff, False, 0, 0, False, 0, 0, legendData, legendList, False, log, True, False, norm, colorList, xsecList, sumWeightsList, lumi, title="comparison", xTitle="pT_{ll}", yTitle="# Events")	






# etaNames = ["", "_EB", "_EE"]
# eleNames = ["leadingEle", "subLeadingEle", "goodLeadingEle", "goodSubLeadingEle"] #,"ele", "goodEle"]

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



