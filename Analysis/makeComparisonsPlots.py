import ROOT
import sys
from array import array
from math import sqrt
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        print "Usage: %s <dir1, dir2, output_dir, legend1, legend2, isdataMC>" % sys.argv[0]  
        sys.exit(1)

#python makeComparisonsPlots.py dir1 dir2 output_dir isDataMC -b


dir1 = str(sys.argv[1])
dir2 = str(sys.argv[2])
output_dir = str(sys.argv[3])
legend1 = str(sys.argv[4])
legend2 = str(sys.argv[5])
dataMC = sys.argv[6] #True

inputfile1 = TFile(dir1+"/distributions_histos.root")
inputfile2 = TFile(dir2+"/distributions_histos.root")



etaNames = ["", "_EB", "_EE"]
etaMassNames = ["", "_EB-EB", "_EE-EE", "_EB-EE", "_noEB-EB"]

vtxHistoNames = ["vtx_x", "vtx_y", "vtx_z"]
objNames = ["ele", "leadingEle", "subLeadingEle", "muons", "leadingMu", "subLeadingMu", "jet", "leadingJet", "subLeadingJet"]
# eleNames = ["ele", "leadingEle", "subLeadingEle", "elePassingHeepId", "leadingElePassingHeepId", "subLeadingElePassingHeepId", "elePassingEleId", "leadingElePassingEleId", "subLeadingElePassingEleId"]
eleNames = ["elePassingHeepId", "leadingElePassingHeepId", "subLeadingElePassingHeepId", "elePassingEleId", "leadingElePassingEleId", "subLeadingElePassingEleId"]
dldjNames = ["dldj","dldjPreselected","dldjPassingHeepId","dldjPassingEleId"]
Znames = ["toEEpassingHeepId", "toEEpassingEleId"]#,"toEE", "toMuMu"]




plot2HistosAndRatio("nvtx_histo", "", inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title="nvtx_histo_comparison", xTitle="nvtx", yTitle="[a.u.]")


for histoName in vtxHistoNames:
	plot2HistosAndRatio(histoName+"_histo", "", inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, False, False, False, False, False, title=histoName+"_comparison", xTitle=histoName, yTitle="[a.u.]")


for objName in objNames:
	plot2HistosAndRatio("pt_"+objName+"_histo", "", inputfile1, inputfile2, output_dir, "", dataMC, True, 0, 250, False, 0, 0, legend1, legend2, False, False, True, False, False, False, False, title=histoName+"_comparison", xTitle="pt_"+objName, yTitle="[a.u.]")
	plot2HistosAndRatio("eta_"+objName+"_histo", "", inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="eta_"+objName, yTitle="[a.u.]")


for etaMassName in etaMassNames:

	for dldjName in dldjNames:
		massDldjHistoNames = ["mass_dileptondijets_histo_"+dldjName, "mass_dileptons_histo_"+dldjName]#, "mass_dijets_histo_"+dldjName]  #,"mass_dijetsLeadingLepton_histo_"+dldjName, "mass_dijetsSubLeadingLepton_histo_"+dldjName] 

		for histoName in massDldjHistoNames:
			if "mass_dileptons_histo_" in histoName:
				plot2HistosAndRatio(histoName, etaMassName, inputfile1, inputfile2, output_dir, "", dataMC, True, 0, 1000, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle=histoName, yTitle="[a.u.]")
			else:		
				plot2HistosAndRatio(histoName, etaMassName, inputfile1, inputfile2, output_dir, "", dataMC, True, 0, 2000, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle=histoName, yTitle="[a.u.]")

	for Zname in Znames:
		ZHistoNames = ["Z"+Zname+"_mass_histo"]

		for histoName in ZHistoNames:
			plot2HistosAndRatio(histoName, etaMassName, inputfile1, inputfile2, output_dir, "", dataMC, True, 70, 110, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle=histoName, yTitle="[a.u.]")


for eleName in eleNames:
	for etaName in etaNames:

		plot2HistosAndRatio("etaSC_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="etaSC_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("isEcalDriven_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="isEcalDriven_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("dEtaIn_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, True, -0.1, 0.1, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="dEtaIn_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("dPhiIn_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, True, -0.1, 0.1, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="dPhiIn_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("hOverE_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, True, 0, 0.3, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="hOverE_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("sigmaIetaIeta_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="sigmaIetaIeta_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("e2x5_e5x5_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="e2x5_e5x5_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("e1x5_e5x5_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="e1x5_e5x5_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("EmHadDepth1Iso_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, True, 0, 20, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="EmHadDepth1Iso_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("missingHits_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="missingHits_"+eleName+etaName, yTitle="[a.u.]")

		plot2HistosAndRatio("dxy_"+eleName+"_histo", etaName, inputfile1, inputfile2, output_dir, "", dataMC, True, 0, 0.05, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title=histoName+"_comparison", xTitle="dxy_"+eleName+etaName, yTitle="[a.u.]")





# histoName, etaRegion, inputfile1, inputfile2, output_name, suff, dataMC, zoomX, xRangeMin, xRangeMax, zoomY, yRangeMin, yRangeMax, leg1_name, leg2_name, leftLegends, profile, log, doRebin, doRebinVariableBinSize, doRestrictedIntegral, MCreweighted, title="", xTitle="", yTitle=""):

