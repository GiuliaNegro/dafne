import ROOT
import sys
from array import array
from math import sqrt
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        print "Usage: %s <dir1, dir2, output_dir, legend1, legend2, isDataMC, log, suff>" % sys.argv[0]  
        sys.exit(1)

#python Analysis/makePlots_data-MC.py DoubleEG/ DYJetsPtBinned Comparisons-Moriond17/DoubleEG-DY data DY True False "" -b
#python Analysis/makePlots_data-MC.py DoubleEG/ DYJetsPtBinned Comparisons-Moriond17/DoubleEG-DY data DY True True "_log" -b

#python Analysis/makePlots_data-MC.py DoubleEG/ DYJetsPtBinned/DY50-1 Comparisons-Moriond17/DoubleEG-DY data DY True False "" -b
#python Analysis/makePlots_data-MC.py DoubleEG/ DYJetsPtBinned/prova Comparisons-Moriond17/DoubleEG-DY/prova data DY True False "" -b


dir1 = str(sys.argv[1])
dir2 = str(sys.argv[2])
output_dir = str(sys.argv[3])
# output_dir = "gnegro@lxplus.cern.ch:/afs/cern.ch/user/g/gnegro/www/cmsWR/"+str(sys.argv[3])
legend1 = str(sys.argv[4])
legend2 = str(sys.argv[5])
dataMC = sys.argv[6] #True
log = str(sys.argv[7])
suff = str(sys.argv[8])

inputDir = "Distributions-Moriond17/"+dir1

inputfile1 = TFile(inputDir+"/distributions_histos.root")
inputfile2 = TFile(inputDir+dir2+"/distributions_histos.root")

lumi = 35.867
CMS_lumi.lumi_13TeV = str(lumi)+" fb^{-1}"




etaMassNames = [""]#, "_EB-EB", "_EE-EE", "_EB-EE", "_noEB-EB"]
regionNames = ["_lowMllCR"]#,"_lowMlljjCR", "_eMuSidebandCR", "_signalRegion"]
objNames = ["leadingEle", "subLeadingEle", "leadingJet", "subLeadingJet"]#"leadingLepton", "subLeadingLepton", "leadingMu", "subLeadingMu"]


for regionName in regionNames:

	for objName in objNames:
		plot2HistosAndRatio("pt_"+objName+regionName+"_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 300, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="pt_"+objName, yTitle="[a.u.]")
		plot2HistosAndRatio("eta_"+objName+regionName+"_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="#eta_"+objName, yTitle="[a.u.]")
		plot2HistosAndRatio("phi_"+objName+regionName+"_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="#phi "+objName, yTitle="[a.u.]")


	for etaMassName in etaMassNames:
		# if "lowMlljjCR" in regionName:
		# 	plot2HistosAndRatio("mass_multiLeptonMultiJets"+regionName, etaMassName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 1000, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="m_{lljj}", yTitle="[a.u.]")
		# else:
		plot2HistosAndRatio("mass_multiLeptonMultiJets"+regionName, etaMassName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 3000, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="m_{lljj}", yTitle="[a.u.]")

		plot2HistosAndRatio("mass_dileptons"+regionName, etaMassName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 800, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="m_{ll}", yTitle="[a.u.]")
		plot2HistosAndRatio("mass_dijets"+regionName, etaMassName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 800, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="m_{ll}", yTitle="[a.u.]")


	plot2HistosAndRatio("pt_dileptons_"+regionName, etaMassName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 800, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="p_{T}(ll)", yTitle="[a.u.]")



# objNames = ["ele", "goodEle", "jet", "goodJet"]#, "muon", "goodMuon", "muonRochCor", "goodMuonRochCor"
# for objName in objNames:
# 	plot2HistosAndRatio("pt_"+objName+"_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 250, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="pt_"+objName, yTitle="[a.u.]")
# 	plot2HistosAndRatio("eta_"+objName+"_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="#eta_"+objName, yTitle="[a.u.]")
# 	plot2HistosAndRatio("phi_"+objName+"_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="#phi "+objName, yTitle="[a.u.]")


plot2HistosAndRatio("nvtx_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="nvtx_histo_comparison", xTitle="nvtx", yTitle="[a.u.]")
plot2HistosAndRatio("rho_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="rho_histo_comparison", xTitle="rho", yTitle="[a.u.]")

# vtxHistoNames = ["vtx_x", "vtx_y", "vtx_z"]
# for histoName in vtxHistoNames:
# 	plot2HistosAndRatio(histoName+"_histo", "", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, log, False, False, False, False, title=histoName+"_comparison", xTitle=histoName, yTitle="[a.u.]")




file1 = TFile(inputDir+"/Zmass_histos.root")
file2 = TFile(inputDir+dir2+"/Zmass_histos.root")

etaMassNames = ["", "_EB-EB", "_EE-EE", "_EB-EE", "_noEB-EB"]
Znames = ["toEE"]#_MLMJelectrons"]#, "toEE_electrons"]
if dir1 == "SingleMuon/":
	Znames = ["toMuMu"]#_MLMJmuons"]#, "toMuMu_muons"

for etaMassName in etaMassNames:
	for Zname in Znames:
 		plot2HistosAndRatio("Z"+Zname+"_mass", etaMassName+"_histo", file1, file2, output_dir, suff, dataMC, True, 70, 110, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title="comparison", xTitle="m_{ll}", yTitle="[a.u.]")








# etaNames = ["", "_EB", "_EE"]
# eleNames = ["ele", "goodEle", "leadingEle", "subLeadingEle", "goodLeadingEle", "goodSubLeadingEle"]

# for eleName in eleNames:
# 	for etaName in etaNames:

# 		plot2HistosAndRatio("etaSC_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="etaSC_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("isEcalDriven_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="isEcalDriven_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("dEtaIn_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, -0.1, 0.1, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="dEtaIn_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("dPhiIn_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, -0.1, 0.1, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="dPhiIn_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("hOverE_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 0.3, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="hOverE_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("sigmaIetaIeta_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="sigmaIetaIeta_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("e2x5_e5x5_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="e2x5_e5x5_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("e1x5_e5x5_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, True, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="e1x5_e5x5_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("EmHadDepth1Iso_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 20, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="EmHadDepth1Iso_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("missingHits_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, False, 0, 0, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="missingHits_"+eleName+etaName, yTitle="[a.u.]")

# 		plot2HistosAndRatio("dxy_"+eleName, etaName+"_histo", inputfile1, inputfile2, output_dir, suff, dataMC, True, 0, 0.05, False, 0, 0, legend1, legend2, False, False, log, False, False, False, False, title=histoName+"_comparison", xTitle="dxy_"+eleName+etaName, yTitle="[a.u.]")





# histoName, etaRegion, inputfile1, inputfile2, output_name, suff, dataMC, zoomX, xRangeMin, xRangeMax, zoomY, yRangeMin, yRangeMax, leg1_name, leg2_name, leftLegends, profile, log, doRebin, doRebinVariableBinSize, doRestrictedIntegral, MCreweighted, title="", xTitle="", yTitle=""):

