import FWCore.ParameterSet.Config as cms


def customizeMultiLeptonMultiJetSystematicsForData(process):
	#EGM scale 1D: put in central value, but omit shifts

	electronScaleBinsData = getattr(process,'electronScaleBinsData',None)
	if hasattr(process,'electronScaleBinsData'):
		print electronScaleBinsData, process.electronScaleBinsData

	for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
		pset.ApplyCentralValue = cms.bool(True) # Turn on central shift 
		pset.NSigmas = cms.vint32() # Do not perform up/down syst shift (1D case)
		if electronScaleBinsData != None: 
				pset.BinList = electronScaleBinsData


def customizeMultiLeptonMultiJetSystematicsForMC(process):
	#EGM smearing 2D: put in central value, but omit shifts

	for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods2D:
		pset.ApplyCentralValue = cms.bool(True)  # Turn on central shift 
		pset.NSigmas = cms.PSet( firstVar = cms.vint32(), secondVar = cms.vint32() ) # Do not perform up/down syst shifts (2D case)





