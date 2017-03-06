import FWCore.ParameterSet.Config as cms
from os import environ


emptyBins = cms.PSet(
	variables = cms.vstring("1"),
	bins = cms.VPSet()
	)


scalesAndSmearingsPrefix = cms.string("EgammaAnalysis/ElectronTools/data/Winter_2016_reReco_v1_ele")


SmearHighR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet2D"),
		Label = cms.string("SmearHighR9EB"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the electron object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)

SmearLowR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet2D"),
		Label = cms.string("SmearLowR9EB"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the electron object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)

SmearHighR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),  
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet2D"),
		Label = cms.string("SmearHighR9EE"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the electron object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)

SmearLowR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet2D"),
		Label = cms.string("SmearLowR9EE"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the electron object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)


ScaleHighR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet"),
		Label = cms.string("ScaleHighR9EB"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)

ScaleLowR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet"),
		Label = cms.string("ScaleLowR9EB"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)

ScaleHighR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet"),
		Label = cms.string("ScaleHighR9EE"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)

ScaleLowR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromMultiLeptonMultiJet"),
		Label = cms.string("ScaleLowR9EE"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)

JEC = cms.PSet( JetMethodName = cms.string("FlashggJetEnergyCorrector"),
		MethodName = cms.string("FlashggJetFromMultiLeptonMultiJet"),
		Label = cms.string("JEC"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("abs(eta)<5.0"),
		Debug = cms.untracked.bool(False),
		ApplyCentralValue = cms.bool(True),
		SetupUncertainties = cms.bool(True),
		# JetCorrectorTag = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
		JetCorrectorTag = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector")
		)


JER = cms.PSet( JetMethodName = cms.string("FlashggJetSmear"),
		MethodName = cms.string("FlashggJetFromMultiLeptonMultiJet"),
		Label = cms.string("JER"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("abs(eta)<5.0"),
		RandomLabel = cms.string("rnd_g_JER"), # for no-match case
		rho = cms.InputTag('fixedGridRhoAll'),
		Debug = cms.untracked.bool(False),
		ApplyCentralValue = cms.bool(True),
		UseTextFiles = cms.bool(True),
		TextFileResolution = cms.string("%s/src/flashgg/Systematics/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt" % environ['CMSSW_BASE']),
		TextFileSF = cms.string("%s/src/flashgg/Systematics/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt" % environ['CMSSW_BASE'])
		)



flashggMultiLeptonMultiJetSystematics = cms.EDProducer('FlashggMultiLeptonMultiJetSystematicProducer',
		src = cms.InputTag("flashggMultiLeptonMultiJet"),
		SystMethods2D = cms.VPSet(
				SmearHighR9EB_EGM,
				SmearLowR9EB_EGM,
				SmearHighR9EE_EGM,
				SmearLowR9EE_EGM			
		),
		# the number of syst methods matches the number of nuisance parameters
		# assumed for a given systematic uncertainty and is NOT required
		# to match 1-to-1 the number of bins above.
		SystMethods = cms.VPSet(
				ScaleHighR9EB_EGM,
				ScaleLowR9EB_EGM,
				ScaleHighR9EE_EGM,
				ScaleLowR9EE_EGM
		)
)





