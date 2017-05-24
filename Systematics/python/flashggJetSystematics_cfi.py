import FWCore.ParameterSet.Config as cms
from os import environ


flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
									JetsTag = cms.InputTag("flashggFinalJets"),
									NCollections = cms.uint32(1)
									)

singleJEC = cms.PSet( MethodName = cms.string("FlashggJetEnergyCorrector"),
			    Label = cms.string("JEC"),
				NSigmas = cms.vint32(-1,1),
				OverallRange = cms.string("abs(eta)<5.0"),
				Debug = cms.untracked.bool(False),
				ApplyCentralValue = cms.bool(True),
				SetupUncertainties = cms.bool(True),
				# JetCorrectorTag = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
				JetCorrectorTag = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector")
				)

singleJER = cms.PSet( MethodName = cms.string("FlashggJetSmear"),
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


flashggJetSystematics = cms.EDProducer('FlashggJetSystematicProducer',
					src = cms.InputTag('flashggUnpackedJets',str(0)),
					SystMethods2D = cms.VPSet(),
					SystMethods = cms.VPSet(singleJEC) #used in both data and MC
)


