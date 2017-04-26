import FWCore.ParameterSet.Config as cms
from dafne.Systematics.flashggMuonFromMLMJSystematics_cfi import binID, binIso, binTrigger, binEffData

flashggMuonSystematics = cms.EDProducer('FlashggMuonSystematicProducer',
					src = cms.InputTag("flashggSelectedMuons"),
					SystMethods2D = cms.VPSet(),
					SystMethods = cms.VPSet(
									cms.PSet( MethodName = cms.string("FlashggMuonWeight"),
											  Label = cms.string("SingleMuonIDSF"),
											  NSigmas = cms.vint32(-1,1),
											  OverallRange = cms.string("abs(eta)<2.4"),
											  BinList = binID,
											  Debug = cms.untracked.bool(False),
											  ApplyCentralValue = cms.bool(True)
									),
									cms.PSet( MethodName = cms.string("FlashggMuonWeight"),
											  Label = cms.string("SingleMuonIsoSF"),
											  NSigmas = cms.vint32(-1,1),
											  OverallRange = cms.string("abs(eta)<2.4"),
											  BinList = binIso,
											  Debug = cms.untracked.bool(False),
											  ApplyCentralValue = cms.bool(True)
									),
									# cms.PSet( MethodName = cms.string("FlashggMuonWeight"),
									# 		  Label = cms.string("SingleMuonTriggerSF"),
									# 		  NSigmas = cms.vint32(-1,1),
									# 		  OverallRange = cms.string("abs(eta)<2.4"),
									# 		  BinList = binTrigger, 
									# 		  BinList2 = binEffData, 
									# 		  Debug = cms.untracked.bool(False),
									# 		  ApplyCentralValue = cms.bool(True)
									# )									
					)
)


