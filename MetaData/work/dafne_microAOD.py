import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100)) #-1 ) ) 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

import os
if os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_2016_TrancheIV_v7')
else:
    raise Exception,"The default setup for microAODstd.py does not support releases other than 80X"

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
process.RandomNumberGeneratorService.flashggRandomizedPhotons = cms.PSet(
          initialSeed = cms.untracked.uint32(16253245)
        )
process.RandomNumberGeneratorService.flashggRandomizedElectrons = cms.PSet(
          initialSeed = cms.untracked.uint32(16253245)
        # engineName = cms.untracked.string('TRandom3') # optional, default to HepJamesRandom if absent
        )

#80x signal
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/RunIISummer16MiniAODv2/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CE51F2A9-81CA-E611-89CE-F04DA27540CA.root"))

#80x data  
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver1-v1/100000/02C07D99-20EB-E611-92B2-3417EBE700D2.root"))


process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now
# process.MessageLogger.suppressWarning.extend(['SimpleMemoryCheck','MemoryCheck']) # this would have been better...

# Uncomment the following if you notice you have a memory leak
# This is a lightweight tool to digg further
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#                                        ignoreTotal = cms.untracked.int32(1),
#                                        monitorPssAndPrivate = cms.untracked.bool(True)
#                                       )
 
process.load("dafne/MicroAOD/flashggMicroAODSequence_MultiLeptonMultiJet_cff")  

# NEEDED FOR ANYTHING PRIOR TO reMiniAOD
#process.weightsCount.pileupInfo = "addPileupInfo"

from dafne.MicroAOD.flashggMicroAODOutputCommands_cff import microAODmultiLeptonMultiJetOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('dafneMicroAOD.root'),
                               outputCommands = microAODmultiLeptonMultiJetOutputCommand
                               )


process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.globalTightHalo2016Filter_cfi')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')

process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)


process.p = cms.Path(process.flashggMicroAODSequenceMultiLeptonMultiJet)
process.e = cms.EndPath(process.out)


from flashgg.MicroAOD.MicroAODCustomize import customize
customize(process)
