import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100) )
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

#80x MC
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
	#signal 2016
	# "/store/mc/RunIISummer16MiniAODv2/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CE51F2A9-81CA-E611-89CE-F04DA27540CA.root"
	# "/store/mc/RunIISummer16MiniAODv2/WRToNuMuToMuMuJJ_MW-1200_MNu-600_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1C38B200-83CE-E611-94C0-0CC47A745298.root"

	#bkg
	"/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00933E2A-A0D5-E611-B2CD-00266CF89130.root"
	# "/store/mc/RunIISummer16MiniAODv2/WZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/60000/002EC119-D4DA-E611-A902-008CFA1979FC.root"
	# "/store/mc/RunIISummer16MiniAODv2/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/02C8895F-E8DA-E611-8CC9-0023AEEEB55F.root"
	# "/store/mc/RunIISummer16MiniAODv2/WW_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/70000/00155A96-32DA-E611-8F20-001E67580BAC.root"
	# "/store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/08D759B1-CBB6-E611-87B3-484D7E8DF0D3.root"
	# "/store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/16E7A136-1FBE-E611-B8DA-0025905C3DD0.root"
	# "/store/mc/RunIISummer16MiniAODv2/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0071A1CE-DDD5-E611-AC2C-001E67397CB5.root"
	# "/store/mc/RunIISummer16MiniAODv2/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/000161D3-22C4-E611-BBCB-002590494DE8.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0073209B-97C9-E611-A6D0-008CFA5D2758.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/006F1B9A-81D0-E611-B9CE-0025905AA9CC.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/60000/003FAFE5-4FCB-E611-99AD-00237DF28460.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0AD9C544-FCD2-E611-B4EE-008CFA197E0C.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/06B39D0E-F0C6-E611-AEC3-001CC47D420E.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1C6C4217-42D0-E611-984C-0CC47A537688.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1C6C4217-42D0-E611-984C-0CC47A537688.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/0419C5FF-03BE-E611-B327-001E67504F55.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/242AE947-8CC0-E611-878C-002590FD5A4C.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/061DAB15-66BE-E611-B8D6-0CC47AA99438.root"
	# "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/0C99DA93-79C1-E611-8D0B-FA163EEA85D4.root"
	# "/store/mc/RunIISummer16MiniAODv2/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/283E8320-42CB-E611-86C6-02163E019B7C.root"
	# "/store/mc/RunIISummer16MiniAODv2/ST_tWnunu_5f_LO_13TeV-MadGraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3201F609-DCC6-E611-9DE2-02163E00B15C.root"
	# "/store/mc/RunIISummer16MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/0A6E003C-CABB-E611-956C-0025905B85AE.root?"
	# "/store/mc/RunIISummer16MiniAODv2/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/3495D426-73C1-E611-B11B-0CC47A4D764A.root"
	))

#80x data
# process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
	# "/store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver1-v1/100000/02C07D99-20EB-E611-92B2-3417EBE700D2.root"
	# "/store/data/Run2016C/SingleMuon/MINIAOD/03Feb2017-v1/50000/001CF316-1AEB-E611-BBBD-0CC47A4C8EE2.root"
	# "/store/data/Run2016D/MuonEG/MINIAOD/03Feb2017-v1/80000/02264DFC-6EEB-E611-95AF-0090FAA572B0.root"
	# ))


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

# from flashgg.MicroAOD.flashggMicroAODOutputCommands_cff import microAODDefaultOutputCommand
# process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myMicroAODOutputFile_DiLeptonDiJet.root'),
#                                outputCommands = microAODDefaultOutputCommand
#                                )

from dafne.MicroAOD.flashggMicroAODOutputCommands_cff import microAODmultiLeptonMultiJetOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myMicroAODOutputFile_MultiLeptonMultiJet.root'),
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

if "DY" in customize.datasetName or "SingleElectron" in customize.datasetName or "DoubleEG" in customize.datasetName:
  customize.customizeHLT(process)