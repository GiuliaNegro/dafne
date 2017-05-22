import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os, sys

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

process = cms.Process("FlashggAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

if os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
	process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_2016_TrancheIV_v7')
else:
	raise Exception,"The default setup does not support releases other than 80X"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )


## input file
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016-Moriond17/dafne-Moriond17/DoubleEG/cmsWR2016-Moriond17-dafne-Moriond17-v0-Run2016F-03Feb2017-v1/170227_174902/0000/dafneMicroAOD_701.root"
				"root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016-Moriond17/dafne-Moriond17/DoubleEG/cmsWR2016-Moriond17-dafne-Moriond17-v0-Run2016D-03Feb2017-v1/170227_174446/0000/dafneMicroAOD_1.root"
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016-Moriond17/dafne-Moriond17/DoubleEG/cmsWR2016-Moriond17-dafne-Moriond17-v0-Run2016D-03Feb2017-v1/170227_174446/0001/dafneMicroAOD_1159.root"

				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016-Moriond17/dafne-Moriond17/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/cmsWR2016-Moriond17-dafne-Moriond17-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/170228_091616/0000/dafneMicroAOD_1.root"
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016-Moriond17-DYinclusive/dafne-Moriond17/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/cmsWR2016-Moriond17-DYinclusive-dafne-Moriond17-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/170403_215703/0000/dafneMicroAOD_103.root"
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016-Moriond17-DYinclusive/dafne-Moriond17/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/cmsWR2016-Moriond17-DYinclusive-dafne-Moriond17-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/170403_215518/0000/dafneMicroAOD_1.root"
))


# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()


## global variables to dump
from flashgg.Taggers.globalVariables_cff import globalVariables


## analyzer
process.analysisTree = cms.EDAnalyzer('EDminiTreeMaker_multiLeptonMultiJet',
										genParticleTag = cms.InputTag( "flashggPrunedGenParticles" ),  
										generatorInfo = cms.InputTag('generator'),  																			
										PileUpTag = cms.InputTag('slimmedAddPileupInfo'),
										VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
										MultiLeptonMultiJetTag=cms.InputTag('flashggMultiLeptonMultiJetSystematics'), 
										JetsTag=cms.InputTag('flashggJetSystematics'),
										GenJetTag=cms.InputTag( "slimmedGenJets"),
										ElectronTag=cms.InputTag('flashggEleSystematics'),
										MuonTag=cms.InputTag('flashggMuonSystematics'),
										triggerBits = cms.InputTag('TriggerResults::HLT'),
										rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll'),	
										lumiWeight=cms.untracked.double(1000.),		
										saveHEEPvariables=cms.untracked.bool(False), 	
										isDoubleEGinSignalRegion=cms.untracked.bool(customize.isDoubleEGinSignalRegion), 
										globalVariables = globalVariables
										)


## Systematics        
## import systs. customize
from dafne.Systematics.SystematicsCustomize_dafne import *

## load syst producer
process.load("dafne.Systematics.flashggMultiLeptonMultiJetSystematics_cfi")
process.load("dafne.Systematics.flashggMuonFromMLMJSystematics_cfi")
process.load("flashgg.Systematics.flashggElectronSystematics_cfi")
process.load("dafne.Systematics.flashggMuonSystematics_cfi")
process.load("dafne.Systematics.flashggJetSystematics_cfi")

## if data apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == "data":
	print 'data' 
	customizeMultiLeptonMultiJetSystematicsForData(process)  # only central value, no syst. shifts 
	customizeSingleObjectSystematicsForData(process)
	print 'customization done'
else:
	print 'mc'
	customizeMultiLeptonMultiJetSystematicsForMC(process)  # only central value, no syst. shifts 
	customizeSingleObjectSystematicsForMC(process)
	print 'customization done'

for pset2 in process.flashggEleSystematics.SystMethods2D:
		pset2.RandomLabel = "rnd_g_E"


debugEle = False
if debugEle: 
	if customize.processType == "data":
		for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
			if pset.Label != "JEC"  and pset.Label != "JER":
				pset.Debug = True
				pset.ExaggerateShiftUp = True
	else:
		for pset2 in process.flashggMultiLeptonMultiJetSystematics.SystMethods2D:
			pset2.Debug = True
			pset2.ExaggerateShiftUp = True

debugJet = False
if debugJet: 
	for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
		if pset.Label == "JEC":
			print "JEC"
			pset.Debug = True
		if pset.Label == "JER":
			print "JER"
			pset.Debug = True
		
debugMuon = False
if debugMuon: 
	for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
		# if pset.Label == "MuonIDSF":
		# 	print "MuonIDSF"
		# 	pset.Debug = True
		# if pset.Label == "MuonIsoSF":
		# 	print "MuonIsoSF"
		# 	pset.Debug = True
		if pset.Label == "MuonTriggerSF":
			print "MuonTriggerSF"
			pset.Debug = True

debugSingleEle = False
if debugSingleEle: 
	if customize.processType == "data":
		for pset in process.flashggEleSystematics.SystMethods:
			pset.Debug = True
			pset.ExaggerateShiftUp = True
	else:
		for pset2 in process.flashggEleSystematics.SystMethods2D:
			pset2.Debug = True
			pset2.ExaggerateShiftUp = True

debugSingleMuon = False
if debugSingleMuon:
	for pset in process.flashggMuonSystematics.SystMethods:
		if pset.Label == "SingleMuonIDSF":
			print "SingleMuonIDSF"
			pset.Debug = True
		# if pset.Label == "SingleMuonIsoSF":
		# 	print "SingleMuonIsoSF"
		# 	pset.Debug = True
		# if pset.Label == "SingleMuonTriggerSF":
		# 	print "SingleMuonTriggerSF"
		# 	pset.Debug = True		

debugSingleJet = False
if debugSingleJet: 
	for pset in process.flashggJetSystematics.SystMethods:
		if pset.Label == "JEC":
			print "JEC"
			pset.Debug = True
		if pset.Label == "JER":
			print "JER"
			pset.Debug = True


process.TFileService = cms.Service("TFileService",
									fileName = cms.string("mytree.root")
									)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )



process.load('dafne.Validation.hltFilters_cff')

if customize.trigger=="signalEE": 
	process.HltSequence = cms.Sequence(process.wReejjHLTFilterMW)
	print "trigger = signalEE"
	print customize.isDoubleEGinSignalRegion
	if customize.isDoubleEGinSignalRegion==True:
		process.HltSequence = cms.Sequence(process.wReejjHLTFilterMW + process.wReejjHLTFilterGsfTrkIdVL) 
elif customize.trigger=="signalMuMu":
	process.HltSequence = cms.Sequence(process.wRmumujjHLTFilter)
	print "trigger = signalMuMu"
elif customize.trigger=="eMuSideband":
	process.HltSequence = cms.Sequence(process.wRemujjHLTFilter) 
	print "trigger = flavorSideband"
elif customize.trigger=="TnPEE":
	process.HltSequence = cms.Sequence(process.tagAndProbeDoubleEleHLTFilter)       
	print "trigger = TnPEE"
elif customize.trigger=="TnPMuMu":
	process.HltSequence = cms.Sequence(process.tagAndProbeDoubleMuHLTFilter)
	print "trigger = TnPMuMu"


process.fullSeq = cms.Sequence( process.flashggEleSystematics * process.flashggMuonSystematics * process.jetCorrectorChain * process.flashggJetSystematics * process.flashggMultiLeptonMultiJetSystematics * process.analysisTree)
process.p = cms.Path(process.HltSequence * process.fullSeq)


## set default options if needed
customize.setDefault("maxEvents", 1000)
customize.setDefault("targetLumi",1e+3)
## call the customization
customize(process)

