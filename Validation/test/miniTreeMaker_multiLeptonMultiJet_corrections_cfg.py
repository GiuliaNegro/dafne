import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

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
				#data
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016/dafne/DoubleEG/cmsWR2016-dafne-v0-Run2016B-PromptReco-v2/161027_122605/0000/dafneMicroAOD_1.root"
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016-ReReco/dafne-v1/DoubleEG/cmsWR2016-ReReco-dafne-v1-v0-Run2016B-23Sep2016-v3/170105_182103/0000/dafneMicroAOD_1.root"

				# mc
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016/dafne/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/cmsWR2016-dafne-v0-gnegro-WR-3200_ToLNu-1600_ToEEJJ_miniAOD_13TeV-2016-b59cb78551aff289588aa5c69db4a3a1/161027_124028/0000/dafneMicroAOD_1.root"
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016/dafne/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/cmsWR2016-dafne-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/161027_121933/0000/dafneMicroAOD_107.root"
				# "/store/user/gnegro/cmsWR/cmsWR2016-signal/dafne-v1/WRToEEJJ_MW-2400_MNu-1200_TuneCUETP8M1_13TeV-pythia8/cmsWR2016-signal-dafne-v1-v0-gnegro-RunIIWinter16_80X_mcRun2_asymptotic_2016_miniAODv2_v1_MINIAODSIM-9232bfa9d7b25477dcde67f5060ed55b/170106_094339/0000/dafneMicroAOD_1.root"				
				# "/store/user/gnegro/cmsWR/cmsWR2016-bkg/dafne-v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/cmsWR2016-bkg-dafne-v1-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/170106_095719/0000/dafneMicroAOD_1.root"
				"file:/afs/cern.ch/work/g/gnegro/NuAnalysis-Moriond17/CMSSW_8_0_26_patch1/src/flashgg/myMicroAODOutputFile_MultiLeptonMultiJet.root"
))


## HLT filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
	  	 	"HLT_DoubleEle33_CaloIdL_MW_v*",
	   		"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
		   	"HLT_Mu50_v*",
		   	"HLT_TkMu50_v*",
		   	"HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v*",
		   	"HLT_Ele27_WPTight_Gsf_v*",
		   	"HLT_IsoMu24_v*",
		   	"HLT_IsoMu27_v*"
			),
			throw = False
			# andOr    = cms.bool(True) # True = or between triggers 
)


# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()


# # flashgg tag sequence (for dipho MVA) and jet collections   
process.load("flashgg/Taggers/flashggTagSequence_cfi")
# process.load("dafne/Taggers/flashggDiLeptonDiJetTagSequence_cfi") #uso questa?
# process.flashggTagSequence.remove(process.flashggUpdatedIdMVADiPhotons) # Needs to be run before systematics
# massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggUpdatedIdMVADiPhotons"),cms.InputTag("flashggDiPhotonSystematics"))

# #remove un-necessary tags ...
# process.flashggTagSequence.remove(process.flashggVBFTag)
# process.flashggTagSequence.remove(process.flashggTTHLeptonicTag)
# process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
# process.flashggTagSequence.remove(process.flashggVHEtTag)
# process.flashggTagSequence.remove(process.flashggVHLooseTag)
# process.flashggTagSequence.remove(process.flashggVHTightTag)
# process.flashggTagSequence.remove(process.flashggTagSorter)

## global variables to dump
from flashgg.Taggers.globalVariables_cff import globalVariables

## analyzer
process.analysisTree = cms.EDAnalyzer('EDminiTreeMaker_multiLeptonMultiJet',
										genParticleTag = cms.InputTag( "flashggPrunedGenParticles" ),  
										generatorInfo = cms.InputTag('generator'),  																			
										PileUpTag = cms.InputTag('slimmedAddPileupInfo'),
										VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
										MultiLeptonMultiJetTag=cms.InputTag('flashggMultiLeptonMultiJet'), 
										JetsTag=cms.InputTag('flashggFinalJets'),
										GenJetTag=cms.InputTag( "slimmedGenJets"),
										ElectronTag=cms.InputTag('flashggSelectedElectrons'),
										MuonTag=cms.InputTag('flashggSelectedMuons'),
										triggerBits = cms.InputTag('TriggerResults::HLT'),
										rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll'),	
										bTag = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),  
										isControlSample = cms.untracked.bool(False),
										lumiWeight=cms.untracked.double(1000.),																		
										globalVariables = globalVariables
										)


process.analysisTree.globalVariables.addTriggerBits = cms.PSet(
		tag = cms.InputTag("TriggerResults::HLT"),
		bits = cms.vstring(
          	"HLT_DoubleEle33_CaloIdL_MW",
			"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW",
			"HLT_Mu50",
			"HLT_TkMu50",
			"HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL",
			"HLT_Ele27_WPTight_Gsf",
			"HLT_IsoMu24",
			"HLT_IsoMu27"
		)
)

process.analysisTree.triggerBits = cms.InputTag('TriggerResults::HLT') 
# process.analysisTree.globalVariables.addTriggerBits.tag = cms.InputTag("TriggerResults::HLT")



## Systematics        
## import systs. customize
# from flashgg.Systematics.SystematicsCustomize import *
from dafne.Systematics.SystematicsCustomize_dafne import *

## load syst producer
process.load("dafne.Systematics.flashggMultiLeptonMultiJetSystematics_cfi")


debug = True

if debug: 
	if customize.processType == "data":
		for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
			pset.Debug = True
			pset.ExaggerateShiftUp = True
	else:
		for pset2 in process.flashggMultiLeptonMultiJetSystematics.SystMethods2D:
			pset2.Debug = True
			pset2.ExaggerateShiftUp = True


## if data apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == "data":
    print 'data' 
    customizeMultiLeptonMultiJetSystematicsForData(process)  # only central value, no syst. shifts 
    print 'customization done'

else:
    print 'mc'
    customizeMultiLeptonMultiJetSystematicsForMC(process)  # only central value, no syst. shifts 
    print 'customization done'




process.TFileService = cms.Service("TFileService",
									fileName = cms.string("mytree.root")
									)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )



if customize.processType == "data":
	print 'data' 
	process.p = cms.Path(process.hltHighLevel*
					# process.flashggTagSequence*
					process.flashggMultiLeptonMultiJetSystematics*
					process.analysisTree)

else : 
	print 'mc'
	process.p = cms.Path(
					# process.flashggTagSequence*
					process.flashggMultiLeptonMultiJetSystematics*
					process.analysisTree)


## set default options if needed
customize.setDefault("maxEvents", 10000)
customize.setDefault("targetLumi",1e+3)
## call the customization
customize(process)
