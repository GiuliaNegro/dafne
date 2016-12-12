import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

process = cms.Process("FlashggAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
		process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v13')
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
		process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_v11')
else:
		raise Exception,"The default setup does not support releases other than 76X and 80X"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )


## input file
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
				#data
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016/dafne/DoubleEG/cmsWR2016-dafne-v0-Run2016B-PromptReco-v2/161027_122605/0000/dafneMicroAOD_1.root"

				# mc
				"root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/cmsWR2016/dafne/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/cmsWR2016-dafne-v0-gnegro-WR-3200_ToLNu-1600_ToEEJJ_miniAOD_13TeV-2016-b59cb78551aff289588aa5c69db4a3a1/161027_124028/0000/dafneMicroAOD_1.root"
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

from dafne.MicroAOD.flashggDiLeptonDiJet_cfi import flashggDiLeptonDiJet

## analyzer
process.analysisTree = cms.EDAnalyzer('EDminiTreeMaker',
										genParticleTag = cms.InputTag( "flashggPrunedGenParticles" ),  
										generatorInfo = cms.InputTag('generator'),  																			
										PileUpTag = cms.InputTag('slimmedAddPileupInfo'),
										VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
										DiLeptonDiJetTag=cms.InputTag('flashggDiLeptonDiJet'), 
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


process.TFileService = cms.Service("TFileService",
									fileName = cms.string("mytree.root")
									)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )



if customize.processType == "data":
	print 'data' 
	process.p = cms.Path(process.hltHighLevel*
					# process.flashggTagSequence*
					process.analysisTree)

else : 
	print 'mc'
	process.p = cms.Path(
					# process.flashggTagSequence*
					process.analysisTree)


## set default options if needed
customize.setDefault("maxEvents",1000)
customize.setDefault("targetLumi",1e+3)
## call the customization
customize(process)
