import FWCore.ParameterSet.Config as cms
from os import environ



def customizeMultiLeptonMultiJetSystematicsForData(process):
	#EGM scale 1D: put in central value, but omit shifts

	electronScaleBinsData = getattr(process,'electronScaleBinsData',None)
	if hasattr(process,'electronScaleBinsData'):
		print electronScaleBinsData, process.electronScaleBinsData

	createJECESource(process)
	process.jec.connect = cms.string('sqlite_file:%s/src/flashgg/Systematics/data/JEC/Summer16_23Sep2016AllV4_DATA.db' % environ['CMSSW_BASE'])
	process.jec.toGet[0].tag = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs')

	process.flashggMultiLeptonMultiJetSystematics.SystMethods.insert(0, process.JEC)

	for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
		pset.ApplyCentralValue = cms.bool(True) # Turn on central shift 
		pset.NSigmas = cms.vint32() # Do not perform up/down syst shift (1D case)

		if electronScaleBinsData != None: 
				pset.BinList = electronScaleBinsData

		if pset.Label.value().count("JEC"):
			pset.SetupUncertainties = False
			


def customizeMultiLeptonMultiJetSystematicsForMC(process):
	#EGM smearing 2D: put in central value, but omit shifts

	for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods2D:
		pset.ApplyCentralValue = cms.bool(True)  # Turn on central shift 
		pset.NSigmas = cms.PSet( firstVar = cms.vint32(), secondVar = cms.vint32() ) # Do not perform up/down syst shifts (2D case)

    # electronSmearBins = getattr(process,'electronSmearBins',None)
    # electronScaleUncertBins = getattr(process,'electronScaleUncertBins',None)
    # for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
    #     if electronSmearBins and pset.Label.value().startswith("MCSmear"):
    #         pset.BinList = electronSmearBins
    #     elif electronScaleUncertBins and pset.Label.value().count("Scale"):
    #         pset.BinList = electronScaleUncertBins

	createJECESource(process)
	# createJERESource(process)

	for isyst in [ process.JEC, process.JER ]:
		process.flashggMultiLeptonMultiJetSystematics.SystMethods.insert(0, isyst)

	for pset in process.flashggMultiLeptonMultiJetSystematics.SystMethods:
		if pset.Label.value().count("JEC") or pset.Label.value().count("JER"):
			pset.NSigmas = cms.vint32() 
			pset.ApplyCentralValue = cms.bool(True) # Turn on central shift 
			if hasattr(pset,"SetupUncertainties"):
				pset.SetupUncertainties = False



def createJECESource(process):

	process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
	# process.jetCorrectorChain = cms.Sequence(process.ak4PFCHSL1FastL2L3CorrectorChain)
	process.jetCorrectorChain = cms.Sequence(process.ak4PFCHSL1FastL2L3ResidualCorrectorChain)

	datadir = "%s/src/flashgg/Systematics/data/JEC" % environ['CMSSW_BASE']
	print "WARNING: we are reading JEC from %s so GRID jobs might not work" % datadir
	process.load("CondCore.DBCommon.CondDBCommon_cfi")
	process.load("CondCore.DBCommon.CondDBSetup_cfi")
	process.jec = cms.ESSource("PoolDBESSource",
						DBParameters = cms.PSet(
							messageLevel = cms.untracked.int32(0)
						),
						timetype = cms.string('runnumber'),
						toGet = cms.VPSet(cms.PSet(
							record = cms.string('JetCorrectionsRecord'),
							tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'),
							label  = cms.untracked.string("AK4PFchs")
						)),
						connect = cms.string('sqlite_file:%s/Summer16_23Sep2016V4_MC.db' % datadir)
			   		)                               
	process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')



def createJERESource(process):
	datadir = "%s/src/flashgg/Systematics/data/JER" % environ['CMSSW_BASE']
	print "WARNING: we are reading JER from %s so GRID jobs might not work" % datadir
	process.load('Configuration.StandardSequences.Services_cff')
	process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
	from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup

	process.jer = cms.ESSource("PoolDBESSource",
						CondDBSetup,
						toGet = cms.VPSet(
							# Resolution
							cms.PSet(
								record = cms.string('JetResolutionRcd'),
								tag    = cms.string('JR_Summer16_23Sep2016V4_MC_PtResolution_AK4PFchs'),
								label  = cms.untracked.string('AK4PFchs_pt')
							),
					
							# Scale factors
							cms.PSet(
								record = cms.string('JetResolutionScaleFactorRcd'),
								tag    = cms.string('JR_Summer16_23Sep2016V4_MC_SF_AK4PFchs'),
								label  = cms.untracked.string('AK4PFchs')
							),
						),  
						connect = cms.string('sqlite_file:%s/Summer16_23Sep2016V4_MC.db' % datadir)
					)
	process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')




