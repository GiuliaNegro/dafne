import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt



wReejjHLTFilterGsfTrkIdVL = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*'
    ]
)

wReejjHLTFilterMW = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_DoubleEle33_CaloIdL_MW_v*'
    ]
)

wRmumujjHLTFilter = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu50_v*',       
        'HLT_TkMu50_v*'
        ]
)

wRemujjHLTFilter =  hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v*',
        'HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*'
        ]
)

tagAndProbeDoubleEleHLTFilter = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Ele27_WPTight_Gsf_v*'
    ]
)

tagAndProbeDoubleMuHLTFilter = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_IsoMu24_v*', 
        'HLT_IsoTkMu24_v*' 
        # "HLT_IsoMu27_v*"
    ]
)


wRHLTFilter_MC =  hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    throw = cms.bool(False),
    HLTPaths = wReejjHLTFilterMW.HLTPaths + wRmumujjHLTFilter.HLTPaths + wRemujjHLTFilter.HLTPaths
)

tagAndProbeHLTFilter_MC = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    throw = cms.bool(False),
    HLTPaths = tagAndProbeDoubleMuHLTFilter.HLTPaths + tagAndProbeDoubleEleHLTFilter.HLTPaths
)

wRHLTFilter_data =  hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    throw = cms.bool(False),
    HLTPaths = wReejjHLTFilterMW.HLTPaths + wRmumujjHLTFilter.HLTPaths + wRemujjHLTFilter.HLTPaths
)

tagAndProbeHLTFilter_data = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    throw = cms.bool(False),
    HLTPaths = tagAndProbeDoubleMuHLTFilter.HLTPaths + tagAndProbeDoubleEleHLTFilter.HLTPaths
)