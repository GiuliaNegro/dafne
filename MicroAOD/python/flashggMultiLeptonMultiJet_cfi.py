import FWCore.ParameterSet.Config as cms

from flashgg.MicroAOD.flashggJets_cfi import maxJetCollections
# maxJetCollections = 1

flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
									JetsTag = cms.InputTag("flashggFinalJets"),
									NCollections = cms.uint32(maxJetCollections)
									)

UnpackedJetCollectionVInputTag = cms.VInputTag()
for i in range(0,maxJetCollections):
	UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))


flashggMultiLeptonMultiJet = cms.EDProducer('FlashggMultiLeptonMultiJetProducer',     
									ElectronTag=cms.InputTag('flashggSelectedElectrons'),  
									MuonTag=cms.InputTag('flashggSelectedMuons'),     
									inputTagJets= UnpackedJetCollectionVInputTag,        
									# TrackTag=cms.InputTag('packedPFCandidates'),
									VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'), 
									##Parameters  
                                    minElectronPt=cms.double(5.),
                                    maxElectronEta=cms.double(2.4),
									minMuonPt=cms.double(5.),
									maxMuonEta=cms.double(2.4),
                                    minJetPt=cms.double(7.),
                                    maxJetEta=cms.double(2.4),
                                    minTrackPt=cms.double(2000.), #to skip tracks
                                    maxTrackEta=cms.double(2.4)
                                  )