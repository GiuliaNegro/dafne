import FWCore.ParameterSet.Config as cms
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag


flashggMultiLeptonMultiJet = cms.EDProducer('FlashggMultiLeptonMultiJetProducer',     
									ElectronTag=cms.InputTag('flashggSelectedElectrons'),  
									MuonTag=cms.InputTag('flashggSelectedMuons'),     
									inputTagJets= UnpackedJetCollectionVInputTag,        
									TrackTag=cms.InputTag('packedPFCandidates'),
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