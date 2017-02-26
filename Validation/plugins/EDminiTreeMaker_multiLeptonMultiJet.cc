#include "FWCore/Framework/interface/MakerMacros.h"

#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"
#include "dafne/Validation/interface/miniTreeMaker_multiLeptonMultiJet.h"


typedef edm::AnalyzerWrapper<miniTreeMaker_multiLeptonMultiJet> EDminiTreeMaker_multiLeptonMultiJet;
DEFINE_FWK_MODULE(EDminiTreeMaker_multiLeptonMultiJet);