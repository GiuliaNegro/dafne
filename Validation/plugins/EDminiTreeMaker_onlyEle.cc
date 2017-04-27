#include "FWCore/Framework/interface/MakerMacros.h"

#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"
#include "dafne/Validation/interface/miniTreeMaker_onlyEle.h"


typedef edm::AnalyzerWrapper<miniTreeMaker_onlyEle> EDminiTreeMaker_onlyEle;
DEFINE_FWK_MODULE(EDminiTreeMaker_onlyEle);