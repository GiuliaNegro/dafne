#include "dafne/Systematics/interface/MuonFromMultiLeptonMultiJetBase.h"

namespace flashgg {

    typedef MuonFromMultiLeptonMultiJetBase<int> MuonFromMultiLeptonMultiJet;

}

DEFINE_EDM_PLUGIN( FlashggSystematicMultiLeptonMultiJetMethodsFactory,
                   flashgg::MuonFromMultiLeptonMultiJet,
                   "FlashggMuonFromMultiLeptonMultiJet" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4