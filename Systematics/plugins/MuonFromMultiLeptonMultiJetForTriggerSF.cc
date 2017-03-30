#include "dafne/Systematics/interface/MuonFromMultiLeptonMultiJetBaseForTriggerSF.h"

namespace flashgg {

    typedef MuonFromMultiLeptonMultiJetBaseForTriggerSF<int> MuonFromMultiLeptonMultiJetForTriggerSF;

}

DEFINE_EDM_PLUGIN( FlashggSystematicMultiLeptonMultiJetMethodsFactory,
                   flashgg::MuonFromMultiLeptonMultiJetForTriggerSF,
                   "FlashggMuonFromMultiLeptonMultiJetForTriggerSF" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4