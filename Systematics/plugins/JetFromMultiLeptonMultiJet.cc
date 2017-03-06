#include "dafne/Systematics/interface/JetFromMultiLeptonMultiJetBase.h"

namespace flashgg {

    typedef JetFromMultiLeptonMultiJetBase<int> JetFromMultiLeptonMultiJet;

}

DEFINE_EDM_PLUGIN( FlashggSystematicMultiLeptonMultiJetMethodsFactory,
                   flashgg::JetFromMultiLeptonMultiJet,
                   "FlashggJetFromMultiLeptonMultiJet" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4