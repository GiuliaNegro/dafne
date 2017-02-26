#include "dafne/Systematics/interface/ElectronFromMultiLeptonMultiJetBase.h"

namespace flashgg {

    typedef ElectronFromMultiLeptonMultiJetBase<int> ElectronFromMultiLeptonMultiJet;

}

DEFINE_EDM_PLUGIN( FlashggSystematicMultiLeptonMultiJetMethodsFactory,
                   flashgg::ElectronFromMultiLeptonMultiJet,
                   "FlashggElectronFromMultiLeptonMultiJet" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4