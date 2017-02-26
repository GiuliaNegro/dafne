#include "dafne/Systematics/interface/ElectronFromMultiLeptonMultiJetBase.h"

namespace flashgg {

    typedef ElectronFromMultiLeptonMultiJetBase<std::pair<int, int> > ElectronFromMultiLeptonMultiJet2D;

}

DEFINE_EDM_PLUGIN( FlashggSystematicMultiLeptonMultiJetMethodsFactory2D,
                   flashgg::ElectronFromMultiLeptonMultiJet2D,
                   "FlashggElectronFromMultiLeptonMultiJet2D" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4