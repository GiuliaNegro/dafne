#ifndef dafne_BaseSystMethod_h
#define dafne_BaseSystMethod_h

#include "flashgg/Systematics/interface/BaseSystMethod.h"
#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"


typedef FlashggSystematicMethodsFactory<flashgg::MultiLeptonMultiJetCandidate, int> FlashggSystematicMultiLeptonMultiJetMethodsFactory;
typedef FlashggSystematicMethodsFactory<flashgg::MultiLeptonMultiJetCandidate, std::pair<int, int> > FlashggSystematicMultiLeptonMultiJetMethodsFactory2D;



#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
