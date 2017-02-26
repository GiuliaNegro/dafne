#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "dafne/Systematics/interface/BaseSystMethod_dafne.h"
#include "dafne/Systematics/interface/ObjectSystematicProducer_dafne.h"

namespace flashgg {

    typedef ObjectSystematicProducer_dafne<MultiLeptonMultiJetCandidate, int, std::vector> MultiLeptonMultiJetSystematicProducer;

}

typedef flashgg::MultiLeptonMultiJetSystematicProducer FlashggMultiLeptonMultiJetSystematicProducer;
DEFINE_FWK_MODULE( FlashggMultiLeptonMultiJetSystematicProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4