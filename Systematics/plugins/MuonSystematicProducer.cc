#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/Systematics/interface/BaseSystMethod.h"
#include "flashgg/Systematics/interface/ObjectSystematicProducer.h"

namespace flashgg {

    typedef ObjectSystematicProducer<Muon, int, std::vector> MuonSystematicProducer;

}

typedef flashgg::MuonSystematicProducer FlashggMuonSystematicProducer;
DEFINE_FWK_MODULE( FlashggMuonSystematicProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4