#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"
#include "dafne/DataFormats/interface/DiEleDiJetTag.h"
#include "dafne/DataFormats/interface/DiMuDiJetTag.h"
#include "dafne/DataFormats/interface/DiEleDiTrackTag.h"
#include "dafne/DataFormats/interface/DiMuDiTrackTag.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetTagCandidate.h"
#include "dafne/DataFormats/interface/TriLeptonsCandidate.h"
#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"


#include <vector>
#include <map>

namespace  {
    struct massiveNuDictionary {

        std::vector<pat::PackedCandidate> vec_track;
        edm::Ptr<pat::PackedCandidate> ptr_track;
        std::vector<edm::Ptr<pat::PackedCandidate> > vec_ptr_track;
        edm::Wrapper<std::vector<edm::Ptr<pat::PackedCandidate> > > wrp_vec_ptr_track;


        flashgg::DiLeptonDiJetCandidate                                       dldj;
        edm::Wrapper<flashgg::DiLeptonDiJetCandidate>                      wrp_dldj;
        std::vector<flashgg::DiLeptonDiJetCandidate>                       vec_dldj;
        edm::Wrapper<std::vector<flashgg::DiLeptonDiJetCandidate> >    wrp_vec_dldj;
        edm::Ptr<flashgg::DiLeptonDiJetCandidate>                          ptr_dldj;
        edm::Wrapper<edm::Ptr<flashgg::DiLeptonDiJetCandidate> >       wrp_ptr_dldj;
        std::vector<edm::Ptr<flashgg::DiLeptonDiJetCandidate> >        vec_ptr_dldj;
        edm::Wrapper<std::vector<edm::Ptr<flashgg::DiLeptonDiJetCandidate> > >   wrp_vec_ptr_dldj;
 
        flashgg::TriLeptonsCandidate                                       tl;
        edm::Wrapper<flashgg::TriLeptonsCandidate>                      wrp_tl;
        std::vector<flashgg::TriLeptonsCandidate>                       vec_tl;
        edm::Wrapper<std::vector<flashgg::TriLeptonsCandidate> >    wrp_vec_tl;
        edm::Ptr<flashgg::TriLeptonsCandidate>                          ptr_tl;
        edm::Wrapper<edm::Ptr<flashgg::TriLeptonsCandidate> >       wrp_ptr_tl;
        std::vector<edm::Ptr<flashgg::TriLeptonsCandidate> >        vec_ptr_tl;
        edm::Wrapper<std::vector<edm::Ptr<flashgg::TriLeptonsCandidate> > >   wrp_vec_ptr_tl;

        flashgg::MultiLeptonMultiJetCandidate                                       mlmj;
        edm::Wrapper<flashgg::MultiLeptonMultiJetCandidate>                      wrp_mlmj;
        std::vector<flashgg::MultiLeptonMultiJetCandidate>                       vec_mlmj;
        edm::Wrapper<std::vector<flashgg::MultiLeptonMultiJetCandidate> >    wrp_vec_mlmj;
        edm::Ptr<flashgg::MultiLeptonMultiJetCandidate>                          ptr_mlmj;
        edm::Wrapper<edm::Ptr<flashgg::MultiLeptonMultiJetCandidate> >       wrp_ptr_mlmj;
        std::vector<edm::Ptr<flashgg::MultiLeptonMultiJetCandidate> >        vec_ptr_mlmj;
        edm::Wrapper<std::vector<edm::Ptr<flashgg::MultiLeptonMultiJetCandidate> > >   wrp_vec_ptr_mlmj;

        flashgg::DiLeptonDiJetTagBase dldj_tagbase;
        std::vector<flashgg::DiLeptonDiJetTagBase> vec_dldj_tagbase;
        edm::Wrapper<std::vector<flashgg::DiLeptonDiJetTagBase> > wrp_vec_dldj_tagbase;
        edm::Ptr<flashgg::DiLeptonDiJetTagBase> Ptr_dldj_tagbase;
        edm::Wrapper<edm::Ptr<flashgg::DiLeptonDiJetTagBase> > wrp_ptr_dldj_tagbase;
        edm::OwnVector<flashgg::DiLeptonDiJetTagBase, edm::ClonePolicy<flashgg::DiLeptonDiJetTagBase> > ownvec_dldj_tagbase;
        edm::Wrapper<edm::OwnVector<flashgg::DiLeptonDiJetTagBase, edm::ClonePolicy<flashgg::DiLeptonDiJetTagBase> > > wrp_ownvec_dldj_tagbase;

        flashgg::DiEleDiJetTag dedj;
        std::vector<flashgg::DiEleDiJetTag> vec_dedj;
        edm::Wrapper<std::vector<flashgg::DiEleDiJetTag> > wrp_vec_dedj;

        flashgg::DiMuDiJetTag dmdj;
        std::vector<flashgg::DiMuDiJetTag> vec_dmdj;
        edm::Wrapper<std::vector<flashgg::DiMuDiJetTag> > wrp_vec_dmdj;

        flashgg::DiEleDiTrackTag dedt;
        std::vector<flashgg::DiEleDiTrackTag> vec_dedt;
        edm::Wrapper<std::vector<flashgg::DiEleDiTrackTag> > wrp_vec_dedt;

        flashgg::DiMuDiTrackTag dmdt;
        std::vector<flashgg::DiMuDiTrackTag> vec_dmdt;
        edm::Wrapper<std::vector<flashgg::DiMuDiTrackTag> > wrp_vec_dmdt;

        flashgg::DiLeptonDiJetTagCandidate                                        dldjtags;
        edm::Wrapper<flashgg::DiLeptonDiJetTagCandidate>                      wrp_dldjtags;
        std::vector<flashgg::DiLeptonDiJetTagCandidate>                       vec_dldjtags;
        edm::Wrapper<std::vector<flashgg::DiLeptonDiJetTagCandidate> >    wrp_vec_dldjtags;

    };
}
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4