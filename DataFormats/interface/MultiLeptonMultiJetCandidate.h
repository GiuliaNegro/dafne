#ifndef FLASHgg_MultiLeptonMultiJetCandidate_h
#define FLASHgg_MultiLeptonMultiJetCandidate_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/DataFormats/interface/WeightedObject.h"

#include "DataFormats/PatCandidates/interface/Electron.h"  
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h" 
#include "flashgg/DataFormats/interface/Jet.h"


using namespace std;

namespace flashgg {

	typedef flashgg::Electron Electron_t;  
	typedef edm::Ptr<flashgg::Electron> Electron_ptr; 

	typedef flashgg::Muon Muon_t;
	typedef edm::Ptr<flashgg::Muon> Muon_ptr;

	typedef flashgg::Jet Jet_t;
	typedef edm::Ptr<flashgg::Jet> Jet_ptr;

	typedef pat::PackedCandidate Track_t;
	typedef edm::Ptr<pat::PackedCandidate> Track_ptr;

	typedef reco::Vertex Vertex_t;
	typedef edm::Ptr<reco::Vertex> Vertex_ptr;

	class MultiLeptonMultiJetCandidate : public reco::LeafCandidate, public WeightedObject
	{
	public:
		enum CandidateType_t { kEEJJ, kMMJJ, kEETT, kMMTT, kEMJJ };

		MultiLeptonMultiJetCandidate();

		MultiLeptonMultiJetCandidate( Electron_ptr, Electron_ptr, Jet_ptr, Jet_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Electron_ptr, Electron_ptr, Track_ptr, Track_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Muon_ptr, Muon_ptr, Jet_ptr, Jet_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Muon_ptr, Muon_ptr, Track_ptr, Track_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Electron_ptr, Muon_ptr, Jet_ptr, Jet_ptr, Vertex_ptr );

		~MultiLeptonMultiJetCandidate();

		const CandidateType_t type() const { return type_; }

		const Vertex_ptr vtx() const { return vertex_; }
		void setVtx( Vertex_ptr val ) { vertex_ = val; }

		const vector<Electron_ptr> & electrons() const { return ptrEle_; }
		const vector<Muon_ptr>  & muons()  const { return ptrMuon_; }
		const vector<Jet_ptr>   & jets()   const { return ptrJet_; }
		const vector<Track_ptr> & tracks() const { return ptrTrack_; }

		const Electron_t *leadingEle() const; 
		const Electron_t *subLeadingEle() const; 

		const Muon_t *leadingMuon() const; 
		const Muon_t *subLeadingMuon() const;

		const Jet_t *leadingJet() const; 
		const Jet_t *subLeadingJet() const; 

		const Track_t *leadingTrack() const;
		const Track_t *subLeadingTrack() const; 

		const reco::Candidate * leadingLepton() const; 
		const reco::Candidate * subLeadingLepton() const; 

		void embedElectrons();
		vector<Electron_t> & embeddedElectrons();

		void embedMuons();
		vector<Muon_t> & embeddedMuons();

		void embedJets();
		vector<Jet_t> & embeddedJets();

		Electron_t &getLeadingElectron();
		Electron_t &getSubLeadingElectron();

		Muon_t &getLeadingMuon();
		Muon_t &getSubLeadingMuon();

		Jet_t &getLeadingJet();
		Jet_t &getSubLeadingJet();

		float sumPt() const;

		bool operator <( const MultiLeptonMultiJetCandidate &b ) const;
		bool operator >( const MultiLeptonMultiJetCandidate &b ) const;

		bool isEEJJ() const { return type_ == kEEJJ; }

		bool isEETT() const { return type_ == kEETT; }

		bool isMMJJ() const { return type_ == kMMJJ; }

		bool isMMTT() const { return type_ == kMMTT; }

		bool isEMJJ() const { return type_ == kEMJJ; }

		MultiLeptonMultiJetCandidate *clone() const { return ( new MultiLeptonMultiJetCandidate( *this ) ); }

	private:
		
		CandidateType_t type_;

		vector<Electron_ptr> ptrEle_;
		vector<Electron_t> ele_;

		vector<Muon_ptr> ptrMuon_;
		vector<Muon_t>   muon_;

		vector<Jet_ptr> ptrJet_;
		vector<Jet_t>   jet_;

		vector<Track_ptr> ptrTrack_;
		vector<Track_t>   track_;

		Vertex_ptr vertex_;
	};
}


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
