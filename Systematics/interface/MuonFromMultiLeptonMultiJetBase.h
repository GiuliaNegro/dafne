#ifndef dafne_MuonFromMultiLeptonMultiJetBase_h
#define dafne_MuonFromMultiLeptonMultiJetBase_h

#include "dafne/Systematics/interface/BaseSystMethod_dafne.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include <memory>

namespace flashgg {


	template <class param_var>
	class MuonFromMultiLeptonMultiJetBase : public BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>
	{

	public:
		MuonFromMultiLeptonMultiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv );

		void applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		float makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		std::string shiftLabel( param_var ) const override;
		void eventInitialize( const edm::Event &iEvent, const edm::EventSetup & iSetup ) override;

		void setRandomEngine( CLHEP::HepRandomEngine &eng ) override
		{
			//            std::cout << " MuonFromMultiLeptonMultiJetBase::setRandomEngine " << std::endl;
			BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::setRandomEngine( eng );
			muon_corr_->setRandomEngine( eng );
		}
		
	protected:
		bool debug_;

	private:
		std::unique_ptr<BaseSystMethod<flashgg::Muon, param_var> > muon_corr_;
		std::unique_ptr<BaseSystMethod<flashgg::Muon, param_var> > muon_corr2_;
	};

	template<class param_var>
	MuonFromMultiLeptonMultiJetBase<param_var>::MuonFromMultiLeptonMultiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv ) :
		BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::BaseSystMethod( conf, std::forward<edm::ConsumesCollector>(iC) ),
		debug_( conf.getUntrackedParameter<bool>( "Debug", false ) )
	{
		std::string MuonMethodName = conf.getParameter<std::string>( "MuonMethodName" );
		muon_corr_.reset( FlashggSystematicMethodsFactory<flashgg::Muon, param_var>::get()->create( MuonMethodName, conf, std::forward<edm::ConsumesCollector>(iC), gv ) );
		if(conf.exists("BinList2"))  //if defined, BinList2 gives bins for sublead, lead uses BinList
			{
				edm::ParameterSet conf2;// =  conf.clone();
				
				conf2.copyFrom(conf,"MuonMethodName");
				conf2.copyFrom(conf,"MethodName");
				conf2.copyFrom(conf,"Label");
				conf2.copyFrom(conf,"NSigmas");
				conf2.copyFrom(conf,"OverallRange");
				conf2.copyFrom(conf,"Debug");
				conf2.copyFrom(conf,"ApplyCentralValue");
				const auto &pset = conf.getParameterSet( "BinList2" );
				conf2.addParameter<edm::ParameterSet>("BinList", pset);
				std::string binListName = "BinList";
				conf2.insertParameterSet(true,binListName, *(conf.retrieveUnknownParameterSet("BinList2")));
				muon_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Muon, param_var>::get()->create( MuonMethodName, conf2, std::forward<edm::ConsumesCollector>(iC),  gv ) );
				
			}
		else { //if BinList2 is not defined, use BinList for both lead and sublead photons
			muon_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Muon, param_var>::get()->create( MuonMethodName, conf, std::forward<edm::ConsumesCollector>(iC),  gv ) );
		}
		this->setMakesWeight( muon_corr_->makesWeight() );
	}

	template<class param_var>
	std::string MuonFromMultiLeptonMultiJetBase<param_var>::shiftLabel( param_var syst_value ) const
	{
		return muon_corr_->shiftLabel( syst_value );
	}

	template<typename param_var>
	float MuonFromMultiLeptonMultiJetBase<param_var>::makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		float weight1 = 1.;
		float weight2 = 1.;

		if (y.isMMJJ() || y.isMMTT()) {
			if( debug_ ) {
				std::cout << "START OF MuonFromMultiLeptonMultiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingMuon()->energy() << " " << y.subLeadingMuon()->energy() << " "
						<< y.leadingMuon()->eta() << " " << y.subLeadingMuon()->eta() 
						<< std::endl;
			}

			weight1 = muon_corr_->makeWeight( *(y.leadingMuon()), syst_shift );
			weight2 = muon_corr2_->makeWeight( *(y.subLeadingMuon()), syst_shift );
		}

		if (y.isEMJJ()) {
			if( debug_ ) {
				std::cout << "START OF MuonFromMultiLeptonMultiJet::makeWeight M PT E1 ETA1 "
						<< y.mass() << " " << y.sumPt() << " " << y.leadingMuon()->energy() << " " 
						<< y.leadingMuon()->eta() << std::endl;
			}
			weight1 = muon_corr_->makeWeight( *(y.leadingMuon()), syst_shift );
		}

		float dimuonweight = weight1*weight2;

		if( debug_ ) {
			std::cout << "END OF MuonFromMultiLeptonMultiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
					<< " weight1=" << weight1 << " weight2=" << weight2 << " dimuonweight=" << dimuonweight << std::endl;
		}

		return dimuonweight;
	}


	template<class param_var> 
	void MuonFromMultiLeptonMultiJetBase<param_var>::applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		if (y.isMMJJ() || y.isMMTT()) {
			if( debug_ ) {
				std::cout << "START OF MuonFromMultiLeptonMultiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingMuon()->energy() << " " << y.subLeadingMuon()->energy() << " "
						<< y.leadingMuon()->eta() << " " << y.subLeadingMuon()->eta() 
						<< std::endl;
			}		

			y.embedMuons();
			muon_corr_->applyCorrection( y.getLeadingMuon(), syst_shift );
			muon_corr2_->applyCorrection( y.getSubLeadingMuon(), syst_shift );

			if ( debug_ ) {
				std::cout << "END OF MuonFromMultiLeptonMultiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingMuon()->energy() << " " << y.subLeadingMuon()->energy() << " "
						<< y.leadingMuon()->eta() << " " << y.subLeadingMuon()->eta() 
						<< std::endl;
			}
			
			if (y.leadingMuon()->pt() < y.subLeadingMuon()->pt()) y.swapMuons();
		}

		if (y.isEMJJ()) {
			if( debug_ ) {
				std::cout << "START OF MuonFromMultiLeptonMultiJet::applyCorrection M PT E1 ETA1 "
						<< y.mass() << " " << y.sumPt() << " " << y.leadingMuon()->energy() << " " 
						<< y.leadingMuon()->eta() << std::endl;
			}		

			y.embedMuons();
			muon_corr_->applyCorrection( y.getLeadingMuon(), syst_shift );

			if( debug_ ) {
				std::cout << "START OF MuonFromMultiLeptonMultiJet::applyCorrection M PT E1 ETA1 "
						<< y.mass() << " " << y.sumPt() << " " << y.leadingMuon()->energy() << " " 
						<< y.leadingMuon()->eta() << std::endl;
			}	
		}
	}

	template<class param_var>
	void MuonFromMultiLeptonMultiJetBase<param_var>::eventInitialize( const edm::Event &ev, const edm::EventSetup & es )
	{
		if( debug_ ) {
			std::cout << "calling event initialize for both muons " << std::endl;
		}
		muon_corr_->eventInitialize( ev, es );
		muon_corr2_->eventInitialize( ev, es );
	}
}

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4