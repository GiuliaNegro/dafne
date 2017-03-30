#ifndef dafne_MuonFromMultiLeptonMultiJetBaseForTriggerSF_h
#define dafne_MuonFromMultiLeptonMultiJetBaseForTriggerSF_h

#include "dafne/Systematics/interface/BaseSystMethod_dafne.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include <memory>

namespace flashgg {


	template <class param_var>
	class MuonFromMultiLeptonMultiJetBaseForTriggerSF : public BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>
	{

	public:
		MuonFromMultiLeptonMultiJetBaseForTriggerSF( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv );

		void applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		float makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		std::string shiftLabel( param_var ) const override;
		void eventInitialize( const edm::Event &iEvent, const edm::EventSetup & iSetup ) override;

		void setRandomEngine( CLHEP::HepRandomEngine &eng ) override
		{
			//            std::cout << " MuonFromMultiLeptonMultiJetBaseForTriggerSF::setRandomEngine " << std::endl;
			BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::setRandomEngine( eng );
			muon_SF_->setRandomEngine( eng );
		}
		
	protected:
		bool debug_;

	private:
		std::unique_ptr<BaseSystMethod<flashgg::Muon, param_var> > muon_SF_;
		std::unique_ptr<BaseSystMethod<flashgg::Muon, param_var> > muon_effData_;
	};

	template<class param_var>
	MuonFromMultiLeptonMultiJetBaseForTriggerSF<param_var>::MuonFromMultiLeptonMultiJetBaseForTriggerSF( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv ) :
		BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::BaseSystMethod( conf, std::forward<edm::ConsumesCollector>(iC) ),
		debug_( conf.getUntrackedParameter<bool>( "Debug", false ) )
	{
		std::string MuonMethodName = conf.getParameter<std::string>( "MuonMethodName" );
		muon_SF_.reset( FlashggSystematicMethodsFactory<flashgg::Muon, param_var>::get()->create( MuonMethodName, conf, std::forward<edm::ConsumesCollector>(iC), gv ) );
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
				muon_effData_.reset( FlashggSystematicMethodsFactory<flashgg::Muon, param_var>::get()->create( MuonMethodName, conf2, std::forward<edm::ConsumesCollector>(iC),  gv ) );
				
			}
		else { //if BinList2 is not defined, use BinList for both lead and sublead photons
			muon_effData_.reset( FlashggSystematicMethodsFactory<flashgg::Muon, param_var>::get()->create( MuonMethodName, conf, std::forward<edm::ConsumesCollector>(iC),  gv ) );
		}
		this->setMakesWeight( muon_SF_->makesWeight() );
	}

	template<class param_var>
	std::string MuonFromMultiLeptonMultiJetBaseForTriggerSF<param_var>::shiftLabel( param_var syst_value ) const
	{
		return muon_SF_->shiftLabel( syst_value );
	}

	template<typename param_var>
	float MuonFromMultiLeptonMultiJetBaseForTriggerSF<param_var>::makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		float wSF1 = 1.;
		float wSF2 = 1.;
		float wEffData1 = 1.;
		float wEffData2 = 1.;
		float dimuonweight = 1.;

		if (y.isMMJJ() || y.isMMTT()) {
			if( debug_ ) {
				std::cout << "START OF MuonFromMultiLeptonMultiJetForTriggerSF::makeWeight M PT1 ETA1 PT2 ETA2 "
						<< y.mass() << " " << y.leadingMuon()->pt() << " " << y.leadingMuon()->eta() << " " 
						<< y.subLeadingMuon()->pt() << " " << y.subLeadingMuon()->eta() << std::endl;
			}

			wSF1 = muon_SF_->makeWeight( *(y.leadingMuon()), syst_shift );
			wSF2 = muon_SF_->makeWeight( *(y.subLeadingMuon()), syst_shift );

			wEffData1 = muon_effData_->makeWeight( *(y.leadingMuon()), syst_shift );
			wEffData2 = muon_effData_->makeWeight( *(y.subLeadingMuon()), syst_shift );

			float wEffData = wEffData1 + wEffData2 - (wEffData1*wEffData2);
			float wEffMC = (wEffData1/wSF1) + (wEffData2/wSF2) - ((wEffData1/wSF1)*(wEffData2/wSF2));

			dimuonweight = wEffData/wEffMC;

			if( debug_ ) {
				std::cout << "END OF MuonFromMultiLeptonMultiJetForTriggerSF::makeWeight M PT1 ETA1 PT2 ETA2 "
						// << " wSF1=" << wSF1 << " wSF2=" << wSF2 
						// << " wEffData1=" << wEffData1 << " wEffData2=" << wEffData2 
						<< " wEffData=" << wEffData << " wEffMC=" << wEffMC 
						<< " dimuonweight=" << dimuonweight << std::endl;
			}
		}

		if (y.isEMJJ()) {
			if( debug_ ) {
				std::cout << "START OF MuonFromMultiLeptonMultiJetForTriggerSF::makeWeight M PT1 ETA1 "
						<< y.mass() << " " << y.leadingMuon()->pt() << " " << y.leadingMuon()->eta() << std::endl;
			}
			wSF1 = muon_SF_->makeWeight( *(y.leadingMuon()), syst_shift );
			dimuonweight = wSF1;

			if( debug_ ) {
				std::cout << "END OF MuonFromMultiLeptonMultiJetForTriggerSF::makeWeight M PT1 ETA1 "
						<< " wSF1=" << wSF1 << " dimuonweight=" << dimuonweight << std::endl;
			}
		}

		return dimuonweight;
	}


	template<class param_var> 
	void MuonFromMultiLeptonMultiJetBaseForTriggerSF<param_var>::applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		// if (y.isMMJJ() || y.isMMTT()) {
		// 	if( debug_ ) {
		// 		std::cout << "START OF MuonFromMultiLeptonMultiJetForTriggerSF::applyCorrection M PT E1 E2 ETA1 ETA2 "
		// 				<< y.mass() << " " << y.sumPt() << " " 
		// 				<< y.leadingMuon()->energy() << " " << y.subLeadingMuon()->energy() << " "
		// 				<< y.leadingMuon()->eta() << " " << y.subLeadingMuon()->eta() 
		// 				<< std::endl;
		// 	}		

		// 	y.embedMuons();
		// 	///FIXME
		// 	muon_SF_->applyCorrection( y.getLeadingMuon(), syst_shift );
		// 	muon_SF_->applyCorrection( y.getSubLeadingMuon(), syst_shift );
		// 	muon_effData_->applyCorrection( y.getLeadingMuon(), syst_shift );
		// 	muon_effData_->applyCorrection( y.getSubLeadingMuon(), syst_shift );
		// 	///

		// 	if ( debug_ ) {
		// 		std::cout << "END OF MuonFromMultiLeptonMultiJetForTriggerSF::applyCorrection M PT E1 E2 ETA1 ETA2 "
		// 				<< y.mass() << " " << y.sumPt() << " " 
		// 				<< y.leadingMuon()->energy() << " " << y.subLeadingMuon()->energy() << " "
		// 				<< y.leadingMuon()->eta() << " " << y.subLeadingMuon()->eta() 
		// 				<< std::endl;
		// 	}
		// }

		// if (y.isEMJJ()) {
		// 	if( debug_ ) {
		// 		std::cout << "START OF MuonFromMultiLeptonMultiJetForTriggerSF::applyCorrection M PT E1 ETA1 "
		// 				<< y.mass() << " " << y.sumPt() << " " << y.leadingMuon()->energy() << " " 
		// 				<< y.leadingMuon()->eta() << std::endl;
		// 	}		

		// 	y.embedMuons();
		// 	muon_SF_->applyCorrection( y.getLeadingMuon(), syst_shift );

		// 	if( debug_ ) {
		// 		std::cout << "START OF MuonFromMultiLeptonMultiJetForTriggerSF::applyCorrection M PT E1 ETA1 "
		// 				<< y.mass() << " " << y.sumPt() << " " << y.leadingMuon()->energy() << " " 
		// 				<< y.leadingMuon()->eta() << std::endl;
		// 	}	
		// }
	}

	template<class param_var>
	void MuonFromMultiLeptonMultiJetBaseForTriggerSF<param_var>::eventInitialize( const edm::Event &ev, const edm::EventSetup & es )
	{
		if( debug_ ) {
			std::cout << "calling event initialize for both muons " << std::endl;
		}
		muon_SF_->eventInitialize( ev, es );
		muon_effData_->eventInitialize( ev, es );
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