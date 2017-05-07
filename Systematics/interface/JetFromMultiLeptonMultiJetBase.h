#ifndef dafne_JetFromMultiLeptonMultiJetBase_h
#define dafne_JetFromMultiLeptonMultiJetBase_h

#include "dafne/Systematics/interface/BaseSystMethod_dafne.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include <memory>

namespace flashgg {


	template <class param_var>
	class JetFromMultiLeptonMultiJetBase : public BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>
	{

	public:
		JetFromMultiLeptonMultiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv );

		void applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		float makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		std::string shiftLabel( param_var ) const override;
		void eventInitialize( const edm::Event &iEvent, const edm::EventSetup & iSetup ) override;

		void setRandomEngine( CLHEP::HepRandomEngine &eng ) override
		{
			//            std::cout << " JetFromMultiLeptonMultiJetBase::setRandomEngine " << std::endl;
			BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::setRandomEngine( eng );
			jet_corr_->setRandomEngine( eng );
		}
		
	protected:
		bool debug_;

	private:
		std::unique_ptr<BaseSystMethod<flashgg::Jet, param_var> > jet_corr_;
		std::unique_ptr<BaseSystMethod<flashgg::Jet, param_var> > jet_corr2_;
	};

	template<class param_var>
	JetFromMultiLeptonMultiJetBase<param_var>::JetFromMultiLeptonMultiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv ) :
		BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::BaseSystMethod( conf, std::forward<edm::ConsumesCollector>(iC) ),
		debug_( conf.getUntrackedParameter<bool>( "Debug", false ) )
	{
		std::string JetMethodName = conf.getParameter<std::string>( "JetMethodName" );
		jet_corr_.reset( FlashggSystematicMethodsFactory<flashgg::Jet, param_var>::get()->create( JetMethodName, conf, std::forward<edm::ConsumesCollector>(iC), gv ) );
		if(conf.exists("BinList2"))  //if defined, BinList2 gives bins for sublead, lead uses BinList
			{
				edm::ParameterSet conf2;// =  conf.clone();
				
				conf2.copyFrom(conf,"JetMethodName");
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
				jet_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Jet, param_var>::get()->create( JetMethodName, conf2, std::forward<edm::ConsumesCollector>(iC),  gv ) );
				
			}
		else { //if BinList2 is not defined, use BinList for both lead and sublead photons
			jet_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Jet, param_var>::get()->create( JetMethodName, conf, std::forward<edm::ConsumesCollector>(iC),  gv ) );
		}
		this->setMakesWeight( jet_corr_->makesWeight() );
	}

	template<class param_var>
	std::string JetFromMultiLeptonMultiJetBase<param_var>::shiftLabel( param_var syst_value ) const
	{
		return jet_corr_->shiftLabel( syst_value );
	}

	template<typename param_var>
	float JetFromMultiLeptonMultiJetBase<param_var>::makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		float weight1 = 1.;
		float weight2 = 1.;

		if (y.isEEJJ() || y.isMMJJ() || y.isEMJJ()) {
			if( debug_ ) {
				std::cout << "START OF JetFromMultiLeptonMultiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingJet()->energy() << " " << y.subLeadingJet()->energy() << " "
						<< y.leadingJet()->eta() << " " << y.subLeadingJet()->eta() 
						<< std::endl;
			}

			weight1 = jet_corr_->makeWeight( *(y.leadingJet()), syst_shift );
			weight2 = jet_corr2_->makeWeight( *(y.subLeadingJet()), syst_shift );
		}

		float dieleweight = weight1*weight2;

		if( debug_ ) {
			std::cout << "END OF JetFromMultiLeptonMultiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
					<< " weight1=" << weight1 << " weight2=" << weight2 << " dieleweight=" << dieleweight << std::endl;
		}

		return dieleweight;
	}


	template<class param_var> 
	void JetFromMultiLeptonMultiJetBase<param_var>::applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		if (y.isEEJJ() || y.isMMJJ() || y.isEMJJ()) {
			if( debug_ ) {
				std::cout << "START OF JetFromMultiLeptonMultiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingJet()->energy() << " " << y.subLeadingJet()->energy() << " "
						<< y.leadingJet()->eta() << " " << y.subLeadingJet()->eta() 
						<< std::endl;
			}		

			y.embedJets();

			jet_corr_->applyCorrection( y.getLeadingJet(), syst_shift );
			jet_corr2_->applyCorrection( y.getSubLeadingJet(), syst_shift );

			if ( debug_ ) {
				std::cout << "END OF JetFromMultiLeptonMultiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingJet()->energy() << " " << y.subLeadingJet()->energy() << " "
						<< y.leadingJet()->eta() << " " << y.subLeadingJet()->eta() 
						<< std::endl;
			}

			if (y.leadingJet()->pt() < y.subLeadingJet()->pt()) y.swapJets();			
		}
	}

	template<class param_var>
	void JetFromMultiLeptonMultiJetBase<param_var>::eventInitialize( const edm::Event &ev, const edm::EventSetup & es )
	{
		if( debug_ ) {
			std::cout << "calling event initialize for both electrons " << std::endl;
		}
		jet_corr_->eventInitialize( ev, es );
		jet_corr2_->eventInitialize( ev, es );
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