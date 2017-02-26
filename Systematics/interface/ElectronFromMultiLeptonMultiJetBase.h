#ifndef dafne_ElectronFromMultiLeptonMultiJetBase_h
#define dafne_ElectronFromMultiLeptonMultiJetBase_h

#include "dafne/Systematics/interface/BaseSystMethod_dafne.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include <memory>

namespace flashgg {


	template <class param_var>
	class ElectronFromMultiLeptonMultiJetBase : public BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>
	{

	public:
		ElectronFromMultiLeptonMultiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv );

		void applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		float makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift ) override;
		std::string shiftLabel( param_var ) const override;
		void eventInitialize( const edm::Event &iEvent, const edm::EventSetup & iSetup ) override;

		void setRandomEngine( CLHEP::HepRandomEngine &eng ) override
		{
			//            std::cout << " ElectronFromMultiLeptonMultiJetBase::setRandomEngine " << std::endl;
			BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::setRandomEngine( eng );
			electron_corr_->setRandomEngine( eng );
		}
		
	protected:
		bool debug_;

	private:
		std::unique_ptr<BaseSystMethod<flashgg::Electron, param_var> > electron_corr_;
		std::unique_ptr<BaseSystMethod<flashgg::Electron, param_var> > electron_corr2_;
	};

	template<class param_var>
	ElectronFromMultiLeptonMultiJetBase<param_var>::ElectronFromMultiLeptonMultiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv ) :
		BaseSystMethod<MultiLeptonMultiJetCandidate, param_var>::BaseSystMethod( conf, std::forward<edm::ConsumesCollector>(iC) ),
		debug_( conf.getUntrackedParameter<bool>( "Debug", false ) )
	{
		std::string electronMethodName = conf.getParameter<std::string>( "ElectronMethodName" );
		electron_corr_.reset( FlashggSystematicMethodsFactory<flashgg::Electron, param_var>::get()->create( electronMethodName, conf, std::forward<edm::ConsumesCollector>(iC), gv ) );
		if(conf.exists("BinList2"))  //if defined, BinList2 gives bins for sublead, lead uses BinList
			{
				edm::ParameterSet conf2;// =  conf.clone();
				
				conf2.copyFrom(conf,"ElectronMethodName");
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
				electron_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Electron, param_var>::get()->create( electronMethodName, conf2, std::forward<edm::ConsumesCollector>(iC),  gv ) );
				
			}
		else { //if BinList2 is not defined, use BinList for both lead and sublead photons
			electron_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Electron, param_var>::get()->create( electronMethodName, conf, std::forward<edm::ConsumesCollector>(iC),  gv ) );
		}
		this->setMakesWeight( electron_corr_->makesWeight() );
	}

	template<class param_var>
	std::string ElectronFromMultiLeptonMultiJetBase<param_var>::shiftLabel( param_var syst_value ) const
	{
		return electron_corr_->shiftLabel( syst_value );
	}

	template<typename param_var>
	float ElectronFromMultiLeptonMultiJetBase<param_var>::makeWeight( const MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		float weight1 = 1.;
		float weight2 = 1.;

		if (y.isEEJJ() || y.isEETT()) {
			if( debug_ ) {
				std::cout << "START OF ElectronFromMultiLeptonMultiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingEle()->energy() << " " << y.subLeadingEle()->energy() << " "
						<< y.leadingEle()->eta() << " " << y.subLeadingEle()->eta() 
						<< std::endl;
			}

			weight1 = electron_corr_->makeWeight( *(y.leadingEle()), syst_shift );
			weight2 = electron_corr2_->makeWeight( *(y.subLeadingEle()), syst_shift );
		}

		if (y.isEMJJ()) {
			if( debug_ ) {
				std::cout << "START OF ElectronFromMultiLeptonMultiJet::makeWeight M PT E1 ETA1 "
						<< y.mass() << " " << y.sumPt() << " " << y.leadingEle()->energy() << " " 
						<< y.leadingEle()->eta() << std::endl;
			}
			weight1 = electron_corr_->makeWeight( *(y.leadingEle()), syst_shift );
		}

		float dieleweight = weight1*weight2;

		if( debug_ ) {
			std::cout << "END OF ElectronFromMultiLeptonMultiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
					<< " weight1=" << weight1 << " weight2=" << weight2 << " dieleweight=" << dieleweight << std::endl;
		}

		return dieleweight;
	}


	template<class param_var> 
	void ElectronFromMultiLeptonMultiJetBase<param_var>::applyCorrection( MultiLeptonMultiJetCandidate &y, param_var syst_shift )
	{
		if (y.isEEJJ() || y.isEETT()) {
			if( debug_ ) {
				std::cout << "START OF ElectronFromMultiLeptonMultiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingEle()->energy() << " " << y.subLeadingEle()->energy() << " "
						<< y.leadingEle()->eta() << " " << y.subLeadingEle()->eta() 
						<< std::endl;
			}		

			y.embedElectrons();
			electron_corr_->applyCorrection( y.getLeadingElectron(), syst_shift );
			electron_corr2_->applyCorrection( y.getSubLeadingElectron(), syst_shift );

			if ( debug_ ) {
				std::cout << "END OF ElectronFromMultiLeptonMultiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
						<< y.mass() << " " << y.sumPt() << " " 
						<< y.leadingEle()->energy() << " " << y.subLeadingEle()->energy() << " "
						<< y.leadingEle()->eta() << " " << y.subLeadingEle()->eta() 
						<< std::endl;
			}
		}

		if (y.isEMJJ()) {
			if( debug_ ) {
				std::cout << "START OF ElectronFromMultiLeptonMultiJet::applyCorrection M PT E1 ETA1 "
						<< y.mass() << " " << y.sumPt() << " " << y.leadingEle()->energy() << " " 
						<< y.leadingEle()->eta() << std::endl;
			}		

			y.embedElectrons();
			electron_corr_->applyCorrection( y.getLeadingElectron(), syst_shift );

			if( debug_ ) {
				std::cout << "START OF ElectronFromMultiLeptonMultiJet::applyCorrection M PT E1 ETA1 "
						<< y.mass() << " " << y.sumPt() << " " << y.leadingEle()->energy() << " " 
						<< y.leadingEle()->eta() << std::endl;
			}	
		}
	}

	template<class param_var>
	void ElectronFromMultiLeptonMultiJetBase<param_var>::eventInitialize( const edm::Event &ev, const edm::EventSetup & es )
	{
		if( debug_ ) {
			std::cout << "calling event initialize for both electrons " << std::endl;
		}
		electron_corr_->eventInitialize( ev, es );
		electron_corr2_->eventInitialize( ev, es );
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