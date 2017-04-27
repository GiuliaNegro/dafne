#include "dafne/Validation/interface/miniTreeMaker_multiLeptonMultiJet.h"


miniTreeMaker_multiLeptonMultiJet::miniTreeMaker_multiLeptonMultiJet( const ParameterSet &iConfig, TFileDirectory& fs, ConsumesCollector && cc ):
	BasicAnalyzer::BasicAnalyzer(iConfig, fs),
	genParticleToken_( cc.consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "genParticleTag" ) ) ),
	genInfoToken_(cc.consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "generatorInfo" ) ) ),
	PileUpToken_(cc.consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
	vertexToken_( cc.consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
	MultiLeptonMultiJetToken_( cc.consumes<View<flashgg::MultiLeptonMultiJetCandidate> >( iConfig.getParameter<InputTag> ( "MultiLeptonMultiJetTag" ) ) ),
	genJetToken_( cc.consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
	triggerBitsToken_( cc.consumes<TriggerResults>( iConfig.getParameter<InputTag>( "triggerBits" ) ) ),
	rhoToken_(cc.consumes<double>(iConfig.getParameter <InputTag>("rhoFixedGridCollection" ) ) )
{
	bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
	lumiWeight_ = iConfig.getUntrackedParameter<double>( "lumiWeight", 1. ); //pb                                                                                                                              
  
	globalVarsDumper_ = new GlobalVariablesDumper( iConfig.getParameter<ParameterSet>( "globalVariables" ), forward<ConsumesCollector>(cc) );
  
	eventTree = fs.make<TTree>( "event", "event" );
}



miniTreeMaker_multiLeptonMultiJet::~miniTreeMaker_multiLeptonMultiJet()
{
}



// ******************************************************************************************
void miniTreeMaker_multiLeptonMultiJet::beginJob()
{
	ngen = 0;
	ngenPre = 0;
	ndldj = 0;
	npre = 0;

  // per-event tree
	eventTree->Branch( "run", &evInfo.run, "run/I" );
	eventTree->Branch( "event", &evInfo.event, "event/I" );
	eventTree->Branch( "lumi", &evInfo.lumi, "lumi/I" );
	eventTree->Branch( "rho", &evInfo.rho, "rho/F" );
	eventTree->Branch( "weight", &evInfo.weight, "weight/F" );
	eventTree->Branch( "puweight", &evInfo.puweight,"puweight/F");
	eventTree->Branch( "nvtx", &evInfo.nvtx, "nvtx/I" );
	eventTree->Branch( "npu", &evInfo.npu, "npu/I" );
	eventTree->Branch( "passEEJJhlt", &evInfo.passEEJJhlt, "passEEJJhlt/I" );
	eventTree->Branch( "passMMJJhlt", &evInfo.passMMJJhlt, "passMMJJhlt/I" );
	eventTree->Branch( "passEMJJhlt", &evInfo.passEMJJhlt, "passEMJJhlt/I" );
	eventTree->Branch( "passTandPEEhlt", &evInfo.passTandPEEhlt, "passTandPEEhlt/I" );
	eventTree->Branch( "passTandPMMhlt", &evInfo.passTandPMMhlt, "passTandPMMhlt/I" );

	eventTree->Branch( "vtx_x", &evInfo.vtx_x );
	eventTree->Branch( "vtx_y", &evInfo.vtx_y );
	eventTree->Branch( "vtx_z", &evInfo.vtx_z );

	eventTree->Branch( "isEEJJ", &evInfo.isEEJJ ); 
	eventTree->Branch( "isEETT", &evInfo.isEETT ); 
	eventTree->Branch( "isMMJJ", &evInfo.isMMJJ ); 
	eventTree->Branch( "isMMTT", &evInfo.isMMTT ); 
	eventTree->Branch( "isEMJJ", &evInfo.isEMJJ ); 

	eventTree->Branch( "isSignalRegion", &evInfo.isSignalRegion );
	eventTree->Branch( "isLowMllCR", &evInfo.isLowMllCR );
	eventTree->Branch( "isLowMlljjCR", &evInfo.isLowMlljjCR );

	eventTree->Branch( "isBB", &evInfo.isBB );
	eventTree->Branch( "isEE", &evInfo.isEE );
	eventTree->Branch( "isEB", &evInfo.isEB );

	eventTree->Branch( "passPreselections", &evInfo.passPreselections );

	eventTree->Branch( "leadingLepton_e", &evInfo.leadingLepton_e ); 
	eventTree->Branch( "leadingLepton_pt", &evInfo.leadingLepton_pt ); 
	eventTree->Branch( "leadingLepton_eta", &evInfo.leadingLepton_eta ); 
	eventTree->Branch( "leadingLepton_phi", &evInfo.leadingLepton_phi ); 
	eventTree->Branch( "leadingLepton_charge", &evInfo.leadingLepton_charge ); 

	eventTree->Branch( "subLeadingLepton_e", &evInfo.subLeadingLepton_e ); 
	eventTree->Branch( "subLeadingLepton_pt", &evInfo.subLeadingLepton_pt ); 
	eventTree->Branch( "subLeadingLepton_eta", &evInfo.subLeadingLepton_eta ); 
	eventTree->Branch( "subLeadingLepton_phi", &evInfo.subLeadingLepton_phi ); 
	eventTree->Branch( "subLeadingLepton_charge", &evInfo.subLeadingLepton_charge ); 

	eventTree->Branch( "leadingJet_e", &evInfo.leadingJet_e );
	eventTree->Branch( "leadingJet_pt", &evInfo.leadingJet_pt ); 
	eventTree->Branch( "leadingJet_eta", &evInfo.leadingJet_eta );
	eventTree->Branch( "leadingJet_phi", &evInfo.leadingJet_phi ); 
	eventTree->Branch( "leadingJet_isMatchedToGen", &evInfo.leadingJet_isMatchedToGen );
	eventTree->Branch( "leadingJet_isThight", &evInfo.leadingJet_isThight );

	eventTree->Branch( "subLeadingJet_e", &evInfo.subLeadingJet_e ); 
	eventTree->Branch( "subLeadingJet_pt", &evInfo.subLeadingJet_pt ); 
	eventTree->Branch( "subLeadingJet_eta", &evInfo.subLeadingJet_eta ); 
	eventTree->Branch( "subLeadingJet_phi", &evInfo.subLeadingJet_phi ); 
	eventTree->Branch( "subLeadingJet_isMatchedToGen", &evInfo.subLeadingJet_isMatchedToGen );
	eventTree->Branch( "subLeadingJet_isThight", &evInfo.subLeadingJet_isThight );

	eventTree->Branch( "dRLeadLeptonLeadJet", &evInfo.dRLeadLeptonLeadJet );
	eventTree->Branch( "dRLeadLeptonSubLeadJet", &evInfo.dRLeadLeptonSubLeadJet );
	eventTree->Branch( "dRSubLeadLeptonLeadJet", &evInfo.dRSubLeadLeptonLeadJet );
	eventTree->Branch( "dRSubLeadLeptonSubLeadJet", &evInfo.dRSubLeadLeptonSubLeadJet );

	eventTree->Branch( "multiLeptonMultiJet_sumPt", &evInfo.multiLeptonMultiJet_sumPt );
	eventTree->Branch( "multiLeptonMultiJet_invMass", &evInfo.multiLeptonMultiJet_invMass ); 
	eventTree->Branch( "diLepton_invMass", &evInfo.diLepton_invMass ); 
	eventTree->Branch( "diLepton_pt", &evInfo.diLepton_pt ); 
	eventTree->Branch( "diJet_invMass", &evInfo.diJet_invMass ); 
	eventTree->Branch( "diJetLeadingLepton_invMass", &evInfo.diJetLeadingLepton_invMass ); 
	eventTree->Branch( "diJetSubLeadingLepton_invMass", &evInfo.diJetSubLeadingLepton_invMass ); 

	eventTree->Branch( "leadingEle_passHEEPId", &evInfo.leadingEle_passHEEPId );
	eventTree->Branch( "leadingEle_HEEPBitMapValues", &evInfo.leadingEle_HEEPBitMapValues );
	eventTree->Branch( "leadingEle_passTightId", &evInfo.leadingEle_passTightId );
	eventTree->Branch( "leadingEle_passMediumId", &evInfo.leadingEle_passMediumId );
	eventTree->Branch( "leadingEle_passLooseId", &evInfo.leadingEle_passLooseId );
	eventTree->Branch( "leadingEle_passVetoId", &evInfo.leadingEle_passVetoId );
	eventTree->Branch( "leadingEle_passMVATightId", &evInfo.leadingEle_passMVATightId );
	eventTree->Branch( "leadingEle_passMVAMediumId", &evInfo.leadingEle_passMVAMediumId );
	eventTree->Branch( "leadingEle_idmva", &evInfo.leadingEle_idmva);
	eventTree->Branch( "leadingEle_iso", &evInfo.leadingEle_iso);
	eventTree->Branch( "leadingEle_isMatchedToGen", &evInfo.leadingEle_isMatchedToGen);
	eventTree->Branch( "leadingEle_etaSC", &evInfo.leadingEle_etaSC );
	eventTree->Branch( "leadingEle_isEcalDriven", &evInfo.leadingEle_isEcalDriven );
	eventTree->Branch( "leadingEle_dEtaIn", &evInfo.leadingEle_dEtaIn );
	eventTree->Branch( "leadingEle_dPhiIn", &evInfo.leadingEle_dPhiIn );
	eventTree->Branch( "leadingEle_hOverE", &evInfo.leadingEle_hOverE );
	eventTree->Branch( "leadingEle_full5x5_r9", &evInfo.leadingEle_full5x5_r9 );
	eventTree->Branch( "leadingEle_full5x5_sigmaIetaIeta", &evInfo.leadingEle_full5x5_sigmaIetaIeta );
	eventTree->Branch( "leadingEle_full5x5_E5x5", &evInfo.leadingEle_full5x5_E5x5 );
	eventTree->Branch( "leadingEle_full5x5_E1x5", &evInfo.leadingEle_full5x5_E1x5 );
	eventTree->Branch( "leadingEle_full5x5_E2x5", &evInfo.leadingEle_full5x5_E2x5 );
	eventTree->Branch( "leadingEle_EmHadDepth1Iso", &evInfo.leadingEle_EmHadDepth1Iso );
	eventTree->Branch( "leadingEle_ptTracksIso", &evInfo.leadingEle_ptTracksIso );
	eventTree->Branch( "leadingEle_innerLayerLostHits", &evInfo.leadingEle_innerLayerLostHits );
	eventTree->Branch( "leadingEle_dxy", &evInfo.leadingEle_dxy );
	eventTree->Branch( "leadingEle_dz", &evInfo.leadingEle_dz);
	eventTree->Branch( "leadingEle_eOverP", &evInfo.leadingEle_eOverP );
	eventTree->Branch( "leadingEle_ecalEnergy", &evInfo.leadingEle_ecalEnergy );
	eventTree->Branch( "leadingEle_hcalOverEcal", &evInfo.leadingEle_hcalOverEcal );

	eventTree->Branch( "subLeadingEle_passHEEPId", &evInfo.subLeadingEle_passHEEPId );
	eventTree->Branch( "subLeadingEle_HEEPBitMapValues", &evInfo.subLeadingEle_HEEPBitMapValues );
	eventTree->Branch( "subLeadingEle_passTightId", &evInfo.subLeadingEle_passTightId );
	eventTree->Branch( "subLeadingEle_passMediumId", &evInfo.subLeadingEle_passMediumId );
	eventTree->Branch( "subLeadingEle_passLooseId", &evInfo.subLeadingEle_passLooseId );
	eventTree->Branch( "subLeadingEle_passVetoId", &evInfo.subLeadingEle_passVetoId );
	eventTree->Branch( "subLeadingEle_passMVATightId", &evInfo.subLeadingEle_passMVATightId );
	eventTree->Branch( "subLeadingEle_passMVAMediumId", &evInfo.subLeadingEle_passMVAMediumId );
	eventTree->Branch( "subLeadingEle_idmva", &evInfo.subLeadingEle_idmva);
	eventTree->Branch( "subLeadingEle_iso", &evInfo.subLeadingEle_iso);
	eventTree->Branch( "subLeadingEle_isMatchedToGen", &evInfo.subLeadingEle_isMatchedToGen);
	eventTree->Branch( "subLeadingEle_etaSC", &evInfo.subLeadingEle_etaSC );
	eventTree->Branch( "subLeadingEle_isEcalDriven", &evInfo.subLeadingEle_isEcalDriven );
	eventTree->Branch( "subLeadingEle_dEtaIn", &evInfo.subLeadingEle_dEtaIn );
	eventTree->Branch( "subLeadingEle_dPhiIn", &evInfo.subLeadingEle_dPhiIn );
	eventTree->Branch( "subLeadingEle_hOverE", &evInfo.subLeadingEle_hOverE );
	eventTree->Branch( "subLeadingEle_full5x5_r9", &evInfo.subLeadingEle_full5x5_r9 );
	eventTree->Branch( "subLeadingEle_full5x5_sigmaIetaIeta", &evInfo.subLeadingEle_full5x5_sigmaIetaIeta );
	eventTree->Branch( "subLeadingEle_full5x5_E5x5", &evInfo.subLeadingEle_full5x5_E5x5 );
	eventTree->Branch( "subLeadingEle_full5x5_E1x5", &evInfo.subLeadingEle_full5x5_E1x5 );
	eventTree->Branch( "subLeadingEle_full5x5_E2x5", &evInfo.subLeadingEle_full5x5_E2x5 );
	eventTree->Branch( "subLeadingEle_EmHadDepth1Iso", &evInfo.subLeadingEle_EmHadDepth1Iso );
	eventTree->Branch( "subLeadingEle_ptTracksIso", &evInfo.subLeadingEle_ptTracksIso );
	eventTree->Branch( "subLeadingEle_innerLayerLostHits", &evInfo.subLeadingEle_innerLayerLostHits );
	eventTree->Branch( "subLeadingEle_dxy", &evInfo.subLeadingEle_dxy );
	eventTree->Branch( "subLeadingEle_dz", &evInfo.subLeadingEle_dz);
	eventTree->Branch( "subLeadingEle_eOverP", &evInfo.subLeadingEle_eOverP );
	eventTree->Branch( "subLeadingEle_ecalEnergy", &evInfo.subLeadingEle_ecalEnergy );
	eventTree->Branch( "subLeadingEle_hcalOverEcal", &evInfo.subLeadingEle_hcalOverEcal );

	eventTree->Branch( "leadingMuon_iso", &evInfo.leadingMuon_iso );
	eventTree->Branch( "leadingMuon_PFiso", &evInfo.leadingMuon_PFiso );
	eventTree->Branch( "leadingMuon_isHighPt", &evInfo.leadingMuon_isHighPt );
	eventTree->Branch( "leadingMuon_isTight", &evInfo.leadingMuon_isTight );
	eventTree->Branch( "leadingMuon_isMedium", &evInfo.leadingMuon_isMedium );
	eventTree->Branch( "leadingMuon_isLoose", &evInfo.leadingMuon_isLoose );
	eventTree->Branch( "leadingMuon_isMatchedToGen", &evInfo.leadingMuon_isMatchedToGen );
	eventTree->Branch( "leadingMuon_dz", &evInfo.leadingMuon_dz );
	eventTree->Branch( "leadingMuon_dxy", &evInfo.leadingMuon_dxy );
	eventTree->Branch( "leadingMuon_RochCor", &evInfo.leadingMuon_RochCor );

	eventTree->Branch( "subLeadingMuon_iso", &evInfo.subLeadingMuon_iso );
	eventTree->Branch( "subLeadingMuon_PFiso", &evInfo.subLeadingMuon_PFiso );
	eventTree->Branch( "subLeadingMuon_isHighPt", &evInfo.subLeadingMuon_isHighPt );
	eventTree->Branch( "subLeadingMuon_isTight", &evInfo.subLeadingMuon_isTight );
	eventTree->Branch( "subLeadingMuon_isMedium", &evInfo.subLeadingMuon_isMedium );
	eventTree->Branch( "subLeadingMuon_isLoose", &evInfo.subLeadingMuon_isLoose );
	eventTree->Branch( "subLeadingMuon_isMatchedToGen", &evInfo.subLeadingMuon_isMatchedToGen );
	eventTree->Branch( "subLeadingMuon_dz", &evInfo.subLeadingMuon_dz );
	eventTree->Branch( "subLeadingMuon_dxy", &evInfo.subLeadingMuon_dxy );
	eventTree->Branch( "subLeadingMuon_RochCor", &evInfo.subLeadingMuon_RochCor );

	// cout << "inizialed branches" << endl;

}
// ******************************************************************************************




// ******************************************************************************************
// analyzer
//
void miniTreeMaker_multiLeptonMultiJet::analyze(const EventBase& evt)
{
	const Event *fullEvent = dynamic_cast<const Event *>(&evt);
	const Event &iEvent = (*fullEvent);  
  

	//-------------- access edm objects

	//--- only if MC
	Handle<View<reco::GenParticle> > genParticles;
	Handle<GenEventInfoProduct> genInfo;
	Handle<View< PileupSummaryInfo> > PileupInfos;
	Handle<View<reco::GenJet> > genJets;

	if ( !iEvent.isRealData() ) {
		// cout << "MC" << endl;
		iEvent.getByToken( genParticleToken_, genParticles );		
		iEvent.getByToken( genInfoToken_, genInfo );
		iEvent.getByToken( PileUpToken_, PileupInfos );
		iEvent.getByToken( genJetToken_, genJets );

		// for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		//    Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
		//    // cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
		// }    
		
		// for( unsigned int i = 0 ; i < genJets->size(); i++ ) {
		//    Ptr<reco::GenJet> genJet = genJets->ptrAt(i);
		// }  
	}

	Handle<View<reco::Vertex> > vertices;
	iEvent.getByToken( vertexToken_, vertices );
  
	Handle<View<flashgg::MultiLeptonMultiJetCandidate> > multiLeptonMultiJets;
	iEvent.getByToken( MultiLeptonMultiJetToken_, multiLeptonMultiJets );
	
	Handle<TriggerResults> triggerBits;
	iEvent.getByToken( triggerBitsToken_, triggerBits );
  
	Handle<double> rhoHandle;
	iEvent.getByToken( rhoToken_, rhoHandle );
	double rho = *( rhoHandle.product() );

	// cout << "Handle done" << endl;


	//-------------- initialize tree
	initEventStructure();
	// cout << "initEventStructure() done" << endl;

	
	//-------------- check if event passes HLT
	const TriggerNames &triggerNames = iEvent.triggerNames( *triggerBits );

	bool passTrigger = false;

	for( unsigned index = 0; index < triggerNames.size(); ++index ) {

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_DoubleEle33_CaloIdL_MW") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passEEJJhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Mu50") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_TkMu50")) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passMMJJhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passEMJJhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Ele27_WPTight_Gsf") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passTandPEEhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_IsoMu24") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_IsoMu27") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passTandPMMhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

	}
	// cout << "arriva a linea " << __LINE__ << endl;


	globalVarsDumper_->fill( iEvent );

	evInfo.run = globalVarsDumper_->cache().run;
	evInfo.event = globalVarsDumper_->cache().event;		
	evInfo.lumi = globalVarsDumper_->cache().lumi;
	evInfo.rho = rho;

	// -- event weight (gen weight x lumi x cross section x pu /nTotEvents) with nTotEvents=sumWeights if genWeight != 1
	float w = 1.;
	// cout << "w = " << w << endl;
	if( ! iEvent.isRealData() ) {
		w = lumiWeight_;
		// cout << "lumiWeight = " << w << endl;
		if( genInfo.isValid() ) {
			const auto &weights = genInfo->weights();
			if( ! weights.empty() ) {
				w *= weights[0];
				// cout << "genWeight*lumiWeight = " << w << endl;
			}
		}

		if( globalVarsDumper_->puReWeight() ) {
			w *= globalVarsDumper_->cache().puweight;
		}

	}
	evInfo.weight = w;
	// cout << "weight = " << w << endl;


	// -- pileup weight
	if( globalVarsDumper_->puReWeight() ) {
		evInfo.puweight = globalVarsDumper_->cache().puweight;
	}
	// cout << "puweight = " << evInfo.puweight << endl;


	// -- number of reco vertices
	evInfo.nvtx = vertices->size() ;
	// cout << vertices->size() << " nvertices " << endl;

	// -- vertices
	for (UInt_t ivtx = 0 ; ivtx < vertices->size(); ivtx++){
		Ptr<reco::Vertex> vtx = vertices->ptrAt( ivtx );
		evInfo.vtx_x.push_back(vtx->x());
		evInfo.vtx_y.push_back(vtx->y());
		evInfo.vtx_z.push_back(vtx->z());
	}


	// -- number of pileup events
	float pu = 0.; 
	if( ! iEvent.isRealData() ) {
		for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
			Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
			if( pu_bunchcrossing == 0 ) {
				pu = PileupInfos->ptrAt( PVI )->getPU_NumInteractions();
			}
		}
	}
	evInfo.npu = pu;
	// cout << "pu " << pu << endl;



	// -- multiLeptonMultiJets
	// cout << "MLMJ size " << multiLeptonMultiJets->size() << endl;
	for ( unsigned int imlmj = 0; imlmj < multiLeptonMultiJets->size(); imlmj++){
		// cout << "enter MLMJ loop" << endl;
		Ptr<flashgg::MultiLeptonMultiJetCandidate> multiLeptonMultiJet = multiLeptonMultiJets->ptrAt( imlmj );        
		// cout << "arriva a linea " << __LINE__ << endl;

		// if ( multiLeptonMultiJet->isEEJJ() ) {
		// 	cout << "e1 = " << multiLeptonMultiJet->leadingEle()->energy() << ", eta1 = " << multiLeptonMultiJet->leadingEle()->eta() << endl;
		// 	cout << "e2 = " << multiLeptonMultiJet->subLeadingEle()->energy() << ", eta2 = " << multiLeptonMultiJet->subLeadingEle()->eta() << endl;
		// }

		// if (multiLeptonMultiJet->isEEJJ()) cout << "isEEJJ" << endl;
		// if (multiLeptonMultiJet->isMMJJ()) cout << "isMMJJ" << endl;
		// if (multiLeptonMultiJet->isEETT()) cout << "isEETT" << endl;
		// if (multiLeptonMultiJet->isMMTT()) cout << "isMMTT" << endl;
		// if (multiLeptonMultiJet->isEMJJ()) cout << "isEMJJ" << endl; 
		// if (!(multiLeptonMultiJet->isEMJJ())) cout << "no EMJJ" << endl;

		// if (multiLeptonMultiJet->leadingLepton()->pt() < 35 || multiLeptonMultiJet->subLeadingLepton()->pt() < 35) continue;
		if (multiLeptonMultiJet->leadingLepton()->pt() < 30 || multiLeptonMultiJet->subLeadingLepton()->pt() < 28) continue;
		if (fabs(multiLeptonMultiJet->leadingLepton()->eta()) > 2.4 || fabs(multiLeptonMultiJet->subLeadingLepton()->eta()) > 2.4) continue;

		if (multiLeptonMultiJet->leadingJet()->pt() < 40 || multiLeptonMultiJet->subLeadingJet()->pt() < 40) continue;
		if (fabs(multiLeptonMultiJet->leadingJet()->eta()) > 2.4 || fabs(multiLeptonMultiJet->subLeadingJet()->eta()) > 2.4) continue;
		// cout << "arriva a linea " << __LINE__ << endl;

		bool leadingElePassHEEPId      = false;
		unsigned leadingEleHeepBitMap  = 0;	
		bool leadingElePassTightId     = false;
		bool leadingElePassMediumId    = false;
		bool leadingElePassLooseId     = false;
		bool leadingElePassVetoId      = false;
		bool leadingElePassMVATightId  = false;
		bool leadingElePassMVAMediumId = false;
		float leadingEleIdmva          = -999.;
		float leadingEleIso            = -999.;
		int leadingEleIsMatchedToGen   = 0;
		float leadingEleEtaSC          = -999.;
		bool leadingEleIsEcalDriven    = false;
		float leadingEleDEtaIn         = -999.;
		float leadingEleDPhiIn         = -999.;
		float leadingEleHoE            = -999.;
		float leadingEleR9             = -999.;
		float leadingEleSigmaIetaIeta  = -999.;
		float leadingEleE5x5           = -999.;
		float leadingEleE1x5           = -999.; 
		float leadingEleE2x5           = -999.;
		float leadingEleEmHadDepth1Iso = -999.;
		float leadingElePtTracksIso    = -999.;
		float leadingEleMissingHits    = -999.;
		float leadingEleDxy            = -999.;
		float leadingEleDz             = -999.;
		float leadingEleEoverP         = -999.;
		float leadingEleEcalEnergy     = -999.;
		float leadingEleHcalOverEcal   = -999.;

		bool subLeadingElePassHEEPId      = false;
		unsigned subLeadingEleHeepBitMap  = 0;	
		bool subLeadingElePassTightId     = false;
		bool subLeadingElePassMediumId    = false;
		bool subLeadingElePassLooseId     = false;
		bool subLeadingElePassVetoId      = false;
		bool subLeadingElePassMVATightId  = false;
		bool subLeadingElePassMVAMediumId = false;	
		float subLeadingEleIdmva          = -999.;
		float subLeadingEleIso            = -999.;
		int subLeadingEleIsMatchedToGen   = 0;
		float subLeadingEleEtaSC          = -999.;
		bool subLeadingEleIsEcalDriven    = false;
		float subLeadingEleDEtaIn         = -999.;
		float subLeadingEleDPhiIn         = -999.;
		float subLeadingEleHoE            = -999.;
		float subLeadingEleR9             = -999.;
		float subLeadingEleSigmaIetaIeta  = -999.;
		float subLeadingEleE5x5           = -999.;
		float subLeadingEleE1x5           = -999.; 
		float subLeadingEleE2x5           = -999.;
		float subLeadingEleEmHadDepth1Iso = -999.;
		float subLeadingElePtTracksIso    = -999.;
		float subLeadingEleMissingHits    = -999.;
		float subLeadingEleDxy            = -999.;
		float subLeadingEleDz             = -999.;
		float subLeadingEleEoverP         = -999.;
		float subLeadingEleEcalEnergy     = -999.;
		float subLeadingEleHcalOverEcal   = -999.;

		float leadingMuonRochCor      = -999.;
		float leadingMuonIso          = -999.;
		float leadingMuonPFiso        = -999.;
		bool leadingMuonIsHighPt      = false;
		bool leadingMuonIsTight       = false;
		bool leadingMuonIsMedium      = false;
		bool leadingMuonIsLoose       = false;
		int leadingMuonIsMatchedToGen = 0;
		float leadingMuonDz           = -999.;
		float leadingMuonDxy          = -999.;

		float subLeadingMuonRochCor      = -999.;
		float subLeadingMuonIso          = -999.;
		float subLeadingMuonPFiso        = -999.;
		bool subLeadingMuonIsHighPt      = false;
		bool subLeadingMuonIsTight       = false;
		bool subLeadingMuonIsMedium      = false;
		bool subLeadingMuonIsLoose       = false;
		int subLeadingMuonIsMatchedToGen = 0;
		float subLeadingMuonDz           = -999.;
		float subLeadingMuonDxy          = -999.;

		math::XYZTLorentzVector leadingLeptonCorr = multiLeptonMultiJet->leadingLepton()->p4();
		math::XYZTLorentzVector subLeadingLeptonCorr = multiLeptonMultiJet->subLeadingLepton()->p4();


		if (! iEvent.isRealData() ) { 
			ngen++;
			if ( passMultiLeptonMultiJetPreselection(multiLeptonMultiJet, leadingLeptonCorr, subLeadingLeptonCorr) ) ngenPre++; 
		} else {
			ndldj++;
			if ( passMultiLeptonMultiJetPreselection(multiLeptonMultiJet, leadingLeptonCorr, subLeadingLeptonCorr) ) npre++; 
		}
		// cout << "arriva a linea " << __LINE__ << endl;


		if (multiLeptonMultiJet->isEEJJ() || multiLeptonMultiJet->isEETT() || multiLeptonMultiJet->isEMJJ()) {
			// cout << "arriva a linea " << __LINE__ << endl;
			Ptr<reco::Vertex> best_vtx_leadEle = chooseBestVtx(vertices->ptrs(), multiLeptonMultiJet->leadingEle());
			int mcMatch_leadEle = -1;
			if( ! iEvent.isRealData() ) mcMatch_leadEle = electronMatchingToGen(multiLeptonMultiJet->leadingEle(), genParticles);

			leadingElePassHEEPId = multiLeptonMultiJet->leadingEle()->passHeepId();
			leadingEleHeepBitMap = multiLeptonMultiJet->leadingEle()->heepBitMap();
			leadingElePassTightId = multiLeptonMultiJet->leadingEle()->passTightId();
			leadingElePassMediumId = multiLeptonMultiJet->leadingEle()->passMediumId();
			leadingElePassLooseId = multiLeptonMultiJet->leadingEle()->passLooseId();
			leadingElePassVetoId = multiLeptonMultiJet->leadingEle()->passVetoId();
			leadingElePassMVATightId = multiLeptonMultiJet->leadingEle()->passMVATightId();
			leadingElePassMVAMediumId = multiLeptonMultiJet->leadingEle()->passMVAMediumId();
			leadingEleIdmva = multiLeptonMultiJet->leadingEle()->nonTrigMVA();
			leadingEleIso = electronIsolation(multiLeptonMultiJet->leadingEle(), rho);
			leadingEleIsMatchedToGen = mcMatch_leadEle;
			leadingEleEtaSC = multiLeptonMultiJet->leadingEle()->superCluster()->eta();
			leadingEleIsEcalDriven = multiLeptonMultiJet->leadingEle()->ecalDrivenSeed();
			leadingEleDEtaIn = multiLeptonMultiJet->leadingEle()->deltaEtaSuperClusterTrackAtVtx();  
			leadingEleDPhiIn = multiLeptonMultiJet->leadingEle()->deltaPhiSuperClusterTrackAtVtx();
			leadingEleHoE = multiLeptonMultiJet->leadingEle()->hadronicOverEm();
			leadingEleR9 = multiLeptonMultiJet->leadingEle()->full5x5_r9();
			leadingEleSigmaIetaIeta = multiLeptonMultiJet->leadingEle()->full5x5_sigmaIetaIeta();
			leadingEleE5x5 = multiLeptonMultiJet->leadingEle()->full5x5_e5x5();
			leadingEleE1x5 = multiLeptonMultiJet->leadingEle()->full5x5_e1x5();
			leadingEleE2x5 = multiLeptonMultiJet->leadingEle()->full5x5_e2x5Max();
			leadingEleEmHadDepth1Iso = multiLeptonMultiJet->leadingEle()->dr03EcalRecHitSumEt() + multiLeptonMultiJet->leadingEle()->dr03HcalDepth1TowerSumEt();
			leadingElePtTracksIso = multiLeptonMultiJet->leadingEle()->dr03TkSumPt();  //TO MODIFY
			leadingEleMissingHits = multiLeptonMultiJet->leadingEle()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
			leadingEleDxy = multiLeptonMultiJet->leadingEle()->gsfTrack()->dxy( best_vtx_leadEle->position() );
			leadingEleDz = multiLeptonMultiJet->leadingEle()->gsfTrack()->dz( best_vtx_leadEle->position() );
			leadingEleEoverP = multiLeptonMultiJet->leadingEle()->eSuperClusterOverP();
			// cout << "arriva a linea " << __LINE__ << endl;
		}

		if  (multiLeptonMultiJet->isEEJJ() || multiLeptonMultiJet->isEETT()) {  
			// cout << "arriva a linea " << __LINE__ << endl;
			Ptr<reco::Vertex> best_vtx_subLeadEle = chooseBestVtx(vertices->ptrs(), multiLeptonMultiJet->subLeadingEle());
			int mcMatch_subLeadEle = -1;
			if( ! iEvent.isRealData() ) mcMatch_subLeadEle = electronMatchingToGen(multiLeptonMultiJet->subLeadingEle(), genParticles);

			subLeadingElePassHEEPId = multiLeptonMultiJet->subLeadingEle()->passHeepId();
			subLeadingEleHeepBitMap = multiLeptonMultiJet->subLeadingEle()->heepBitMap();
			subLeadingElePassTightId = multiLeptonMultiJet->subLeadingEle()->passTightId();
			subLeadingElePassMediumId = multiLeptonMultiJet->subLeadingEle()->passMediumId();
			subLeadingElePassLooseId = multiLeptonMultiJet->subLeadingEle()->passLooseId();
			subLeadingElePassVetoId = multiLeptonMultiJet->subLeadingEle()->passVetoId();
			subLeadingElePassMVATightId = multiLeptonMultiJet->subLeadingEle()->passMVATightId();
			subLeadingElePassMVAMediumId = multiLeptonMultiJet->subLeadingEle()->passMVAMediumId();
			subLeadingEleIdmva = multiLeptonMultiJet->subLeadingEle()->nonTrigMVA();
			subLeadingEleIso = electronIsolation(multiLeptonMultiJet->subLeadingEle(), rho);
			subLeadingEleIsMatchedToGen = mcMatch_subLeadEle;
			subLeadingEleEtaSC = multiLeptonMultiJet->subLeadingEle()->superCluster()->eta();
			subLeadingEleIsEcalDriven = multiLeptonMultiJet->subLeadingEle()->ecalDrivenSeed();
			subLeadingEleDEtaIn = multiLeptonMultiJet->subLeadingEle()->deltaEtaSuperClusterTrackAtVtx();  
			subLeadingEleDPhiIn = multiLeptonMultiJet->subLeadingEle()->deltaPhiSuperClusterTrackAtVtx();
			subLeadingEleHoE = multiLeptonMultiJet->subLeadingEle()->hadronicOverEm();
			subLeadingEleR9 = multiLeptonMultiJet->subLeadingEle()->full5x5_r9();
			subLeadingEleSigmaIetaIeta = multiLeptonMultiJet->subLeadingEle()->full5x5_sigmaIetaIeta();
			subLeadingEleE5x5 = multiLeptonMultiJet->subLeadingEle()->full5x5_e5x5();
			subLeadingEleE1x5 = multiLeptonMultiJet->subLeadingEle()->full5x5_e1x5();
			subLeadingEleE2x5 = multiLeptonMultiJet->subLeadingEle()->full5x5_e2x5Max();
			subLeadingEleEmHadDepth1Iso = multiLeptonMultiJet->subLeadingEle()->dr03EcalRecHitSumEt() + multiLeptonMultiJet->subLeadingEle()->dr03HcalDepth1TowerSumEt();
			subLeadingElePtTracksIso = multiLeptonMultiJet->subLeadingEle()->dr03TkSumPt();  //TO MODIFY
			subLeadingEleMissingHits = multiLeptonMultiJet->subLeadingEle()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
			subLeadingEleDxy = multiLeptonMultiJet->subLeadingEle()->gsfTrack()->dxy( best_vtx_subLeadEle->position() );
			subLeadingEleDz = multiLeptonMultiJet->subLeadingEle()->gsfTrack()->dz( best_vtx_subLeadEle->position() );
			subLeadingEleEoverP = multiLeptonMultiJet->subLeadingEle()->eSuperClusterOverP();		
			// cout << "arriva a linea " << __LINE__ << endl;
		} 


		if (multiLeptonMultiJet->isMMJJ() || multiLeptonMultiJet->isMMTT() || multiLeptonMultiJet->isEMJJ()) {
			// cout << "arriva a linea " << __LINE__ << endl;
			Ptr<reco::Vertex> leadingMuonVtx = chooseBestMuonVtx(vertices->ptrs(), multiLeptonMultiJet->leadingMuon());
			int mcMatch_leadMuon =  -1;
			if( ! iEvent.isRealData() ) mcMatch_leadMuon = muonMatchingToGen(multiLeptonMultiJet->leadingMuon(), genParticles); 

			// muon ID and isolation: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
			float muIso_leadMuon = multiLeptonMultiJet->leadingMuon()->isolationR03().sumPt/multiLeptonMultiJet->leadingMuon()->pt();
			float muPFCombRelIso_leadMuon = ( multiLeptonMultiJet->leadingMuon()->pfIsolationR04().sumChargedHadronPt + max( 0.,multiLeptonMultiJet->leadingMuon()->pfIsolationR04().sumNeutralHadronEt + multiLeptonMultiJet->leadingMuon()->pfIsolationR04().sumPhotonEt - 0.5 * multiLeptonMultiJet->leadingMuon()->pfIsolationR04().sumPUPt ) ) / ( multiLeptonMultiJet->leadingMuon()->pt() );

			float dz = -999;
			float dxy = -999;
			if ( !(!multiLeptonMultiJet->leadingMuon()->innerTrack()) ) {
				dz = multiLeptonMultiJet->leadingMuon()->innerTrack()->dz( leadingMuonVtx->position() );
				dxy = multiLeptonMultiJet->leadingMuon()->innerTrack()->dxy( leadingMuonVtx->position() ); 			
			}

			leadingMuonIso = muIso_leadMuon;
			leadingMuonPFiso = muPFCombRelIso_leadMuon;
			leadingMuonIsHighPt = multiLeptonMultiJet->leadingMuon()->isHighPtMuon( *leadingMuonVtx );
			leadingMuonIsTight = multiLeptonMultiJet->leadingMuon()->isTightMuon( *leadingMuonVtx );
			leadingMuonIsMedium = multiLeptonMultiJet->leadingMuon()->isMediumMuon();
			leadingMuonIsLoose = multiLeptonMultiJet->leadingMuon()->isLooseMuon();
			leadingMuonIsMatchedToGen = mcMatch_leadMuon;
			leadingMuonDz = dz;
			leadingMuonDxy = dxy;
			leadingMuonRochCor = RochesterCorrection(multiLeptonMultiJet->leadingMuon(), genParticles, iEvent.isRealData());

			float leadingMuonPtRochCor = multiLeptonMultiJet->leadingMuon()->pt()*leadingMuonRochCor;
			math::PtEtaPhiELorentzVector lMuonRochCorVec( leadingMuonPtRochCor, multiLeptonMultiJet->leadingMuon()->eta(), multiLeptonMultiJet->leadingMuon()->phi(),  multiLeptonMultiJet->leadingMuon()->energy());
			leadingLeptonCorr = lMuonRochCorVec;
			// cout << "arriva a linea " << __LINE__ << endl;
		}

		if (multiLeptonMultiJet->isMMJJ() || multiLeptonMultiJet->isMMTT()) {
			// cout << "arriva a linea " << __LINE__ << endl;
			Ptr<reco::Vertex> subLeadingMuonVtx = chooseBestMuonVtx(vertices->ptrs(), multiLeptonMultiJet->subLeadingMuon());			
			int mcMatch_subLeadMuon =  -1;
			if( ! iEvent.isRealData() ) mcMatch_subLeadMuon = muonMatchingToGen(multiLeptonMultiJet->subLeadingMuon(), genParticles); 

			// muon ID and isolation: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
			float muIso_subLeadMuon = multiLeptonMultiJet->subLeadingMuon()->isolationR03().sumPt/multiLeptonMultiJet->subLeadingMuon()->pt();
			float muPFCombRelIso_subLeadMuon = ( multiLeptonMultiJet->subLeadingMuon()->pfIsolationR04().sumChargedHadronPt + max( 0.,multiLeptonMultiJet->subLeadingMuon()->pfIsolationR04().sumNeutralHadronEt + multiLeptonMultiJet->subLeadingMuon()->pfIsolationR04().sumPhotonEt - 0.5 * multiLeptonMultiJet->subLeadingMuon()->pfIsolationR04().sumPUPt ) ) / ( multiLeptonMultiJet->subLeadingMuon()->pt() );

			float dz = -999;
			float dxy = -999;
			if ( !(!multiLeptonMultiJet->subLeadingMuon()->innerTrack()) ) {
				dz = multiLeptonMultiJet->subLeadingMuon()->innerTrack()->dz( subLeadingMuonVtx->position() );
				dxy = multiLeptonMultiJet->subLeadingMuon()->innerTrack()->dxy( subLeadingMuonVtx->position() ); 			
			}

			subLeadingMuonIso = muIso_subLeadMuon;
			subLeadingMuonPFiso = muPFCombRelIso_subLeadMuon;
			subLeadingMuonIsHighPt = multiLeptonMultiJet->subLeadingMuon()->isHighPtMuon( *subLeadingMuonVtx );
			subLeadingMuonIsTight = multiLeptonMultiJet->subLeadingMuon()->isTightMuon( *subLeadingMuonVtx );
			subLeadingMuonIsMedium = multiLeptonMultiJet->subLeadingMuon()->isMediumMuon();
			subLeadingMuonIsLoose = multiLeptonMultiJet->subLeadingMuon()->isLooseMuon();
			subLeadingMuonIsMatchedToGen = mcMatch_subLeadMuon;
			subLeadingMuonDz = dz;
			subLeadingMuonDxy = dxy;
			subLeadingMuonRochCor = RochesterCorrection(multiLeptonMultiJet->subLeadingMuon(), genParticles, iEvent.isRealData());

			float subLeadingMuonPtRochCor = multiLeptonMultiJet->subLeadingMuon()->pt()*leadingMuonRochCor;
			math::PtEtaPhiELorentzVector slMuonRochCorVec( subLeadingMuonPtRochCor, multiLeptonMultiJet->subLeadingMuon()->eta(), multiLeptonMultiJet->subLeadingMuon()->phi(),  multiLeptonMultiJet->subLeadingMuon()->energy());
			subLeadingLeptonCorr = slMuonRochCorVec;
			// cout << "arriva a linea " << __LINE__ << endl;
		}


		// cout << "arriva a linea " << __LINE__ << endl;
		float massCorr = (leadingLeptonCorr + subLeadingLeptonCorr + multiLeptonMultiJet->leadingJet()->p4() + multiLeptonMultiJet->subLeadingJet()->p4()).mass();
		float sumPtCorr = leadingLeptonCorr.pt() + subLeadingLeptonCorr.pt() + multiLeptonMultiJet->leadingJet()->pt() + multiLeptonMultiJet->subLeadingJet()->pt();

		float diLeptonInvMass = (leadingLeptonCorr + subLeadingLeptonCorr).mass();
		float diLeptonPt = (leadingLeptonCorr + subLeadingLeptonCorr).pt();
		// cout << "diLeptonInvMass = " << diLeptonInvMass << ", diLeptonPt = " << diLeptonPt << endl;

		float diJetInvMass = (multiLeptonMultiJet->leadingJet()->p4() + multiLeptonMultiJet->subLeadingJet()->p4()).mass();
		float diJetLeadingLeptonInvMass = (multiLeptonMultiJet->leadingJet()->p4() + multiLeptonMultiJet->subLeadingJet()->p4() + leadingLeptonCorr).mass();
		float diJetSubLeadingLeptonInvMass = (multiLeptonMultiJet->leadingJet()->p4() + multiLeptonMultiJet->subLeadingJet()->p4() + subLeadingLeptonCorr).mass();

		int leadingJetIsMatchedToGen =  -1;
		if( ! iEvent.isRealData() ) leadingJetIsMatchedToGen = jetMatchingToGen(multiLeptonMultiJet->leadingJet(), genJets); 
		int subLeadingJetIsMatchedToGen =  -1;
		if( ! iEvent.isRealData() ) subLeadingJetIsMatchedToGen = jetMatchingToGen(multiLeptonMultiJet->subLeadingJet(), genJets); 

		bool leadingJetIsTight = (multiLeptonMultiJet->leadingJet()->neutralHadronEnergyFraction()<0.90 && multiLeptonMultiJet->leadingJet()->neutralEmEnergyFraction()<0.9 && (multiLeptonMultiJet->leadingJet()->chargedMultiplicity()+multiLeptonMultiJet->leadingJet()->neutralMultiplicity())>1 && multiLeptonMultiJet->leadingJet()->muonEnergyFraction()<0.8) && ((fabs(multiLeptonMultiJet->leadingJet()->eta())<=2.4 && multiLeptonMultiJet->leadingJet()->chargedHadronEnergyFraction()>0 && multiLeptonMultiJet->leadingJet()->chargedMultiplicity()>0 && multiLeptonMultiJet->leadingJet()->chargedEmEnergyFraction()<0.90) || fabs(multiLeptonMultiJet->leadingJet()->eta())>2.4);
		bool subLeadingJetIsTight = (multiLeptonMultiJet->subLeadingJet()->neutralHadronEnergyFraction()<0.90 && multiLeptonMultiJet->subLeadingJet()->neutralEmEnergyFraction()<0.9 && (multiLeptonMultiJet->subLeadingJet()->chargedMultiplicity()+multiLeptonMultiJet->subLeadingJet()->neutralMultiplicity())>1 && multiLeptonMultiJet->subLeadingJet()->muonEnergyFraction()<0.8) && ((fabs(multiLeptonMultiJet->subLeadingJet()->eta())<=2.4 && multiLeptonMultiJet->subLeadingJet()->chargedHadronEnergyFraction()>0 && multiLeptonMultiJet->subLeadingJet()->chargedMultiplicity()>0 && multiLeptonMultiJet->subLeadingJet()->chargedEmEnergyFraction()<0.90) || fabs(multiLeptonMultiJet->subLeadingJet()->eta())>2.4);


		evInfo.isEEJJ.push_back(multiLeptonMultiJet->isEEJJ());  
		evInfo.isEETT.push_back(multiLeptonMultiJet->isEETT());
		evInfo.isMMJJ.push_back(multiLeptonMultiJet->isMMJJ());
		evInfo.isMMTT.push_back(multiLeptonMultiJet->isMMTT());
		evInfo.isEMJJ.push_back(multiLeptonMultiJet->isEMJJ());  

		evInfo.isSignalRegion.push_back( isSignalRegion(massCorr, diLeptonInvMass) );   
		evInfo.isLowMllCR.push_back( isLowMllCR(diLeptonInvMass) );   
		evInfo.isLowMlljjCR.push_back( isLowMlljjCR(massCorr, diLeptonInvMass) );   

		evInfo.isBB.push_back( isBB(multiLeptonMultiJet) );
		evInfo.isEE.push_back( isEE(multiLeptonMultiJet) );
		evInfo.isEB.push_back( isEB(multiLeptonMultiJet) );

		evInfo.passPreselections.push_back( passMultiLeptonMultiJetPreselection(multiLeptonMultiJet, leadingLeptonCorr, subLeadingLeptonCorr) );  

		evInfo.leadingLepton_e.push_back(multiLeptonMultiJet->leadingLepton()->energy());  
		evInfo.leadingLepton_pt.push_back(leadingLeptonCorr.pt());   
		evInfo.leadingLepton_eta.push_back(multiLeptonMultiJet->leadingLepton()->eta());  
		evInfo.leadingLepton_phi.push_back(multiLeptonMultiJet->leadingLepton()->phi());  
		evInfo.leadingLepton_charge.push_back(multiLeptonMultiJet->leadingLepton()->charge()); 

		evInfo.subLeadingLepton_e.push_back(multiLeptonMultiJet->subLeadingLepton()->energy());  
		evInfo.subLeadingLepton_pt.push_back(subLeadingLeptonCorr.pt());   
		evInfo.subLeadingLepton_eta.push_back(multiLeptonMultiJet->subLeadingLepton()->eta());   
		evInfo.subLeadingLepton_phi.push_back(multiLeptonMultiJet->subLeadingLepton()->phi());   
		evInfo.subLeadingLepton_charge.push_back(multiLeptonMultiJet->subLeadingLepton()->charge());  

		evInfo.leadingJet_e.push_back(multiLeptonMultiJet->leadingJet()->energy());
		evInfo.leadingJet_pt.push_back(multiLeptonMultiJet->leadingJet()->pt());
		evInfo.leadingJet_eta.push_back(multiLeptonMultiJet->leadingJet()->eta());
		evInfo.leadingJet_phi.push_back(multiLeptonMultiJet->leadingJet()->phi());
		evInfo.leadingJet_isMatchedToGen.push_back(leadingJetIsMatchedToGen);
		evInfo.leadingJet_isThight.push_back(leadingJetIsTight);

		evInfo.subLeadingJet_e.push_back(multiLeptonMultiJet->subLeadingJet()->energy());
		evInfo.subLeadingJet_pt.push_back(multiLeptonMultiJet->subLeadingJet()->pt());
		evInfo.subLeadingJet_eta.push_back(multiLeptonMultiJet->subLeadingJet()->eta());
		evInfo.subLeadingJet_phi.push_back(multiLeptonMultiJet->subLeadingJet()->phi());
		evInfo.subLeadingJet_isMatchedToGen.push_back(subLeadingJetIsMatchedToGen);
		evInfo.subLeadingJet_isThight.push_back(subLeadingJetIsTight);

		evInfo.dRLeadLeptonLeadJet.push_back(deltaR(multiLeptonMultiJet->leadingLepton()->eta(), multiLeptonMultiJet->leadingLepton()->phi(), multiLeptonMultiJet->leadingJet()->eta(), multiLeptonMultiJet->leadingJet()->phi()));
		evInfo.dRLeadLeptonSubLeadJet.push_back(deltaR(multiLeptonMultiJet->leadingLepton()->eta(), multiLeptonMultiJet->leadingLepton()->phi(), multiLeptonMultiJet->subLeadingJet()->eta(), multiLeptonMultiJet->subLeadingJet()->phi())); 
		evInfo.dRSubLeadLeptonLeadJet.push_back(deltaR(multiLeptonMultiJet->subLeadingLepton()->eta(), multiLeptonMultiJet->subLeadingLepton()->phi(), multiLeptonMultiJet->leadingJet()->eta(), multiLeptonMultiJet->leadingJet()->phi()));  
		evInfo.dRSubLeadLeptonSubLeadJet.push_back(deltaR(multiLeptonMultiJet->subLeadingLepton()->eta(), multiLeptonMultiJet->subLeadingLepton()->phi(), multiLeptonMultiJet->subLeadingJet()->eta(), multiLeptonMultiJet->subLeadingJet()->phi()));  

		evInfo.multiLeptonMultiJet_sumPt.push_back(sumPtCorr);  
		evInfo.multiLeptonMultiJet_invMass.push_back(massCorr);  
		evInfo.diLepton_invMass.push_back(diLeptonInvMass);  
		evInfo.diLepton_pt.push_back(diLeptonPt);  
		evInfo.diJet_invMass.push_back(diJetInvMass);
		evInfo.diJetLeadingLepton_invMass.push_back(diJetLeadingLeptonInvMass);  
		evInfo.diJetSubLeadingLepton_invMass.push_back(diJetSubLeadingLeptonInvMass);  

		evInfo.leadingEle_passHEEPId.push_back(leadingElePassHEEPId);
		evInfo.leadingEle_HEEPBitMapValues.push_back(leadingEleHeepBitMap);
		evInfo.leadingEle_passTightId.push_back(leadingElePassTightId);
		evInfo.leadingEle_passMediumId.push_back(leadingElePassMediumId);
		evInfo.leadingEle_passLooseId.push_back(leadingElePassLooseId);
		evInfo.leadingEle_passVetoId.push_back(leadingElePassVetoId);
		evInfo.leadingEle_passMVATightId.push_back(leadingElePassMVATightId);
		evInfo.leadingEle_passMVAMediumId.push_back(leadingElePassMVAMediumId);
		evInfo.leadingEle_idmva.push_back(leadingEleIdmva);
		evInfo.leadingEle_iso.push_back(leadingEleIso);
		evInfo.leadingEle_isMatchedToGen.push_back(leadingEleIsMatchedToGen);
		evInfo.leadingEle_etaSC.push_back(leadingEleEtaSC);
		evInfo.leadingEle_isEcalDriven.push_back(leadingEleIsEcalDriven);
		evInfo.leadingEle_dEtaIn.push_back(leadingEleDEtaIn);
		evInfo.leadingEle_dPhiIn.push_back(leadingEleDPhiIn);
		evInfo.leadingEle_hOverE.push_back(leadingEleHoE);
		evInfo.leadingEle_full5x5_r9.push_back(leadingEleR9); 
		evInfo.leadingEle_full5x5_sigmaIetaIeta.push_back(leadingEleSigmaIetaIeta);
		evInfo.leadingEle_full5x5_E5x5.push_back(leadingEleE5x5);
		evInfo.leadingEle_full5x5_E1x5.push_back(leadingEleE1x5);
		evInfo.leadingEle_full5x5_E2x5.push_back(leadingEleE2x5);
		evInfo.leadingEle_EmHadDepth1Iso.push_back(leadingEleEmHadDepth1Iso);
		evInfo.leadingEle_ptTracksIso.push_back(leadingElePtTracksIso);
		evInfo.leadingEle_innerLayerLostHits.push_back(leadingEleMissingHits);
		evInfo.leadingEle_dxy.push_back(leadingEleDxy);
		evInfo.leadingEle_dz.push_back(leadingEleDz);
		evInfo.leadingEle_eOverP.push_back(leadingEleEoverP);
		evInfo.leadingEle_ecalEnergy.push_back(leadingEleEcalEnergy);
		evInfo.leadingEle_hcalOverEcal.push_back(leadingEleHcalOverEcal);

		evInfo.subLeadingEle_passHEEPId.push_back(subLeadingElePassHEEPId);
		evInfo.subLeadingEle_HEEPBitMapValues.push_back(subLeadingEleHeepBitMap);		
		evInfo.subLeadingEle_passTightId.push_back(subLeadingElePassTightId);
		evInfo.subLeadingEle_passMediumId.push_back(subLeadingElePassMediumId);
		evInfo.subLeadingEle_passLooseId.push_back(subLeadingElePassLooseId);
		evInfo.subLeadingEle_passVetoId.push_back(subLeadingElePassVetoId);
		evInfo.subLeadingEle_passMVATightId.push_back(subLeadingElePassMVATightId);
		evInfo.subLeadingEle_passMVAMediumId.push_back(subLeadingElePassMVAMediumId);
		evInfo.subLeadingEle_idmva.push_back(subLeadingEleIdmva);
		evInfo.subLeadingEle_iso.push_back(subLeadingEleIso);
		evInfo.subLeadingEle_isMatchedToGen.push_back(subLeadingEleIsMatchedToGen);
		evInfo.subLeadingEle_etaSC.push_back(subLeadingEleEtaSC);		
		evInfo.subLeadingEle_isEcalDriven.push_back(subLeadingEleIsEcalDriven);
		evInfo.subLeadingEle_dEtaIn.push_back(subLeadingEleDEtaIn);
		evInfo.subLeadingEle_dPhiIn.push_back(subLeadingEleDPhiIn);
		evInfo.subLeadingEle_hOverE.push_back(subLeadingEleHoE);
		evInfo.subLeadingEle_full5x5_r9.push_back(subLeadingEleR9); 
		evInfo.subLeadingEle_full5x5_sigmaIetaIeta.push_back(subLeadingEleSigmaIetaIeta);
		evInfo.subLeadingEle_full5x5_E5x5.push_back(subLeadingEleE5x5);
		evInfo.subLeadingEle_full5x5_E1x5.push_back(subLeadingEleE1x5);
		evInfo.subLeadingEle_full5x5_E2x5.push_back(subLeadingEleE2x5);
		evInfo.subLeadingEle_EmHadDepth1Iso.push_back(subLeadingEleEmHadDepth1Iso);
		evInfo.subLeadingEle_ptTracksIso.push_back(subLeadingElePtTracksIso);
		evInfo.subLeadingEle_innerLayerLostHits.push_back(subLeadingEleMissingHits);
		evInfo.subLeadingEle_dxy.push_back(subLeadingEleDxy);
		evInfo.subLeadingEle_dz.push_back(subLeadingEleDz);
		evInfo.subLeadingEle_eOverP.push_back(subLeadingEleEoverP);
		evInfo.subLeadingEle_ecalEnergy.push_back(subLeadingEleEcalEnergy);
		evInfo.subLeadingEle_hcalOverEcal.push_back(subLeadingEleHcalOverEcal);

		evInfo.leadingMuon_iso.push_back(leadingMuonIso);
		evInfo.leadingMuon_PFiso.push_back(leadingMuonPFiso);
		evInfo.leadingMuon_isHighPt.push_back(leadingMuonIsHighPt);
		evInfo.leadingMuon_isTight.push_back(leadingMuonIsTight);
		evInfo.leadingMuon_isMedium.push_back(leadingMuonIsMedium);
		evInfo.leadingMuon_isLoose.push_back(leadingMuonIsLoose);
		evInfo.leadingMuon_isMatchedToGen.push_back(leadingMuonIsMatchedToGen);
		evInfo.leadingMuon_dz.push_back(leadingMuonDz);
		evInfo.leadingMuon_dxy.push_back(leadingMuonDxy);
		evInfo.leadingMuon_RochCor.push_back(leadingMuonRochCor);

		evInfo.subLeadingMuon_iso.push_back(subLeadingMuonIso);
		evInfo.subLeadingMuon_PFiso.push_back(subLeadingMuonPFiso);
		evInfo.subLeadingMuon_isHighPt.push_back(subLeadingMuonIsHighPt);
		evInfo.subLeadingMuon_isTight.push_back(subLeadingMuonIsTight);
		evInfo.subLeadingMuon_isMedium.push_back(subLeadingMuonIsMedium);
		evInfo.subLeadingMuon_isLoose.push_back(subLeadingMuonIsLoose);
		evInfo.subLeadingMuon_isMatchedToGen.push_back(subLeadingMuonIsMatchedToGen);
		evInfo.subLeadingMuon_dz.push_back(subLeadingMuonDz);
		evInfo.subLeadingMuon_dxy.push_back(subLeadingMuonDxy);
		evInfo.subLeadingMuon_RochCor.push_back(subLeadingMuonRochCor);

		// cout << "arriva a linea " << __LINE__ << endl;

	}
	// cout << "exit MLMJ loop" << endl;


	// --- fill the tree  
	if ( passTrigger ) eventTree->Fill(); 

	// cout << "fillo tree" << endl;

}
// ******************************************************************************************



// ******************************************************************************************
void miniTreeMaker_multiLeptonMultiJet::endJob() {
	cout << "Total number of generated MLMJ before preselection = "<< ngen << endl;
	cout << "Number of generated MLMJ after preselection        = "<< ngenPre << endl;
	cout << "Total number of MLMJ before preselection           = "<< ndldj << endl;
	cout << "Number of MLMJ after preselection                  = "<< npre << endl; 
} // end of endJob
// ******************************************************************************************



// ******************************************************************************************
void miniTreeMaker_multiLeptonMultiJet::initEventStructure() {
	// per-event tree:
	evInfo.run = -999;
	evInfo.event = -999.;
	evInfo.lumi = -999.;
	evInfo.rho = -999.;
	evInfo.weight = -999.;
	evInfo.puweight = -999.;
	evInfo.nvtx = -999;
	evInfo.npu = -999;
	evInfo.passEEJJhlt = -1;
	evInfo.passMMJJhlt = -1;
	evInfo.passEMJJhlt = -1;
	evInfo.passTandPEEhlt = -1;
	evInfo.passTandPMMhlt = -1;

	evInfo.vtx_x .clear();
	evInfo.vtx_y .clear();
	evInfo.vtx_z .clear();
	
	evInfo.isEEJJ .clear();
	evInfo.isEETT .clear();
	evInfo.isMMJJ .clear();
	evInfo.isMMTT .clear();
	evInfo.isEMJJ .clear();

	evInfo.isSignalRegion .clear();
	evInfo.isLowMllCR .clear();
	evInfo.isLowMlljjCR .clear();

	evInfo.isBB .clear();
	evInfo.isEE .clear();
	evInfo.isEB .clear();

	evInfo.passPreselections .clear();

	evInfo.leadingLepton_e .clear(); 
	evInfo.leadingLepton_pt .clear();
	evInfo.leadingLepton_eta .clear(); 
	evInfo.leadingLepton_phi .clear();  
	evInfo.leadingLepton_charge .clear();

	evInfo.subLeadingLepton_e .clear();
	evInfo.subLeadingLepton_pt .clear();
	evInfo.subLeadingLepton_eta .clear(); 
	evInfo.subLeadingLepton_phi .clear();  
	evInfo.subLeadingLepton_charge .clear();

	evInfo.leadingJet_e .clear(); 
	evInfo.leadingJet_pt .clear();
	evInfo.leadingJet_eta .clear(); 
	evInfo.leadingJet_phi .clear();
	evInfo.leadingJet_isMatchedToGen .clear();
	evInfo.leadingJet_isThight .clear(); 

	evInfo.subLeadingJet_e .clear();
	evInfo.subLeadingJet_pt .clear();
	evInfo.subLeadingJet_eta .clear();
	evInfo.subLeadingJet_phi .clear();  
	evInfo.subLeadingJet_isMatchedToGen .clear();
	evInfo.subLeadingJet_isThight .clear(); 

	evInfo.dRLeadLeptonLeadJet .clear();
	evInfo.dRLeadLeptonSubLeadJet .clear();
	evInfo.dRSubLeadLeptonLeadJet .clear();
	evInfo.dRSubLeadLeptonSubLeadJet .clear();

	evInfo.multiLeptonMultiJet_sumPt .clear();
	evInfo.multiLeptonMultiJet_invMass .clear();
	evInfo.diLepton_invMass .clear();
	evInfo.diLepton_pt .clear();
	evInfo.diJet_invMass .clear();
	evInfo.diJetLeadingLepton_invMass .clear();
	evInfo.diJetSubLeadingLepton_invMass .clear();

	evInfo.leadingEle_passHEEPId .clear();
	evInfo.leadingEle_HEEPBitMapValues .clear();
	evInfo.leadingEle_passTightId .clear();
	evInfo.leadingEle_passMediumId .clear();
	evInfo.leadingEle_passLooseId .clear();
	evInfo.leadingEle_passVetoId .clear();
	evInfo.leadingEle_passMVATightId .clear();
	evInfo.leadingEle_passMVAMediumId .clear();
	evInfo.leadingEle_idmva .clear();
	evInfo.leadingEle_iso .clear();
	evInfo.leadingEle_isMatchedToGen .clear();
	evInfo.leadingEle_etaSC .clear();
	evInfo.leadingEle_isEcalDriven .clear();
	evInfo.leadingEle_dEtaIn .clear();  
	evInfo.leadingEle_dPhiIn .clear();
	evInfo.leadingEle_hOverE .clear();
	evInfo.leadingEle_full5x5_r9 .clear();  
	evInfo.leadingEle_full5x5_sigmaIetaIeta .clear();
	evInfo.leadingEle_full5x5_E5x5 .clear();
	evInfo.leadingEle_full5x5_E1x5 .clear();
	evInfo.leadingEle_full5x5_E2x5 .clear();
	evInfo.leadingEle_EmHadDepth1Iso .clear();
	evInfo.leadingEle_ptTracksIso .clear();
	evInfo.leadingEle_innerLayerLostHits .clear();
	evInfo.leadingEle_dxy .clear();
	evInfo.leadingEle_dz .clear();
	evInfo.leadingEle_eOverP .clear();
	evInfo.leadingEle_ecalEnergy .clear();
	evInfo.leadingEle_hcalOverEcal .clear();

	evInfo.subLeadingEle_passHEEPId .clear();
	evInfo.subLeadingEle_HEEPBitMapValues .clear();
	evInfo.subLeadingEle_passTightId .clear();
	evInfo.subLeadingEle_passMediumId .clear();
	evInfo.subLeadingEle_passLooseId .clear();
	evInfo.subLeadingEle_passVetoId .clear();
	evInfo.subLeadingEle_passMVATightId .clear();
	evInfo.subLeadingEle_passMVAMediumId .clear();
	evInfo.subLeadingEle_idmva . clear();
	evInfo.subLeadingEle_iso .clear();
	evInfo.subLeadingEle_isMatchedToGen .clear();
	evInfo.subLeadingEle_etaSC .clear();
	evInfo.subLeadingEle_isEcalDriven .clear();
	evInfo.subLeadingEle_dEtaIn .clear();  
	evInfo.subLeadingEle_dPhiIn .clear();
	evInfo.subLeadingEle_hOverE .clear();
	evInfo.subLeadingEle_full5x5_r9 .clear();  
	evInfo.subLeadingEle_full5x5_sigmaIetaIeta .clear();
	evInfo.subLeadingEle_full5x5_E5x5 .clear();
	evInfo.subLeadingEle_full5x5_E1x5 .clear();
	evInfo.subLeadingEle_full5x5_E2x5 .clear();
	evInfo.subLeadingEle_EmHadDepth1Iso .clear();
	evInfo.subLeadingEle_ptTracksIso .clear();
	evInfo.subLeadingEle_innerLayerLostHits .clear();
	evInfo.subLeadingEle_dxy .clear();
	evInfo.subLeadingEle_dz .clear();
	evInfo.subLeadingEle_eOverP .clear();
	evInfo.subLeadingEle_ecalEnergy .clear();	
	evInfo.subLeadingEle_hcalOverEcal .clear();

	evInfo.leadingMuon_iso .clear();
	evInfo.leadingMuon_PFiso .clear();
	evInfo.leadingMuon_isHighPt .clear();
	evInfo.leadingMuon_isTight .clear();
	evInfo.leadingMuon_isMedium .clear();
	evInfo.leadingMuon_isLoose .clear();
	evInfo.leadingMuon_isMatchedToGen .clear();
	evInfo.leadingMuon_dz .clear();
	evInfo.leadingMuon_dxy .clear();
	evInfo.leadingMuon_RochCor .clear();

	evInfo.subLeadingMuon_iso .clear();
	evInfo.subLeadingMuon_PFiso .clear();
	evInfo.subLeadingMuon_isHighPt .clear();
	evInfo.subLeadingMuon_isTight .clear();
	evInfo.subLeadingMuon_isMedium .clear();
	evInfo.subLeadingMuon_isLoose .clear();
	evInfo.subLeadingMuon_isMatchedToGen .clear();
	evInfo.subLeadingMuon_dz .clear();
	evInfo.subLeadingMuon_dxy .clear();
	evInfo.subLeadingMuon_RochCor .clear();

}
// ******************************************************************************************
