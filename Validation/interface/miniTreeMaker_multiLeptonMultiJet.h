#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "flashgg/Taggers/interface/GlobalVariablesDumper.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "dafne/Systematics/interface/RochCor.h"


using namespace std;
using namespace edm;
using namespace flashgg;


// define the structures used to create tree branches and fill the trees
struct eventInfo {

	int run;
	int event;
	int lumi;
	float rho;
	float weight;
	float puweight;
	int nvtx;
	int npu;

	vector<float> vtx_x;
	vector<float> vtx_y;
	vector<float> vtx_z;

	vector<float> ele_e;
	vector<float> ele_pt;
	vector<float> ele_eta;
	vector<float> ele_phi;
	vector<bool> ele_passHEEPId; 
	vector<bool> ele_passMediumId;
	vector<float> ele_iso;
	vector<float> ele_dz;
	vector<float> ele_d0;
	vector<int> ele_isMatchedToGen;
	vector<int> ele_charge;
	vector<float> ele_etaSC;
	vector<bool> ele_isEcalDriven;
	vector<float> ele_dEtaIn; 
	vector<float> ele_dPhiIn;
	vector<float> ele_hOverE;
	vector<float> ele_full5x5_r9;	
	vector<float> ele_full5x5_sigmaIetaIeta;
	vector<float> ele_full5x5_E5x5;
	vector<float> ele_full5x5_E1x5;
	vector<float> ele_full5x5_E2x5;  
	vector<float> ele_EmHadDepth1Iso;
	vector<float> ele_ptTracksIso;
	vector<int> ele_innerLayerLostHits;
	vector<float> ele_dxy;
	vector<float> ele_eOverP;
	vector<float> ele_ecalEnergy;
	vector<float> ele_hcalOverEcal;

	vector<float> mu_e;
	vector<float> mu_pt;
	vector<float> mu_eta;
	vector<float> mu_phi;
	vector<float> mu_iso;
	vector<float> mu_PFiso;
	vector<bool> mu_isHighPt;
	vector<int> mu_isMatchedToGen;
	vector<int> mu_charge;
	vector<float> mu_dz;
	vector<float> mu_dxy;
	vector<float> mu_RochCor;

	vector<float> jet_e;
	vector<float> jet_pt;
	vector<float> jet_eta;
	vector<float> jet_phi;
	vector<int>   jet_isMatchedToGen;
	vector<bool>  jet_isThight;

	vector<bool> isEEJJ;
	vector<bool> isEETT;
	vector<bool> isMMJJ;
	vector<bool> isMMTT;
	vector<bool> isEMJJ;

	vector<bool> isSignalRegion;
	vector<bool> isLowMllCR;
	vector<bool> isLowMlljjCR;

	vector<bool> isBB;
	vector<bool> isEE;
	vector<bool> isEB;

	vector<bool> passPreselections;

	vector<float> leadingLepton_e;
	vector<float> leadingLepton_pt;
	vector<float> leadingLepton_eta;  
	vector<float> leadingLepton_phi; 
	vector<float> leadingLepton_charge;

	vector<float> subLeadingLepton_e;
	vector<float> subLeadingLepton_pt;
	vector<float> subLeadingLepton_eta;  
	vector<float> subLeadingLepton_phi;  
	vector<float> subLeadingLepton_charge;

	vector<float> leadingJet_e;
	vector<float> leadingJet_pt;
	vector<float> leadingJet_eta;  
	vector<float> leadingJet_phi;  
	vector<int>   leadingJet_isMatchedToGen;
	vector<bool>  leadingJet_isThight;

	vector<float> subLeadingJet_e;
	vector<float> subLeadingJet_pt;
	vector<float> subLeadingJet_eta;  
	vector<float> subLeadingJet_phi;  
	vector<int>   subLeadingJet_isMatchedToGen;
	vector<bool>  subLeadingJet_isThight;

	vector<float> dRLeadLeptonLeadJet;
	vector<float> dRLeadLeptonSubLeadJet;
	vector<float> dRSubLeadLeptonLeadJet;
	vector<float> dRSubLeadLeptonSubLeadJet;

	vector<float> multiLeptonMultiJet_sumPt;
	vector<float> multiLeptonMultiJet_invMass;
	vector<float> diLepton_invMass;
	vector<float> diLepton_pt;
	vector<float> diJet_invMass;
	vector<float> diJetLeadingLepton_invMass;
	vector<float> diJetSubLeadingLepton_invMass;

	vector<bool> leadingEle_passHEEPId;
	vector<bool> leadingEle_passMediumId;
	vector<float> leadingEle_iso;
	vector<int> leadingEle_isMatchedToGen;
	vector<float> leadingEle_etaSC;
	vector<bool> leadingEle_isEcalDriven;
	vector<float> leadingEle_dEtaIn;  
	vector<float> leadingEle_dPhiIn;
	vector<float> leadingEle_hOverE;
	vector<float> leadingEle_full5x5_r9;  
	vector<float> leadingEle_full5x5_sigmaIetaIeta;
	vector<float> leadingEle_full5x5_E5x5;
	vector<float> leadingEle_full5x5_E1x5;
	vector<float> leadingEle_full5x5_E2x5;
	vector<float> leadingEle_EmHadDepth1Iso;
	vector<float> leadingEle_ptTracksIso;
	vector<int> leadingEle_innerLayerLostHits;
	vector<float> leadingEle_dxy;
	vector<float> leadingEle_dz;
	vector<float> leadingEle_eOverP;
	vector<float> leadingEle_ecalEnergy;
	vector<float> leadingEle_hcalOverEcal;

	vector<bool> subLeadingEle_passHEEPId;	
	vector<bool> subLeadingEle_passMediumId;
	vector<float> subLeadingEle_iso;
	vector<int> subLeadingEle_isMatchedToGen;
	vector<float> subLeadingEle_etaSC;			
	vector<bool> subLeadingEle_isEcalDriven;
	vector<float> subLeadingEle_dEtaIn;  
	vector<float> subLeadingEle_dPhiIn;
	vector<float> subLeadingEle_hOverE;
	vector<float> subLeadingEle_full5x5_r9; 
	vector<float> subLeadingEle_full5x5_sigmaIetaIeta;
	vector<float> subLeadingEle_full5x5_E5x5;
	vector<float> subLeadingEle_full5x5_E1x5;
	vector<float> subLeadingEle_full5x5_E2x5;
	vector<float> subLeadingEle_EmHadDepth1Iso;	
	vector<float> subLeadingEle_ptTracksIso;	
	vector<int> subLeadingEle_innerLayerLostHits;
	vector<float> subLeadingEle_dxy;	
	vector<float> subLeadingEle_dz;
	vector<float> subLeadingEle_eOverP;
	vector<float> subLeadingEle_ecalEnergy;
	vector<float> subLeadingEle_hcalOverEcal;

	vector<float> leadingMuon_iso;
	vector<float> leadingMuon_PFiso;
	vector<bool> leadingMuon_isHighPt;
	vector<int> leadingMuon_isMatchedToGen;
	vector<float> leadingMuon_dz;
	vector<float> leadingMuon_dxy;
	vector<float> leadingMuon_RochCor;

	vector<float> subLeadingMuon_iso;
	vector<float> subLeadingMuon_PFiso;
	vector<bool> subLeadingMuon_isHighPt;
	vector<int> subLeadingMuon_isMatchedToGen;
	vector<float> subLeadingMuon_dz;
	vector<float> subLeadingMuon_dxy;
	vector<float> subLeadingMuon_RochCor;

};


// ************************** 
int electronMatchingToGen(Ptr<flashgg::Electron> electron,  Handle<View<reco::GenParticle> > genParticles){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
		if ( fabs(gen->pdgId()) != 11 ) continue;
		if ( !(gen->isPromptFinalState())) continue;
		float dR = deltaR( electron->eta(), electron->phi(), gen->eta(), gen->phi() );
		if (dR < 0.1){ //??? 0.1 ok???
			mcmatch = 1;
		}
	}
	return (mcmatch);
}
// ******************************************************************************************

// ************************** 
int electronMatchingToGen(const flashgg::Electron* electron,  Handle<View<reco::GenParticle> > genParticles){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
		if ( fabs(gen->pdgId()) != 11 ) continue;
		if ( !(gen->isPromptFinalState())) continue;
		float dR = deltaR( electron->eta(), electron->phi(), gen->eta(), gen->phi() );
		if (dR < 0.1){ 
			mcmatch = 1;
		}
	}
	return (mcmatch);
}
// ******************************************************************************************



// *********************** 
int muonMatchingToGen(Ptr<flashgg::Muon> muon, Handle<View<reco::GenParticle> > genParticles){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
			if ( fabs(gen->pdgId()) != 13 ) continue;
			if ( !(gen)->isPromptFinalState()) continue;
			float dR = deltaR( muon->eta(), muon->phi(), gen->eta(), gen->phi() );
			//cout << " *** Found muon:  ***"<<endl;
			//cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
			//cout << "dR = " << dR <<endl;
			if (dR < 0.1){ //??? 0.1 ok???
				mcmatch = 1;
			}
	}
	return (mcmatch);
}
// ******************************************************************************************

// ************************** 
int muonMatchingToGen(const flashgg::Muon* muon, Handle<View<reco::GenParticle> > genParticles){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
			if ( fabs(gen->pdgId()) != 13 ) continue;
			if ( !(gen)->isPromptFinalState()) continue;
			float dR = deltaR( muon->eta(), muon->phi(), gen->eta(), gen->phi() );
			//cout << " *** Found muon:  ***"<<endl;
			//cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
			//cout << "dR = " << dR <<endl;
			if (dR < 0.1){ 
				mcmatch = 1;
			}
	}
	return (mcmatch);
}
// ******************************************************************************************



// // ************************** 
int jetMatchingToGen(Ptr<flashgg::Jet> jet,  Handle<View<reco::GenJet> > genJets){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genJets->size(); i++ ) {
		Ptr<reco::GenJet> genJet = genJets->ptrAt(i);
		float dR = deltaR( jet->eta(), jet->phi(), genJet->eta(), genJet->phi() );
		if (dR > 0.4) continue;
		mcmatch = 1;
	}
	return (mcmatch);
}
// ******************************************************************************************

// ************************** 
int jetMatchingToGen(const flashgg::Jet* jet,  Handle<View<reco::GenJet> > genJets){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genJets->size(); i++ ) {
		Ptr<reco::GenJet> genJet = genJets->ptrAt(i);
		float dR = deltaR( jet->eta(), jet->phi(), genJet->eta(), genJet->phi() );
		if (dR > 0.4) continue;
		mcmatch = 1;
	}
	return (mcmatch);
}
// ******************************************************************************************



// **************** 
Ptr<reco::Vertex> chooseBestVtx(const vector<Ptr<reco::Vertex> > &vertices, Ptr<flashgg::Electron> electron){
	double vtx_dz = 1000000.;
	unsigned int min_dz_vtx = -1;			
	for( unsigned int vtxi = 0; vtxi < vertices.size(); vtxi++ ) {            
		Ptr<reco::Vertex> vtx = vertices[vtxi];            
		if( vtx_dz > fabs(electron->gsfTrack()->dz( vtx->position() )) ) {                
			vtx_dz = fabs(electron->gsfTrack()->dz( vtx->position() ) );
			min_dz_vtx = vtxi;
		}
	}					
	return vertices[min_dz_vtx];
}
// *****************************************************************************************

// ************************** 
Ptr<reco::Vertex> chooseBestVtx(const vector<Ptr<reco::Vertex> > &vertices, const flashgg::Electron* electron){
	double vtx_dz = 1000000.;
	unsigned int min_dz_vtx = -1;			
	for( unsigned int vtxi = 0; vtxi < vertices.size(); vtxi++ ) {            
		Ptr<reco::Vertex> vtx = vertices[vtxi];            
		if( vtx_dz > fabs(electron->gsfTrack()->dz( vtx->position() )) ) {                
			vtx_dz = fabs(electron->gsfTrack()->dz( vtx->position() ) );
			min_dz_vtx = vtxi;
		}
	}					
	return vertices[min_dz_vtx];
}
// ******************************************************************************************



// **************** 
Ptr<reco::Vertex> chooseBestMuonVtx(const vector<Ptr<reco::Vertex> > &vertices, Ptr<flashgg::Muon> muon){
	int vtxInd = 0;
	double dzmin = 9999;
	for( size_t ivtx = 0 ; ivtx < vertices.size(); ivtx++ ) {
		Ptr<reco::Vertex> vtx = vertices[ivtx];
		if( !muon->innerTrack() ) { continue; }
		if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {
			dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
			vtxInd = ivtx;
		}
	}
	return vertices[vtxInd];
}
// ******************************************************************************************

// ************************** 
Ptr<reco::Vertex> chooseBestMuonVtx(const vector<Ptr<reco::Vertex> > &vertices, const flashgg::Muon* muon){
	int vtxInd = 0;
	double dzmin = 9999;
	for( size_t ivtx = 0 ; ivtx < vertices.size(); ivtx++ ) {
		Ptr<reco::Vertex> vtx = vertices[ivtx];
		if( !muon->innerTrack() ) { continue; }
		if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {
			dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
			vtxInd = ivtx;
		}
	}
	return vertices[vtxInd];
}
// ******************************************************************************************



// **************** 
float RochesterCorrection(Ptr<flashgg::Muon> muon, Handle<View<reco::GenParticle> > genParticles, bool isData){

	RoccoR rc("/afs/cern.ch/user/g/gnegro/work/NuAnalysis-Moriond17/CMSSW_8_0_26_patch1/src/dafne/data/rcdata.2016.v3");
	// cout << "Muon Pt Before = " << muon->pt() << ", Muon Eta Before = " << muon->eta() << endl;

	TRandom *gRandom = new TRandom();
	gRandom->SetSeed(1); 	
	float u1 = gRandom->Rndm();
	float u2 = gRandom->Rndm();
	int nl = muon->bestTrack()->hitPattern().trackerLayersWithMeasurement();

	double RochCor = 1.;

	if( isData ) {
		RochCor = rc.kScaleDT(muon->charge(), muon->pt(), muon->eta(), muon->phi(), 0, 0);
	} else {

		float genMuPt = 0.;
		for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
	   		Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
		   if( fabs(gen->pdgId()) == 13 && deltaR(gen->eta(), gen->phi(), muon->eta(), muon->phi()) < 0.2 ) genMuPt = gen->pt();
		}    

		if (genMuPt != 0) RochCor = rc.kScaleFromGenMC(muon->charge(), muon->pt(), muon->eta(), muon->phi(), nl, genMuPt, u1, 0, 0);
		else RochCor = rc.kScaleAndSmearMC(muon->charge(), muon->pt(), muon->eta(), muon->phi(), nl, u1, u2, 0, 0);
	}	
	// cout << "Muon Pt After = " << muon->pt()*RochCor << ", Muon Eta After = " << muon->eta() << endl;

	return RochCor;
}
// ******************************************************************************************

// ************************** 
float RochesterCorrection(const flashgg::Muon* muon, Handle<View<reco::GenParticle> > genParticles, bool isData){

	RoccoR rc("/afs/cern.ch/user/g/gnegro/work/NuAnalysis-Moriond17/CMSSW_8_0_26_patch1/src/dafne/data/rcdata.2016.v3");
	// cout << "Muon Pt Before = " << muon->pt() << ", Muon Eta Before = " << muon->eta() << endl;

	TRandom *gRandom = new TRandom();
	gRandom->SetSeed(1); 	
	float u1 = gRandom->Rndm();
	float u2 = gRandom->Rndm();
	int nl = muon->bestTrack()->hitPattern().trackerLayersWithMeasurement();

	double RochCor = 1.;

	if( isData ) {
		RochCor = rc.kScaleDT(muon->charge(), muon->pt(), muon->eta(), muon->phi(), 0, 0);
	} else {

		float genMuPt = 0.;
		for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		   Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
		   if( fabs(gen->pdgId()) == 13 && deltaR(gen->eta(), gen->phi(), muon->eta(), muon->phi()) < 0.2 ) genMuPt = gen->pt();
		}    

		if (genMuPt != 0) RochCor = rc.kScaleFromGenMC(muon->charge(), muon->pt(), muon->eta(), muon->phi(), nl, genMuPt, u1, 0, 0);
		else RochCor = rc.kScaleAndSmearMC(muon->charge(), muon->pt(), muon->eta(), muon->phi(), nl, u1, u2, 0, 0);
	}	
	// cout << "Muon Pt After = " << muon->pt()*RochCor << ", Muon Eta After = " << muon->eta() << endl;

	return RochCor;
}
// ******************************************************************************************



// ***************************** 
float electronIsolation(Ptr<flashgg::Electron> electron, double rho){
	// -- compute combined relative isolation: IsoCh + max( 0.0, IsoNh + IsoPh - PU ) )/pT, PU = rho * Aeff 
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	// effective areas:  https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt

	float Aeff = 0;
	float eta = fabs(electron->eta());
	if( eta <  1.0 )                  { Aeff = 0.1703; }
	if( eta >= 1.0   && eta < 1.479 ) { Aeff = 0.1715; }
	if( eta >= 1.479 && eta < 2.0 )   { Aeff = 0.1213; }
	if( eta >= 2.0   && eta < 2.2 )   { Aeff = 0.1230; }
	if( eta >= 2.2   && eta < 2.3 )   { Aeff = 0.1635; }
	if( eta >= 2.3   && eta < 2.4 )   { Aeff = 0.1937; }
	if( eta >= 2.4 )                  { Aeff = 0.2393; }

	//float iso = electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ); 
	reco::GsfElectron::PflowIsolationVariables pfIso = electron->pfIsolationVariables();
	float iso = pfIso.sumChargedHadronPt + max( pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * Aeff, 0. );

	//cout << electron->chargedHadronIso() << "  " <<  pfIso.sumChargedHadronPt << "   pt = " << electron->pt() << endl; 
	//cout << electron->neutralHadronIso() << "  " << pfIso.sumNeutralHadronEt << endl;
	//cout << electron->photonIso() << "  " << pfIso.sumPhotonEt <<endl;
	//cout << electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ) << "   "<< iso<< endl;

	return (iso/ electron->pt());	
}
// ******************************************************************************************

// ************************** 
float electronIsolation(const flashgg::Electron* electron, double rho){
	// -- compute combined relative isolation: IsoCh + max( 0.0, IsoNh + IsoPh - PU ) )/pT, PU = rho * Aeff 
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	// effective areas:  https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt

	float Aeff = 0;
	float eta = fabs(electron->eta());
	if( eta <  1.0 )                  { Aeff = 0.1703; }
	if( eta >= 1.0   && eta < 1.479 ) { Aeff = 0.1715; }
	if( eta >= 1.479 && eta < 2.0 )   { Aeff = 0.1213; }
	if( eta >= 2.0   && eta < 2.2 )   { Aeff = 0.1230; }
	if( eta >= 2.2   && eta < 2.3 )   { Aeff = 0.1635; }
	if( eta >= 2.3   && eta < 2.4 )   { Aeff = 0.1937; }
	if( eta >= 2.4 )                  { Aeff = 0.2393; }

	//float iso = electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ); 
	reco::GsfElectron::PflowIsolationVariables pfIso = electron->pfIsolationVariables();
	float iso = pfIso.sumChargedHadronPt + max( pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * Aeff, 0. );

	//cout << electron->chargedHadronIso() << "  " <<  pfIso.sumChargedHadronPt << "   pt = " << electron->pt() << endl; 
	//cout << electron->neutralHadronIso() << "  " << pfIso.sumNeutralHadronEt << endl;
	//cout << electron->photonIso() << "  " << pfIso.sumPhotonEt <<endl;
	//cout << electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ) << "   "<< iso<< endl;

	return (iso/ electron->pt());	
}
// ******************************************************************************************


// ************************* 
bool passMultiLeptonMultiJetPreselection(Ptr<flashgg::MultiLeptonMultiJetCandidate> mlmj, math::XYZTLorentzVector leadingLepton, math::XYZTLorentzVector subLeadingLepton){
	if (leadingLepton.pt() < 60.) return false;
	if (subLeadingLepton.pt() < 53.) return false;
	if (mlmj->leadingJet()->pt() < 40.) return false;
	if (mlmj->subLeadingJet()->pt() < 40.) return false;

	if (fabs(mlmj->leadingLepton()->eta()) > 2.4) return false;
	if (fabs(mlmj->subLeadingLepton()->eta()) > 2.4) return false;
	if (fabs(mlmj->leadingJet()->eta()) > 2.4) return false;
	if (fabs(mlmj->subLeadingJet()->eta()) > 2.4) return false;


	float dRLeadLeptonLeadJet = deltaR(mlmj->leadingLepton()->eta(), mlmj->leadingLepton()->phi(), mlmj->leadingJet()->eta(), mlmj->leadingJet()->phi());
	float dRLeadLeptonSubLeadJet = deltaR(mlmj->leadingLepton()->eta(), mlmj->leadingLepton()->phi(), mlmj->subLeadingJet()->eta(), mlmj->subLeadingJet()->phi());
	float dRSubLeadLeptonLeadJet = deltaR(mlmj->subLeadingLepton()->eta(), mlmj->subLeadingLepton()->phi(), mlmj->leadingJet()->eta(), mlmj->leadingJet()->phi());
	float dRSubLeadLeptonSubLeadJet = deltaR(mlmj->subLeadingLepton()->eta(), mlmj->subLeadingLepton()->phi(), mlmj->subLeadingJet()->eta(), mlmj->subLeadingJet()->phi());

	if( dRLeadLeptonLeadJet < 0.4 )  return false;
	if( dRLeadLeptonSubLeadJet < 0.4 )  return false;
	if( dRSubLeadLeptonLeadJet < 0.4 )  return false;
	if( dRSubLeadLeptonSubLeadJet < 0.4 )  return false;

	return true;
}
// ******************************************************************************************


// ************************** 
bool isBB(Ptr<flashgg::MultiLeptonMultiJetCandidate> mlmj){
	if ( (fabs(mlmj->leadingLepton()->eta())<1.4442) 
		&& (fabs(mlmj->subLeadingLepton()->eta())<1.4442) 
	) return true;

	return false;
}
// ******************************************************************************************


// ************************** 
bool isEE(Ptr<flashgg::MultiLeptonMultiJetCandidate> mlmj){
	if ( (fabs(mlmj->leadingLepton()->eta())>1.566 && fabs(mlmj->leadingLepton()->eta())<2.5)  
		&& (fabs(mlmj->subLeadingLepton()->eta())>1.566 && fabs(mlmj->subLeadingLepton()->eta())<2.5)
	) return true;
	
	return false;
}
// ******************************************************************************************


// ************************** 
bool isEB(Ptr<flashgg::MultiLeptonMultiJetCandidate> mlmj){
	if ( ((fabs(mlmj->leadingLepton()->eta())<1.4442) && (fabs(mlmj->subLeadingLepton()->eta())>1.566 && fabs(mlmj->subLeadingLepton()->eta())<2.5)) 
		|| ((fabs(mlmj->leadingLepton()->eta())>1.566 && fabs(mlmj->leadingLepton()->eta())<2.5) && (fabs(mlmj->subLeadingLepton()->eta())<1.4442))
	) return true;

	return false;
}
// ******************************************************************************************


// ************************** 
bool isSignalRegion(float mlmjMass, float diLeptonInvMass){
	if (mlmjMass > 600. && diLeptonInvMass > 200.) {
		return true;
	}
	return false;
}
// ******************************************************************************************


// ************************** 
bool isLowMllCR(float diLeptonInvMass){
	if (diLeptonInvMass < 200.) return true;
	return false;
}
// ******************************************************************************************


// ************************** 
bool isLowMlljjCR(float mlmjMass, float diLeptonInvMass){
	if (mlmjMass < 600. && diLeptonInvMass > 200.) return true;
	return false;
}
// ******************************************************************************************




class miniTreeMaker_multiLeptonMultiJet : public BasicAnalyzer 
{
 public:
	miniTreeMaker_multiLeptonMultiJet( const ParameterSet & iConfig, TFileDirectory& fs, ConsumesCollector && cc);
	virtual ~miniTreeMaker_multiLeptonMultiJet();
	void beginJob();
	void analyze( const EventBase& event );
	void endJob();
 
 private:
	void initEventStructure();

	TTree *eventTree;
	eventInfo evInfo;
	int ngen;
	int ngenPre;
	int ndldj;
	int npre;
	int nEvents, nEventsPassingTrigger;

	EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
	EDGetTokenT<GenEventInfoProduct> genInfoToken_;
	EDGetTokenT<View<PileupSummaryInfo> >  PileUpToken_;
	EDGetTokenT<View<reco::Vertex> > vertexToken_;
	EDGetTokenT<View<flashgg::MultiLeptonMultiJetCandidate> > MultiLeptonMultiJetToken_; 
	EDGetTokenT<View<flashgg::Jet> > jetsToken_;
	EDGetTokenT<View<reco::GenJet> > genJetToken_;
	EDGetTokenT<View<Electron> > electronToken_;
	EDGetTokenT<View<Muon> > muonToken_;
	EDGetTokenT<TriggerResults> triggerBitsToken_;
	EDGetTokenT<double> rhoToken_;
	double lumiWeight_;
	bool saveHEEPvariables_;
	bool isDoubleEGinSignalRegion_;
	GlobalVariablesDumper *globalVarsDumper_;
};
// ******************************************************************************************                                                                                                                     


