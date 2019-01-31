#ifndef _makePlots_multiLeptonMultiJet_h_
#define _makePlots_multiLeptonMultiJet_h_
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <vector>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include "dafne/Validation/interface/functions.h"
#include "dafne/Validation/interface/functions_multiLeptonMultiJet.h"



using namespace std;


class makePlots_multiLeptonMultiJet {
	public:
	// TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		TChain          *fChain; 
		Int_t           fCurrent; //!current Tree number in a TChain

	// Declaration of leaf types
		Int_t           run;
		Int_t           event;
		Int_t           lumi;
		Float_t         rho;
		Float_t         weight;
		Float_t         puweight;
		Int_t           nvtx;
		Int_t           npu;
		vector<float>   *vtx_x;
		vector<float>   *vtx_y;
		vector<float>   *vtx_z;  
		vector<float>   *ele_e;
		vector<float>   *ele_pt;
		vector<float>   *ele_eta;
		vector<float>   *ele_phi;
		vector<bool>    *ele_passHEEPId;
		vector<bool>    *ele_passMediumId;
		vector<float>   *ele_iso;
		vector<float>   *ele_dz;
		vector<float>   *ele_d0;
		vector<int>     *ele_isMatchedToGen;
		vector<int>     *ele_charge;
		vector<float>   *ele_etaSC;
		vector<bool>    *ele_isEcalDriven;
		vector<float>   *ele_dEtaIn;
		vector<float>   *ele_dPhiIn;
		vector<float>   *ele_hOverE;
		vector<float>   *ele_full5x5_r9;
		vector<float>   *ele_full5x5_sigmaIetaIeta;
		vector<float>   *ele_full5x5_E5x5;
		vector<float>   *ele_full5x5_E1x5;
		vector<float>   *ele_full5x5_E2x5;
		vector<float>   *ele_EmHadDepth1Iso;
		vector<float>   *ele_ptTracksIso;
		vector<int>     *ele_innerLayerLostHits;
		vector<float>   *ele_dxy;
		vector<float>   *ele_eOverP;
		vector<float>   *ele_ecalEnergy;
		vector<float>   *ele_hcalOverEcal;
		vector<float>   *mu_e;
		vector<float>   *mu_pt;
		vector<float>   *mu_eta;
		vector<float>   *mu_phi;
		vector<float>   *mu_iso;
		vector<float>   *mu_PFiso;
		vector<bool>    *mu_isHighPt;
		vector<int>     *mu_isMatchedToGen;
		vector<int>     *mu_charge;
		vector<float>   *mu_dz;
		vector<float>   *mu_dxy;
		vector<float>   *mu_RochCor;
		vector<float>   *jet_e;
		vector<float>   *jet_pt;
		vector<float>   *jet_eta;
		vector<float>   *jet_phi;
		vector<float>   *jet_bdiscriminant;
		vector<int>     *jet_partonFlavour;
		vector<int>     *jet_hadronFlavour;
		vector<int>     *jet_isMatchedToGen;
		vector<bool>    *jet_isThight;
		vector<bool>    *isEEJJ;
		vector<bool>    *isEETT;
		vector<bool>    *isMMJJ;
		vector<bool>    *isMMTT;
		vector<bool>    *isEMJJ;
		vector<bool>    *isSignalRegion;
		vector<bool>    *isLowMllCR;
		vector<bool>    *isLowMlljjCR;
		vector<bool>    *isBB;
		vector<bool>    *isEE;
		vector<bool>    *isEB;
		vector<bool>    *passPreselections;
		vector<float>   *leadingLepton_e;
		vector<float>   *leadingLepton_pt;
		vector<float>   *leadingLepton_eta;
		vector<float>   *leadingLepton_phi;
		vector<float>   *leadingLepton_charge;
		vector<float>   *subLeadingLepton_e;
		vector<float>   *subLeadingLepton_pt;
		vector<float>   *subLeadingLepton_eta;
		vector<float>   *subLeadingLepton_phi;
		vector<float>   *subLeadingLepton_charge;
		vector<float>   *leadingJet_e;
		vector<float>   *leadingJet_pt;
		vector<float>   *leadingJet_eta;
		vector<float>   *leadingJet_phi;
		vector<int>     *leadingJet_isMatchedToGen;
		vector<bool>    *leadingJet_isThight;
		vector<float>   *subLeadingJet_e;
		vector<float>   *subLeadingJet_pt;
		vector<float>   *subLeadingJet_eta;
		vector<float>   *subLeadingJet_phi;
		vector<int>     *subLeadingJet_isMatchedToGen;
		vector<bool>    *subLeadingJet_isThight;
		vector<float>   *dRLeadLeptonLeadJet;
		vector<float>   *dRLeadLeptonSubLeadJet;
		vector<float>   *dRSubLeadLeptonLeadJet;
		vector<float>   *dRSubLeadLeptonSubLeadJet;
		vector<float>   *multiLeptonMultiJet_sumPt;
		vector<float>   *multiLeptonMultiJet_invMass;
		vector<float>   *diLepton_invMass;
		vector<float>   *diLepton_pt;
		vector<float>   *diJet_invMass;
		vector<float>   *diJetLeadingLepton_invMass;
		vector<float>   *diJetSubLeadingLepton_invMass;
		vector<bool>    *leadingEle_passHEEPId;
		vector<bool>    *leadingEle_passMediumId;
		vector<float>   *leadingEle_iso;
		vector<int>     *leadingEle_isMatchedToGen;
		vector<float>   *leadingEle_etaSC;
		vector<bool>    *leadingEle_isEcalDriven;
		vector<float>   *leadingEle_dEtaIn;
		vector<float>   *leadingEle_dPhiIn;
		vector<float>   *leadingEle_hOverE;
		vector<float>   *leadingEle_full5x5_r9;
		vector<float>   *leadingEle_full5x5_sigmaIetaIeta;
		vector<float>   *leadingEle_full5x5_E5x5;
		vector<float>   *leadingEle_full5x5_E1x5;
		vector<float>   *leadingEle_full5x5_E2x5;
		vector<float>   *leadingEle_EmHadDepth1Iso;
		vector<float>   *leadingEle_ptTracksIso;
		vector<int>     *leadingEle_innerLayerLostHits;
		vector<float>   *leadingEle_dxy;
		vector<float>   *leadingEle_dz;
		vector<float>   *leadingEle_eOverP;
		vector<float>   *leadingEle_ecalEnergy;
		vector<float>   *leadingEle_hcalOverEcal;
		vector<bool>    *subLeadingEle_passHEEPId;
		vector<bool>    *subLeadingEle_passMediumId;
		vector<float>   *subLeadingEle_iso;
		vector<int>     *subLeadingEle_isMatchedToGen;
		vector<float>   *subLeadingEle_etaSC;
		vector<bool>    *subLeadingEle_isEcalDriven;
		vector<float>   *subLeadingEle_dEtaIn;
		vector<float>   *subLeadingEle_dPhiIn;
		vector<float>   *subLeadingEle_hOverE;
		vector<float>   *subLeadingEle_full5x5_r9;
		vector<float>   *subLeadingEle_full5x5_sigmaIetaIeta;
		vector<float>   *subLeadingEle_full5x5_E5x5;
		vector<float>   *subLeadingEle_full5x5_E1x5;
		vector<float>   *subLeadingEle_full5x5_E2x5;
		vector<float>   *subLeadingEle_EmHadDepth1Iso;
		vector<float>   *subLeadingEle_ptTracksIso;
		vector<int>     *subLeadingEle_innerLayerLostHits;
		vector<float>   *subLeadingEle_dxy; 
		vector<float>   *subLeadingEle_dz;
		vector<float>   *subLeadingEle_eOverP;
		vector<float>   *subLeadingEle_ecalEnergy;
		vector<float>   *subLeadingEle_hcalOverEcal;
		vector<float>   *leadingMuon_iso;
		vector<float>   *leadingMuon_PFiso;
		vector<bool>    *leadingMuon_isHighPt;
		vector<int>     *leadingMuon_isMatchedToGen;
		vector<float>   *leadingMuon_dz;
		vector<float>   *leadingMuon_dxy;
		vector<float>   *leadingMuon_RochCor;
		vector<float>   *subLeadingMuon_iso;
		vector<float>   *subLeadingMuon_PFiso;
		vector<bool>    *subLeadingMuon_isHighPt;
		vector<int>     *subLeadingMuon_isMatchedToGen;
		vector<float>   *subLeadingMuon_dz;
		vector<float>   *subLeadingMuon_dxy;
		vector<float>   *subLeadingMuon_RochCor;

	// List of branches
		TBranch        *b_run;   //!
		TBranch        *b_event;   //!
		TBranch        *b_lumi;   //!
		TBranch        *b_rho;   //!
		TBranch        *b_weight;   //!
		TBranch        *b_puweight;   //!
		TBranch        *b_nvtx;   //!
		TBranch        *b_npu;   //!
		TBranch        *b_vtx_x;   //!
		TBranch        *b_vtx_y;   //!
		TBranch        *b_vtx_z;   //!
		TBranch        *b_ele_e;   //!
		TBranch        *b_ele_pt;   //!
		TBranch        *b_ele_eta;   //!
		TBranch        *b_ele_phi;   //!
		TBranch        *b_ele_passHEEPId;   //!
		TBranch        *b_ele_passMediumId;   //!
		TBranch        *b_ele_iso;   //!
		TBranch        *b_ele_dz;   //!
		TBranch        *b_ele_d0;   //!
		TBranch        *b_ele_isMatchedToGen;   //!
		TBranch        *b_ele_charge;   //!
		TBranch        *b_ele_etaSC;   //!
		TBranch        *b_ele_isEcalDriven;   //!
		TBranch        *b_ele_dEtaIn;   //!
		TBranch        *b_ele_dPhiIn;   //!
		TBranch        *b_ele_hOverE;   //!
		TBranch        *b_ele_full5x5_r9;   //!
		TBranch        *b_ele_full5x5_sigmaIetaIeta;   //!
		TBranch        *b_ele_full5x5_E5x5;   //!
		TBranch        *b_ele_full5x5_E1x5;   //!
		TBranch        *b_ele_full5x5_E2x5;   //!
		TBranch        *b_ele_EmHadDepth1Iso;   //!
		TBranch        *b_ele_ptTracksIso;   //!
		TBranch        *b_ele_innerLayerLostHits;   //!
		TBranch        *b_ele_dxy;   //!
		TBranch        *b_ele_eOverP;   //!
		TBranch        *b_ele_ecalEnergy;   //!
		TBranch        *b_ele_hcalOverEcal;   //!
		TBranch        *b_mu_e;   //!
		TBranch        *b_mu_pt;   //!
		TBranch        *b_mu_eta;   //!
		TBranch        *b_mu_phi;   //!
		TBranch        *b_mu_iso;   //!
		TBranch        *b_mu_PFiso;   //!
		TBranch        *b_mu_isHighPt;   //!
		TBranch        *b_mu_isMatchedToGen;   //!
		TBranch        *b_mu_charge;   //!
		TBranch        *b_mu_dz;   //!
		TBranch        *b_mu_dxy;   //!
		TBranch        *b_mu_RochCor;   //!
		TBranch        *b_jet_e;   //!
		TBranch        *b_jet_pt;   //!
		TBranch        *b_jet_eta;   //!
		TBranch        *b_jet_phi;   //!
		TBranch        *b_jet_bdiscriminant;   //!
		TBranch        *b_jet_partonFlavour;   //!
		TBranch        *b_jet_hadronFlavour;   //!
		TBranch        *b_jet_isMatchedToGen;   //!
		TBranch        *b_jet_isThight;   //!
		TBranch        *b_isEEJJ;   //!
		TBranch        *b_isEETT;   //!
		TBranch        *b_isMMJJ;   //!
		TBranch        *b_isMMTT;   //!
		TBranch        *b_isEMJJ;   //!
		TBranch        *b_isSignalRegion;   //!
		TBranch        *b_isLowMllCR;   //!
		TBranch        *b_isLowMlljjCR;   //!
		TBranch        *b_isBB;   //!
		TBranch        *b_isEE;   //!
		TBranch        *b_isEB;   //!
		TBranch        *b_passPreselections;   //!
		TBranch        *b_leadingLepton_e;   //!
		TBranch        *b_leadingLepton_pt;   //!
		TBranch        *b_leadingLepton_eta;   //!
		TBranch        *b_leadingLepton_phi;   //!
		TBranch        *b_leadingLepton_charge;   //!
		TBranch        *b_subLeadingLepton_e;   //!
		TBranch        *b_subLeadingLepton_pt;   //!
		TBranch        *b_subLeadingLepton_eta;   //!
		TBranch        *b_subLeadingLepton_phi;   //!
		TBranch        *b_subLeadingLepton_charge;   //!
		TBranch        *b_leadingJet_e;   //!
		TBranch        *b_leadingJet_pt;   //!
		TBranch        *b_leadingJet_eta;   //!
		TBranch        *b_leadingJet_phi;   //!
		TBranch        *b_leadingJet_isMatchedToGen;   //!
		TBranch        *b_leadingJet_isThight;   //!
		TBranch        *b_subLeadingJet_e;   //!
		TBranch        *b_subLeadingJet_pt;   //!
		TBranch        *b_subLeadingJet_eta;   //!
		TBranch        *b_subLeadingJet_phi;   //!
		TBranch        *b_subLeadingJet_isMatchedToGen;   //!
		TBranch        *b_subLeadingJet_isThight;   //!
		TBranch        *b_dRLeadLeptonLeadJet;   //!
		TBranch        *b_dRLeadLeptonSubLeadJet;   //!
		TBranch        *b_dRSubLeadLeptonLeadJet;   //!
		TBranch        *b_dRSubLeadLeptonSubLeadJet;   //!
		TBranch        *b_multiLeptonMultiJet_sumPt;   //!
		TBranch        *b_multiLeptonMultiJet_invMass;   //!
		TBranch        *b_diLepton_invMass;   //!
		TBranch        *b_diLepton_pt;   //!
		TBranch        *b_diJet_invMass;   //!
		TBranch        *b_diJetLeadingLepton_invMass;   //!
		TBranch        *b_diJetSubLeadingLepton_invMass;   //!
		TBranch        *b_leadingEle_passHEEPId;   //!
		TBranch        *b_leadingEle_passMediumId;   //!
		TBranch        *b_leadingEle_iso;   //!
		TBranch        *b_leadingEle_isMatchedToGen;   //!
		TBranch        *b_leadingEle_etaSC;   //!
		TBranch        *b_leadingEle_isEcalDriven;   //!
		TBranch        *b_leadingEle_dEtaIn;   //!
		TBranch        *b_leadingEle_dPhiIn;   //!
		TBranch        *b_leadingEle_hOverE;   //!
		TBranch        *b_leadingEle_full5x5_r9;   //!
		TBranch        *b_leadingEle_full5x5_sigmaIetaIeta;   //!
		TBranch        *b_leadingEle_full5x5_E5x5;   //!
		TBranch        *b_leadingEle_full5x5_E1x5;   //!
		TBranch        *b_leadingEle_full5x5_E2x5;   //!
		TBranch        *b_leadingEle_EmHadDepth1Iso;   //!
		TBranch        *b_leadingEle_ptTracksIso;   //!
		TBranch        *b_leadingEle_innerLayerLostHits;   //!
		TBranch        *b_leadingEle_dxy;   //!
		TBranch        *b_leadingEle_dz;   //!
		TBranch        *b_leadingEle_eOverP;   //!
		TBranch        *b_leadingEle_ecalEnergy;   //!
		TBranch        *b_leadingEle_hcalOverEcal;   //!
		TBranch        *b_subLeadingEle_passHEEPId;   //!
		TBranch        *b_subLeadingEle_passMediumId;   //!
		TBranch        *b_subLeadingEle_iso;   //!
		TBranch        *b_subLeadingEle_isMatchedToGen;   //!
		TBranch        *b_subLeadingEle_etaSC;   //!
		TBranch        *b_subLeadingEle_isEcalDriven;   //!
		TBranch        *b_subLeadingEle_dEtaIn;   //!
		TBranch        *b_subLeadingEle_dPhiIn;   //!
		TBranch        *b_subLeadingEle_hOverE;   //!
		TBranch        *b_subLeadingEle_full5x5_r9;   //!
		TBranch        *b_subLeadingEle_full5x5_sigmaIetaIeta;   //!
		TBranch        *b_subLeadingEle_full5x5_E5x5;   //!
		TBranch        *b_subLeadingEle_full5x5_E1x5;   //!
		TBranch        *b_subLeadingEle_full5x5_E2x5;   //!
		TBranch        *b_subLeadingEle_EmHadDepth1Iso;   //!
		TBranch        *b_subLeadingEle_ptTracksIso;   //!
		TBranch        *b_subLeadingEle_innerLayerLostHits;   //!
		TBranch        *b_subLeadingEle_dxy;   //!
		TBranch        *b_subLeadingEle_dz;   //!
		TBranch        *b_subLeadingEle_eOverP;   //!
		TBranch        *b_subLeadingEle_ecalEnergy;   //!
		TBranch        *b_subLeadingEle_hcalOverEcal;   //!
		TBranch        *b_leadingMuon_iso;   //!
		TBranch        *b_leadingMuon_PFiso;   //!
		TBranch        *b_leadingMuon_isHighPt;   //!
		TBranch        *b_leadingMuon_isTight;   //!
		TBranch        *b_leadingMuon_isMedium;   //!
		TBranch        *b_leadingMuon_isLoose;   //!
		TBranch        *b_leadingMuon_isMatchedToGen;   //!
		TBranch        *b_leadingMuon_dz;   //!
		TBranch        *b_leadingMuon_dxy;   //!
		TBranch        *b_leadingMuon_RochCor;   //!
		TBranch        *b_subLeadingMuon_iso;   //!
		TBranch        *b_subLeadingMuon_PFiso;   //!
		TBranch        *b_subLeadingMuon_isHighPt;   //!
		TBranch        *b_subLeadingMuon_isTight;   //!
		TBranch        *b_subLeadingMuon_isMedium;   //!
		TBranch        *b_subLeadingMuon_isLoose;   //!
		TBranch        *b_subLeadingMuon_isMatchedToGen;   //!
		TBranch        *b_subLeadingMuon_dz;   //!
		TBranch        *b_subLeadingMuon_dxy;   //!
		TBranch        *b_subLeadingMuon_RochCor;   //!


		static const int nEtaRegions = 3;
		static const int nEtaMassRegions = 5;

		static const int nObjects = 16;
		static const int nEle = 6;

		static const int nRegions = 3;
		static const int nZregions = 6; //4;


	//definition of histos
		TH1D *rho_histo;
		TH1D *nvtx_histo;

		TH1D *vtx_x_histo;
		TH1D *vtx_y_histo;
		TH1D *vtx_z_histo;

		TH1D *muon_dxy_histo;

		TH1D *nLeadingLeptons_histo;
		TH1D *nSubLeadingLeptons_histo;

		TH1D *pt_histo[nObjects][nRegions];
		TH1D *eta_histo[nObjects][nRegions];
		TH1D *phi_histo[nObjects][nRegions];

		TH1D *etaSC_histo[nEle][nEtaRegions];
		TH1D *isEcalDriven_histo[nEle][nEtaRegions];  
		TH1D *dEtaIn_histo[nEle][nEtaRegions];
		TH1D *dPhiIn_histo[nEle][nEtaRegions];
		TH1D *hOverE_histo[nEle][nEtaRegions];
		TH1D *full5x5_r9_histo[nEle][nEtaRegions];
		TH1D *full5x5_sigmaIetaIeta_histo[nEle][nEtaRegions];
		TH1D *full5x5_E5x5_histo[nEle][nEtaRegions];
		TH1D *full5x5_E1x5_histo[nEle][nEtaRegions];
		TH1D *full5x5_E2x5_histo[nEle][nEtaRegions];
		TH1D *full5x5_E2x5_Over_E5x5_histo[nEle][nEtaRegions];
		TH1D *full5x5_E1x5_Over_E5x5_histo[nEle][nEtaRegions];
		TH1D *EmHadDepth1Iso_histo[nEle][nEtaRegions];
		TH1D *ptTracksIso_histo[nEle][nEtaRegions];
		TH1D *innerLayerLostHits_histo[nEle][nEtaRegions];
		TH1D *dxy_histo[nEle][nEtaRegions];
		TH1D *eOverP[nEle][nEtaRegions];

		TH1D *mass_dldj_histo[nRegions][nEtaMassRegions];
		TH1D *mass_dl_histo[nRegions][nEtaMassRegions];
		TH1D *mass_dj_histo[nRegions][nEtaMassRegions];
		TH1D *mass_djLl_histo[nRegions][nEtaMassRegions];
		TH1D *mass_djSLl_histo[nRegions][nEtaMassRegions];
		TH1D *pt_dl_histo[nRegions];

		TH1D *Zmass_histo[nZregions][nEtaMassRegions];
		TH1D *Zpt_histo[nZregions];


		makePlots_multiLeptonMultiJet(TString filename_, TString outputdir_, bool MC_, bool signalEE_, bool signalMuMu_, bool eMuSideband_, bool TnPEE_, bool TnPMuMu_, int nToDivideWeight_);

		void     Init();  
		Int_t    GetEntry(Long64_t entry);
		Long64_t LoadTree(Long64_t entry);

		TH1D* newTH1D(string, string, string, int, double, double);
		TH1D* newTH1DvarBinSize(string, string, string, int, const double *);
		TH2D* newTH2D(string, string, string, string, int, double, double, int, double, double);
		TProfile* newTProfile(string, string, string, string, int, double, double, double, double);
		void SetHistos();

		void getElectrons(vector<eleStruct>& electrons);//, vector<eleStruct>& goodElectrons);
		void getLeadingElectrons(vector<eleStruct>& leadingElectrons);
		void getSubLeadingElectrons(vector<eleStruct>& subLeadingElectrons);

		void getMuons(vector<muonStruct>& muons, vector<muonStruct>& muonsRochCorr);//, vector<muonStruct>& goodMuons, vector<muonStruct>& goodMuonsRochCorr);
		void getLeadingMuons(vector<muonStruct>& leadingMuons);
		void getSubLeadingMuons(vector<muonStruct>& subLeadingMuons);

		void getJets(vector<jetStruct>& jets);//, vector<jetStruct>& goodJets);
		void getLeadingJets(vector<jetStruct>& leadingJets);
		void getSubLeadingJets(vector<jetStruct>& subLeadingJets);

		void getMultiLeptonMultiJets(vector<multiLeptonMultiJetStruct>& multiLeptonMultiJets);

		void doEleDistributionsPlots(vector<eleStruct>& electrons, const int eleIdx, const int variableIdx, const int regionIdx);
		void doMuonDistributionsPlots(vector<muonStruct>& muons, const int muIdx, const int variableIdx, const int regionIdx);
		void doJetsDistributionsPlots(vector<jetStruct>& jets, const int jetIdx, const int variableIdx, const int regionIdx);
		
		void doElePlots(vector<eleStruct>& electrons, const int eleIdx, const int variableIdx, const int etaIdx);

		void doMassPlots(vector<multiLeptonMultiJetStruct>& multiLeptonMultiJets, const int dldjIdx, const int regionIdx, const int etaMassIdx);
		void doZPlots(vector<multiLeptonMultiJetStruct>& multiLeptonMultiJets, const int dldjIdx, const int variableIdx);
		void doZPlots(TLorentzVector Llepton, TLorentzVector SLlepton, const int variableIdx);

		void Loop();
		void saveHistosAndOutputFile(TString& outputdir);

		TString filename, outputdir; 
		bool MC, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu;
		int nToDivideWeight;		
		vector<TH1*> listOfHistograms;

		unsigned int nEvents=0, nEventsWithAtLeast1DLDJ=0, nEventsWithRightLeptonPair=0, nEventsWithDLDJpassingPreselections=0, nEventsWithDLDJpassingSelections=0;
		unsigned int nEventsWithDLDJpassingSelectionsInSignalRegion=0, nEventsWithDLDJpassingSelectionsInLowMlljjCR=0, nEventsWithDLDJpassingSelectionsInLowMllCR=0, nEventsWithDLDJpassingSelectionsInFlavourSidebandCR=0;
		unsigned int nEventsWithLeadingElePassingSelections=0, nEventsWithSubLeadingElePassingSelections=0;
		unsigned int nEventsWithLeadingMuonPassingSelections=0, nEventsWithSubLeadingMuonPassingSelections=0;
		unsigned int nEventsWithLeadingJetPassingSelections=0, nEventsWithSubLeadingJetPassingSelections=0;
		
		unsigned int nEventsWithEEJJpassingLoosePreselections=0, nEventsWithMMJJpassingLoosePreselections=0, nEventsWithEEJJpassingLoosePreselectionsAndCharge=0, nEventsWithMMJJpassingLoosePreselectionsAndCharge=0;
		unsigned int nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId=0, nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt=0;
		unsigned int nEventsWith2elePassingLoosePreselections=0, nEventsWith2elePassingLoosePreselectionsAndCharge=0, nEventsWith2elePassingLoosePreselectionsAndChargeAndHeepId=0; 
		unsigned int nEventsWith2muonsPassingLoosePreselections=0, nEventsWith2muonsPassingLoosePreselectionsAndCharge=0, nEventsWith2muonsPassingLoosePreselectionsAndChargeAndHighPt=0; 
		
		unsigned int nLeadingElePassingPreselections=0, nLeadingElePassingIsEcalDriven=0, nLeadingElePassingdEtaIn=0, nLeadingElePassingdPhiIn=0, nLeadingElePassingE2x5OverE5x5=0, nLeadingElePassingEmHadDepth1Iso=0, nLeadingElePassingMissingHits=0, nLeadingElePassingDxy=0, nLeadingElePassingHeepId=0;
		unsigned int nLeadingElePassingPreselectionsInEB=0, nLeadingEleInEBPassingIsEcalDriven=0, nLeadingEleInEBPassingdEtaIn=0, nLeadingEleInEBPassingdPhiIn=0, nLeadingEleInEBPassingE2x5OverE5x5=0, nLeadingEleInEBPassingEmHadDepth1Iso=0, nLeadingEleInEBPassingMissingHits=0, nLeadingEleInEBPassingDxy=0, nLeadingElePassingHeepIdInEB=0;
		unsigned int nLeadingElePassingPreselectionsInEE=0, nLeadingEleInEEPassingIsEcalDriven=0, nLeadingEleInEEPassingdEtaIn=0, nLeadingEleInEEPassingdPhiIn=0, nLeadingEleInEEPassingE2x5OverE5x5=0, nLeadingEleInEEPassingEmHadDepth1Iso=0, nLeadingEleInEEPassingMissingHits=0, nLeadingEleInEEPassingDxy=0, nLeadingElePassingHeepIdInEE=0;
		unsigned int nSubLeadingElePassingPreselections=0, nSubLeadingElePassingIsEcalDriven=0, nSubLeadingElePassingdEtaIn=0, nSubLeadingElePassingdPhiIn=0, nSubLeadingElePassingE2x5OverE5x5=0, nSubLeadingElePassingEmHadDepth1Iso=0, nSubLeadingElePassingMissingHits=0, nSubLeadingElePassingDxy=0, nSubLeadingElePassingHeepId=0;
		unsigned int nSubLeadingElePassingPreselectionsInEB=0, nSubLeadingEleInEBPassingIsEcalDriven=0, nSubLeadingEleInEBPassingdEtaIn=0, nSubLeadingEleInEBPassingdPhiIn=0, nSubLeadingEleInEBPassingE2x5OverE5x5=0, nSubLeadingEleInEBPassingEmHadDepth1Iso=0, nSubLeadingEleInEBPassingMissingHits=0, nSubLeadingEleInEBPassingDxy=0, nSubLeadingElePassingHeepIdInEB=0;
		unsigned int nSubLeadingElePassingPreselectionsInEE=0, nSubLeadingEleInEEPassingIsEcalDriven=0, nSubLeadingEleInEEPassingdEtaIn=0, nSubLeadingEleInEEPassingdPhiIn=0, nSubLeadingEleInEEPassingE2x5OverE5x5=0, nSubLeadingEleInEEPassingEmHadDepth1Iso=0, nSubLeadingEleInEEPassingMissingHits=0, nSubLeadingEleInEEPassingDxy=0, nSubLeadingElePassingHeepIdInEE=0;

		unsigned int nEventsWithLandSLleptonsRight=0,nEventsWithLandSLleptonsInverted=0,nEventsWithLandSLjetsRight=0,nEventsWithLandSLjetsInverted=0;
		unsigned int nEventsFillingZtoEE=0,nEventsFillingZtoEE_MLMJ=0;

		float w=1, sumWeights = 0; 

		bool saveHEEPvariables_ = false;

		string etaName[nEtaRegions] = {"", "_EB", "_EE"};
		string etaMassName[nEtaMassRegions] = {"", "_EB-EB", "_EE-EE", "_EB-EE", "_noEB-EB"};

		string objName[nObjects] = {"leadingLepton", "subLeadingLepton", "leadingEle", 
		"subLeadingEle", "leadingMu", "subLeadingMu", "leadingJet", "subLeadingJet", 
		"ele", "goodEle", "muon", "goodMuon", "jet", "goodJet", "muonRochCor", "goodMuonRochCor"};

		string eleName[nEle] = {"leadingEle", "subLeadingEle", "goodLeadingEle", 
		"goodSubLeadingEle", "ele", "goodEle"};

		string triggerName;

		int numberRegions;
		string regionName[nRegions];
		string signalRegionNames[3] = {"_signalRegion","_lowMlljjCR","_lowMllCR"};
		string eMuSidebandName[1] = {"_eMuSidebandCR"};
		string TnPCRName[1] = {"_TnPCR"};

		string Zname[nZregions] = {"toEE_electrons", "toMuMu_muons", "toEE_MLMJelectrons", "toMuMu_MLMJmuons", "toEE_electronsAndJets", "toMuMu_muonsAndJets"};

};

#endif

