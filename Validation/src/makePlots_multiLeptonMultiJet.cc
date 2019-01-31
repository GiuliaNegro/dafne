#include "dafne/Validation/interface/makePlots_multiLeptonMultiJet.h"
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <limits>
#include <iostream>
#include <sstream>
// #include "tdrstyle.C"
// #include "CMS_lumi.C"

using namespace std;



// **************** 
void makePlots_multiLeptonMultiJet::Init(){

	// Set object pointer
	vtx_x = 0;
	vtx_y = 0;
	vtx_z = 0;
	ele_e = 0;
	ele_pt = 0;
	ele_eta = 0;
	ele_phi = 0;
	ele_passHEEPId = 0;
	ele_passMediumId = 0;
	ele_iso = 0;
	ele_dz = 0;
	ele_d0 = 0;
	ele_isMatchedToGen = 0;
	ele_charge = 0;
	if (saveHEEPvariables_) {
		ele_etaSC = 0;
		ele_isEcalDriven = 0;
		ele_dEtaIn = 0;
		ele_dPhiIn = 0;
		ele_hOverE = 0;
		ele_full5x5_r9 = 0;
		ele_full5x5_sigmaIetaIeta = 0;
		ele_full5x5_E5x5 = 0;
		ele_full5x5_E1x5 = 0;
		ele_full5x5_E2x5 = 0;
		ele_EmHadDepth1Iso = 0;
		ele_ptTracksIso = 0;
		ele_innerLayerLostHits = 0;
		ele_dxy = 0;
		ele_eOverP = 0;
		ele_ecalEnergy = 0;
		ele_hcalOverEcal = 0;
	}
	mu_e = 0;
	mu_pt = 0;
	mu_eta = 0;
	mu_phi = 0;
	mu_iso = 0;
	mu_PFiso = 0;
	mu_isHighPt = 0;
	mu_isMatchedToGen = 0;
	mu_charge = 0;
	mu_dz = 0;
	mu_dxy = 0;
	mu_RochCor = 0;
	jet_e = 0;
	jet_pt = 0;
	jet_eta = 0;
	jet_phi = 0;
	jet_isMatchedToGen = 0;
	jet_isThight = 0;
	isEEJJ = 0;
	isEETT = 0;
	isMMJJ = 0;
	isMMTT = 0;
	isEMJJ = 0;
	isSignalRegion = 0;
	isLowMllCR = 0;
	isLowMlljjCR = 0;
	isBB = 0;
	isEE = 0;
	isEB = 0;
	passPreselections = 0;
	leadingLepton_e = 0;
	leadingLepton_pt = 0;
	leadingLepton_eta = 0;
	leadingLepton_phi = 0;
	leadingLepton_charge = 0;
	subLeadingLepton_e = 0;
	subLeadingLepton_pt = 0;
	subLeadingLepton_eta = 0;
	subLeadingLepton_phi = 0;
	subLeadingLepton_charge = 0;
	leadingJet_e = 0;
	leadingJet_pt = 0;
	leadingJet_eta = 0;
	leadingJet_phi = 0;
	leadingJet_isMatchedToGen = 0;
	leadingJet_isThight = 0;
	subLeadingJet_e = 0;
	subLeadingJet_pt = 0;
	subLeadingJet_eta = 0;
	subLeadingJet_phi = 0;
	subLeadingJet_isMatchedToGen = 0;
	subLeadingJet_isThight = 0;
	dRLeadLeptonLeadJet = 0;
	dRLeadLeptonSubLeadJet = 0;
	dRSubLeadLeptonLeadJet = 0;
	dRSubLeadLeptonSubLeadJet = 0;
	multiLeptonMultiJet_sumPt = 0;
	multiLeptonMultiJet_invMass = 0;
	diLepton_invMass = 0;
	diLepton_pt = 0;
	diJet_invMass = 0;
	diJetLeadingLepton_invMass = 0;
	diJetSubLeadingLepton_invMass = 0;
	leadingEle_passHEEPId = 0;
	leadingEle_passMediumId = 0;
	leadingEle_iso = 0;
	leadingEle_isMatchedToGen = 0;
	if (saveHEEPvariables_) {
		leadingEle_etaSC = 0;
		leadingEle_isEcalDriven = 0;
		leadingEle_dEtaIn = 0;
		leadingEle_dPhiIn = 0;
		leadingEle_hOverE = 0;
		leadingEle_full5x5_r9 = 0;
		leadingEle_full5x5_sigmaIetaIeta = 0;
		leadingEle_full5x5_E5x5 = 0;
		leadingEle_full5x5_E1x5 = 0;
		leadingEle_full5x5_E2x5 = 0;
		leadingEle_EmHadDepth1Iso = 0;
		leadingEle_ptTracksIso = 0;
		leadingEle_innerLayerLostHits = 0;
		leadingEle_dxy = 0;
		leadingEle_dz = 0;
		leadingEle_eOverP = 0;
		leadingEle_ecalEnergy = 0;
		leadingEle_hcalOverEcal = 0;
	}
	subLeadingEle_passHEEPId = 0;
	subLeadingEle_passMediumId = 0;
	subLeadingEle_iso = 0;
	subLeadingEle_isMatchedToGen = 0;
	if (saveHEEPvariables_) {
		subLeadingEle_etaSC = 0;
		subLeadingEle_isEcalDriven = 0;
		subLeadingEle_dEtaIn = 0;
		subLeadingEle_dPhiIn = 0;
		subLeadingEle_hOverE = 0;
		subLeadingEle_full5x5_r9 = 0;
		subLeadingEle_full5x5_sigmaIetaIeta = 0;
		subLeadingEle_full5x5_E5x5 = 0;
		subLeadingEle_full5x5_E1x5 = 0;
		subLeadingEle_full5x5_E2x5 = 0;
		subLeadingEle_EmHadDepth1Iso = 0;
		subLeadingEle_ptTracksIso = 0;
		subLeadingEle_innerLayerLostHits = 0;
		subLeadingEle_dxy = 0;
		subLeadingEle_dz = 0;
		subLeadingEle_eOverP = 0;
		subLeadingEle_ecalEnergy = 0;
		subLeadingEle_hcalOverEcal = 0;
	}
	leadingMuon_iso = 0;
	leadingMuon_PFiso = 0;
	leadingMuon_isHighPt = 0;
	leadingMuon_isMatchedToGen = 0;
	leadingMuon_dz = 0;
	leadingMuon_dxy = 0;
	leadingMuon_RochCor = 0;
	subLeadingMuon_iso = 0;
	subLeadingMuon_PFiso = 0;
	subLeadingMuon_isHighPt = 0;
	subLeadingMuon_isMatchedToGen = 0;
	subLeadingMuon_dz = 0;
	subLeadingMuon_dxy = 0;
	subLeadingMuon_RochCor = 0;

	// Set branch addresses and branch pointers
	// if (!tree) return;
	// fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("run", &run, &b_run);
	fChain->SetBranchAddress("event", &event, &b_event);
	fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
	fChain->SetBranchAddress("rho", &rho, &b_rho);
	fChain->SetBranchAddress("weight", &weight, &b_weight);
	fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
	fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
	fChain->SetBranchAddress("npu", &npu, &b_npu);
	fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
	fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
	fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
	fChain->SetBranchAddress("ele_e", &ele_e, &b_ele_e);
	fChain->SetBranchAddress("ele_pt", &ele_pt, &b_ele_pt);
	fChain->SetBranchAddress("ele_eta", &ele_eta, &b_ele_eta);
	fChain->SetBranchAddress("ele_phi", &ele_phi, &b_ele_phi);
	fChain->SetBranchAddress("ele_passHEEPId", &ele_passHEEPId, &b_ele_passHEEPId);
	fChain->SetBranchAddress("ele_passMediumId", &ele_passMediumId, &b_ele_passMediumId);
	fChain->SetBranchAddress("ele_iso", &ele_iso, &b_ele_iso);
	fChain->SetBranchAddress("ele_dz", &ele_dz, &b_ele_dz);
	fChain->SetBranchAddress("ele_d0", &ele_d0, &b_ele_d0);
	fChain->SetBranchAddress("ele_isMatchedToGen", &ele_isMatchedToGen, &b_ele_isMatchedToGen);
	fChain->SetBranchAddress("ele_charge", &ele_charge, &b_ele_charge);
	if (saveHEEPvariables_) {
		fChain->SetBranchAddress("ele_etaSC", &ele_etaSC, &b_ele_etaSC);
		fChain->SetBranchAddress("ele_isEcalDriven", &ele_isEcalDriven, &b_ele_isEcalDriven);
		fChain->SetBranchAddress("ele_dEtaIn", &ele_dEtaIn, &b_ele_dEtaIn);
		fChain->SetBranchAddress("ele_dPhiIn", &ele_dPhiIn, &b_ele_dPhiIn);
		fChain->SetBranchAddress("ele_hOverE", &ele_hOverE, &b_ele_hOverE);
		fChain->SetBranchAddress("ele_full5x5_r9", &ele_full5x5_r9, &b_ele_full5x5_r9);
		fChain->SetBranchAddress("ele_full5x5_sigmaIetaIeta", &ele_full5x5_sigmaIetaIeta, &b_ele_full5x5_sigmaIetaIeta);
		fChain->SetBranchAddress("ele_full5x5_E5x5", &ele_full5x5_E5x5, &b_ele_full5x5_E5x5);
		fChain->SetBranchAddress("ele_full5x5_E1x5", &ele_full5x5_E1x5, &b_ele_full5x5_E1x5);
		fChain->SetBranchAddress("ele_full5x5_E2x5", &ele_full5x5_E2x5, &b_ele_full5x5_E2x5);
		fChain->SetBranchAddress("ele_EmHadDepth1Iso", &ele_EmHadDepth1Iso, &b_ele_EmHadDepth1Iso);
		fChain->SetBranchAddress("ele_ptTracksIso", &ele_ptTracksIso, &b_ele_ptTracksIso);
		fChain->SetBranchAddress("ele_innerLayerLostHits", &ele_innerLayerLostHits, &b_ele_innerLayerLostHits);
		fChain->SetBranchAddress("ele_dxy", &ele_dxy, &b_ele_dxy);
		fChain->SetBranchAddress("ele_eOverP", &ele_eOverP, &b_ele_eOverP);
		fChain->SetBranchAddress("ele_ecalEnergy", &ele_ecalEnergy, &b_ele_ecalEnergy);
		fChain->SetBranchAddress("ele_hcalOverEcal", &ele_hcalOverEcal, &b_ele_hcalOverEcal);
	}
	fChain->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
	fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
	fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
	fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
	fChain->SetBranchAddress("mu_iso", &mu_iso, &b_mu_iso);
	fChain->SetBranchAddress("mu_PFiso", &mu_PFiso, &b_mu_PFiso);
	fChain->SetBranchAddress("mu_isHighPt", &mu_isHighPt, &b_mu_isHighPt);
	fChain->SetBranchAddress("mu_isMatchedToGen", &mu_isMatchedToGen, &b_mu_isMatchedToGen);
	fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
	fChain->SetBranchAddress("mu_dz", &mu_dz, &b_mu_dz);
	fChain->SetBranchAddress("mu_dxy", &mu_dxy, &b_mu_dxy);
	fChain->SetBranchAddress("mu_RochCor", &mu_RochCor, &b_mu_RochCor);
	fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
	fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
	fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
	fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
	fChain->SetBranchAddress("jet_isMatchedToGen", &jet_isMatchedToGen, &b_jet_isMatchedToGen);
	fChain->SetBranchAddress("jet_isThight", &jet_isThight, &b_jet_isThight);
	fChain->SetBranchAddress("isEEJJ", &isEEJJ, &b_isEEJJ);
	fChain->SetBranchAddress("isEETT", &isEETT, &b_isEETT);
	fChain->SetBranchAddress("isMMJJ", &isMMJJ, &b_isMMJJ);
	fChain->SetBranchAddress("isMMTT", &isMMTT, &b_isMMTT);
	fChain->SetBranchAddress("isEMJJ", &isEMJJ, &b_isEMJJ);
	fChain->SetBranchAddress("isSignalRegion", &isSignalRegion, &b_isSignalRegion);
	fChain->SetBranchAddress("isLowMllCR", &isLowMllCR, &b_isLowMllCR);
	fChain->SetBranchAddress("isLowMlljjCR", &isLowMlljjCR, &b_isLowMlljjCR);
	fChain->SetBranchAddress("isBB", &isBB, &b_isBB);
	fChain->SetBranchAddress("isEE", &isEE, &b_isEE);
	fChain->SetBranchAddress("isEB", &isEB, &b_isEB);
	fChain->SetBranchAddress("passPreselections", &passPreselections, &b_passPreselections);
	fChain->SetBranchAddress("leadingLepton_e", &leadingLepton_e, &b_leadingLepton_e);
	fChain->SetBranchAddress("leadingLepton_pt", &leadingLepton_pt, &b_leadingLepton_pt);
	fChain->SetBranchAddress("leadingLepton_eta", &leadingLepton_eta, &b_leadingLepton_eta);
	fChain->SetBranchAddress("leadingLepton_phi", &leadingLepton_phi, &b_leadingLepton_phi);
	fChain->SetBranchAddress("leadingLepton_charge", &leadingLepton_charge, &b_leadingLepton_charge);
	fChain->SetBranchAddress("subLeadingLepton_e", &subLeadingLepton_e, &b_subLeadingLepton_e);
	fChain->SetBranchAddress("subLeadingLepton_pt", &subLeadingLepton_pt, &b_subLeadingLepton_pt);
	fChain->SetBranchAddress("subLeadingLepton_eta", &subLeadingLepton_eta, &b_subLeadingLepton_eta);
	fChain->SetBranchAddress("subLeadingLepton_phi", &subLeadingLepton_phi, &b_subLeadingLepton_phi);
	fChain->SetBranchAddress("subLeadingLepton_charge", &subLeadingLepton_charge, &b_subLeadingLepton_charge);
	fChain->SetBranchAddress("leadingJet_e", &leadingJet_e, &b_leadingJet_e);
	fChain->SetBranchAddress("leadingJet_pt", &leadingJet_pt, &b_leadingJet_pt);
	fChain->SetBranchAddress("leadingJet_eta", &leadingJet_eta, &b_leadingJet_eta);
	fChain->SetBranchAddress("leadingJet_phi", &leadingJet_phi, &b_leadingJet_phi);
	fChain->SetBranchAddress("leadingJet_isMatchedToGen", &leadingJet_isMatchedToGen, &b_leadingJet_isMatchedToGen);
	fChain->SetBranchAddress("leadingJet_isThight", &leadingJet_isThight, &b_leadingJet_isThight);
	fChain->SetBranchAddress("subLeadingJet_e", &subLeadingJet_e, &b_subLeadingJet_e);
	fChain->SetBranchAddress("subLeadingJet_pt", &subLeadingJet_pt, &b_subLeadingJet_pt);
	fChain->SetBranchAddress("subLeadingJet_eta", &subLeadingJet_eta, &b_subLeadingJet_eta);
	fChain->SetBranchAddress("subLeadingJet_phi", &subLeadingJet_phi, &b_subLeadingJet_phi);
	fChain->SetBranchAddress("subLeadingJet_isMatchedToGen", &subLeadingJet_isMatchedToGen, &b_subLeadingJet_isMatchedToGen);
	fChain->SetBranchAddress("subLeadingJet_isThight", &subLeadingJet_isThight, &b_subLeadingJet_isThight);
	fChain->SetBranchAddress("dRLeadLeptonLeadJet", &dRLeadLeptonLeadJet, &b_dRLeadLeptonLeadJet);
	fChain->SetBranchAddress("dRLeadLeptonSubLeadJet", &dRLeadLeptonSubLeadJet, &b_dRLeadLeptonSubLeadJet);
	fChain->SetBranchAddress("dRSubLeadLeptonLeadJet", &dRSubLeadLeptonLeadJet, &b_dRSubLeadLeptonLeadJet);
	fChain->SetBranchAddress("dRSubLeadLeptonSubLeadJet", &dRSubLeadLeptonSubLeadJet, &b_dRSubLeadLeptonSubLeadJet);
	fChain->SetBranchAddress("multiLeptonMultiJet_sumPt", &multiLeptonMultiJet_sumPt, &b_multiLeptonMultiJet_sumPt);
	fChain->SetBranchAddress("multiLeptonMultiJet_invMass", &multiLeptonMultiJet_invMass, &b_multiLeptonMultiJet_invMass);
	fChain->SetBranchAddress("diLepton_invMass", &diLepton_invMass, &b_diLepton_invMass);
	fChain->SetBranchAddress("diLepton_pt", &diLepton_pt, &b_diLepton_pt);
	fChain->SetBranchAddress("diJet_invMass", &diJet_invMass, &b_diJet_invMass);
	fChain->SetBranchAddress("diJetLeadingLepton_invMass", &diJetLeadingLepton_invMass, &b_diJetLeadingLepton_invMass);
	fChain->SetBranchAddress("diJetSubLeadingLepton_invMass", &diJetSubLeadingLepton_invMass, &b_diJetSubLeadingLepton_invMass);
	fChain->SetBranchAddress("leadingEle_passHEEPId", &leadingEle_passHEEPId, &b_leadingEle_passHEEPId);
	fChain->SetBranchAddress("leadingEle_passMediumId", &leadingEle_passMediumId, &b_leadingEle_passMediumId);
	fChain->SetBranchAddress("leadingEle_iso", &leadingEle_iso, &b_leadingEle_iso);
	fChain->SetBranchAddress("leadingEle_isMatchedToGen", &leadingEle_isMatchedToGen, &b_leadingEle_isMatchedToGen);
	if (saveHEEPvariables_) {
		fChain->SetBranchAddress("leadingEle_etaSC", &leadingEle_etaSC, &b_leadingEle_etaSC);
		fChain->SetBranchAddress("leadingEle_isEcalDriven", &leadingEle_isEcalDriven, &b_leadingEle_isEcalDriven);
		fChain->SetBranchAddress("leadingEle_dEtaIn", &leadingEle_dEtaIn, &b_leadingEle_dEtaIn);
		fChain->SetBranchAddress("leadingEle_dPhiIn", &leadingEle_dPhiIn, &b_leadingEle_dPhiIn);
		fChain->SetBranchAddress("leadingEle_hOverE", &leadingEle_hOverE, &b_leadingEle_hOverE);
		fChain->SetBranchAddress("leadingEle_full5x5_r9", &leadingEle_full5x5_r9, &b_leadingEle_full5x5_r9);
		fChain->SetBranchAddress("leadingEle_full5x5_sigmaIetaIeta", &leadingEle_full5x5_sigmaIetaIeta, &b_leadingEle_full5x5_sigmaIetaIeta);
		fChain->SetBranchAddress("leadingEle_full5x5_E5x5", &leadingEle_full5x5_E5x5, &b_leadingEle_full5x5_E5x5);
		fChain->SetBranchAddress("leadingEle_full5x5_E1x5", &leadingEle_full5x5_E1x5, &b_leadingEle_full5x5_E1x5);
		fChain->SetBranchAddress("leadingEle_full5x5_E2x5", &leadingEle_full5x5_E2x5, &b_leadingEle_full5x5_E2x5);
		fChain->SetBranchAddress("leadingEle_EmHadDepth1Iso", &leadingEle_EmHadDepth1Iso, &b_leadingEle_EmHadDepth1Iso);
		fChain->SetBranchAddress("leadingEle_ptTracksIso", &leadingEle_ptTracksIso, &b_leadingEle_ptTracksIso);
		fChain->SetBranchAddress("leadingEle_innerLayerLostHits", &leadingEle_innerLayerLostHits, &b_leadingEle_innerLayerLostHits);
		fChain->SetBranchAddress("leadingEle_dxy", &leadingEle_dxy, &b_leadingEle_dxy);
		fChain->SetBranchAddress("leadingEle_dz", &leadingEle_dz, &b_leadingEle_dz);
		fChain->SetBranchAddress("leadingEle_eOverP", &leadingEle_eOverP, &b_leadingEle_eOverP);
		fChain->SetBranchAddress("leadingEle_ecalEnergy", &leadingEle_ecalEnergy, &b_leadingEle_ecalEnergy);
		fChain->SetBranchAddress("leadingEle_hcalOverEcal", &leadingEle_hcalOverEcal, &b_leadingEle_hcalOverEcal);
	}
	fChain->SetBranchAddress("subLeadingEle_passHEEPId", &subLeadingEle_passHEEPId, &b_subLeadingEle_passHEEPId);
	fChain->SetBranchAddress("subLeadingEle_passMediumId", &subLeadingEle_passMediumId, &b_subLeadingEle_passMediumId);
	fChain->SetBranchAddress("subLeadingEle_iso", &subLeadingEle_iso, &b_subLeadingEle_iso);
	fChain->SetBranchAddress("subLeadingEle_isMatchedToGen", &subLeadingEle_isMatchedToGen, &b_subLeadingEle_isMatchedToGen);
	if (saveHEEPvariables_) {
		fChain->SetBranchAddress("subLeadingEle_etaSC", &subLeadingEle_etaSC, &b_subLeadingEle_etaSC);
		fChain->SetBranchAddress("subLeadingEle_isEcalDriven", &subLeadingEle_isEcalDriven, &b_subLeadingEle_isEcalDriven);
		fChain->SetBranchAddress("subLeadingEle_dEtaIn", &subLeadingEle_dEtaIn, &b_subLeadingEle_dEtaIn);
		fChain->SetBranchAddress("subLeadingEle_dPhiIn", &subLeadingEle_dPhiIn, &b_subLeadingEle_dPhiIn);
		fChain->SetBranchAddress("subLeadingEle_hOverE", &subLeadingEle_hOverE, &b_subLeadingEle_hOverE);
		fChain->SetBranchAddress("subLeadingEle_full5x5_r9", &subLeadingEle_full5x5_r9, &b_subLeadingEle_full5x5_r9);
		fChain->SetBranchAddress("subLeadingEle_full5x5_sigmaIetaIeta", &subLeadingEle_full5x5_sigmaIetaIeta, &b_subLeadingEle_full5x5_sigmaIetaIeta);
		fChain->SetBranchAddress("subLeadingEle_full5x5_E5x5", &subLeadingEle_full5x5_E5x5, &b_subLeadingEle_full5x5_E5x5);
		fChain->SetBranchAddress("subLeadingEle_full5x5_E1x5", &subLeadingEle_full5x5_E1x5, &b_subLeadingEle_full5x5_E1x5);
		fChain->SetBranchAddress("subLeadingEle_full5x5_E2x5", &subLeadingEle_full5x5_E2x5, &b_subLeadingEle_full5x5_E2x5);
		fChain->SetBranchAddress("subLeadingEle_EmHadDepth1Iso", &subLeadingEle_EmHadDepth1Iso, &b_subLeadingEle_EmHadDepth1Iso);
		fChain->SetBranchAddress("subLeadingEle_ptTracksIso", &subLeadingEle_ptTracksIso, &b_subLeadingEle_ptTracksIso);
		fChain->SetBranchAddress("subLeadingEle_innerLayerLostHits", &subLeadingEle_innerLayerLostHits, &b_subLeadingEle_innerLayerLostHits);
		fChain->SetBranchAddress("subLeadingEle_dxy", &subLeadingEle_dxy, &b_subLeadingEle_dxy);
		fChain->SetBranchAddress("subLeadingEle_dz", &subLeadingEle_dz, &b_subLeadingEle_dz);
		fChain->SetBranchAddress("subLeadingEle_eOverP", &subLeadingEle_eOverP, &b_subLeadingEle_eOverP);
		fChain->SetBranchAddress("subLeadingEle_ecalEnergy", &subLeadingEle_ecalEnergy, &b_subLeadingEle_ecalEnergy);
		fChain->SetBranchAddress("subLeadingEle_hcalOverEcal", &subLeadingEle_hcalOverEcal, &b_subLeadingEle_hcalOverEcal);
	}
	fChain->SetBranchAddress("leadingMuon_iso", &leadingMuon_iso, &b_leadingMuon_iso);
	fChain->SetBranchAddress("leadingMuon_PFiso", &leadingMuon_PFiso, &b_leadingMuon_PFiso);
	fChain->SetBranchAddress("leadingMuon_isHighPt", &leadingMuon_isHighPt, &b_leadingMuon_isHighPt);
	fChain->SetBranchAddress("leadingMuon_isMatchedToGen", &leadingMuon_isMatchedToGen, &b_leadingMuon_isMatchedToGen);
	fChain->SetBranchAddress("leadingMuon_dz", &leadingMuon_dz, &b_leadingMuon_dz);
	fChain->SetBranchAddress("leadingMuon_dxy", &leadingMuon_dxy, &b_leadingMuon_dxy);
	fChain->SetBranchAddress("leadingMuon_RochCor", &leadingMuon_RochCor, &b_leadingMuon_RochCor);
	fChain->SetBranchAddress("subLeadingMuon_iso", &subLeadingMuon_iso, &b_subLeadingMuon_iso);
	fChain->SetBranchAddress("subLeadingMuon_PFiso", &subLeadingMuon_PFiso, &b_subLeadingMuon_PFiso);
	fChain->SetBranchAddress("subLeadingMuon_isHighPt", &subLeadingMuon_isHighPt, &b_subLeadingMuon_isHighPt);
	fChain->SetBranchAddress("subLeadingMuon_isMatchedToGen", &subLeadingMuon_isMatchedToGen, &b_subLeadingMuon_isMatchedToGen);
	fChain->SetBranchAddress("subLeadingMuon_dz", &subLeadingMuon_dz, &b_subLeadingMuon_dz);
	fChain->SetBranchAddress("subLeadingMuon_dxy", &subLeadingMuon_dxy, &b_subLeadingMuon_dxy);
	fChain->SetBranchAddress("subLeadingMuon_RochCor", &subLeadingMuon_RochCor, &b_subLeadingMuon_RochCor);

	cout << "Branches are properly initialized." << endl;
}
// ******************************************************************************************


// **************** 
Int_t makePlots_multiLeptonMultiJet::GetEntry(Long64_t entry){
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
// ******************************************************************************************


// **************** 
Long64_t makePlots_multiLeptonMultiJet::LoadTree(Long64_t entry){
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
	}
	return centry;
}
// ******************************************************************************************



// **************** 
TH1D* makePlots_multiLeptonMultiJet::newTH1D(string name, string title, string xTitle, int nBins, double xLow, double xUp){
	TH1D* hist = new TH1D(name.c_str(), title.c_str(), nBins, xLow, xUp);
	hist->GetXaxis()->SetTitle(xTitle.c_str());
	hist->GetYaxis()->SetTitle("# Events");
	hist->GetYaxis()->SetTitleOffset(1.9);
	hist->SetOption("HIST");           
	listOfHistograms.push_back(hist);
	return hist;
}
// ******************************************************************************************



// **************** 
TH1D* makePlots_multiLeptonMultiJet::newTH1DvarBinSize(string name, string title, string xTitle, int nBins, const double *xbins){
	// Int_t  nBins = sizeof(xbins) / sizeof(Float_t) - 1;
	// cout << "nBins = " << nBins << endl;
	TH1D* hist = new TH1D(name.c_str(), title.c_str(), nBins, xbins);
	hist->GetXaxis()->SetTitle(xTitle.c_str());
	hist->GetYaxis()->SetTitle("# Events");
	hist->GetYaxis()->SetTitleOffset(1.9);
	hist->SetOption("HIST");           
	listOfHistograms.push_back(hist);
	return hist;
}
// ******************************************************************************************



// **************** 
TH2D* makePlots_multiLeptonMultiJet::newTH2D(string name, string title, string xTitle, string yTitle, int nBinsX, double xLow, double xUp, int nBinsY, double yLow, double yUp){
	TH2D* hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xLow, xUp, nBinsY, yLow, yUp);
	hist->GetXaxis()->SetTitle(xTitle.c_str());
	hist->GetYaxis()->SetTitle(yTitle.c_str());
	hist->GetYaxis()->SetTitleOffset(1.9); 
	listOfHistograms.push_back(hist);
	return hist;
}
// ******************************************************************************************


// **************** 
TProfile* makePlots_multiLeptonMultiJet::newTProfile(string name, string title, string xTitle, string yTitle, int nBinsX, double xLow, double xUp, double yLow, double yUp){
	TProfile* hist = new TProfile(name.c_str(), title.c_str(), nBinsX, xLow, xUp, yLow, yUp);
	hist->GetXaxis()->SetTitle(xTitle.c_str());
	hist->GetYaxis()->SetTitle(yTitle.c_str());
	hist->GetYaxis()->SetTitleOffset(1.9);
	listOfHistograms.push_back(hist);
	return hist;
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::SetHistos(){

	if (signalEE || signalMuMu) {
		numberRegions = 3;
		for (int i(0); i < numberRegions; i++) regionName[i] = signalRegionNames[i];		
	}
	if (eMuSideband) {
		numberRegions = 1;
		for (int i(0); i < numberRegions; i++) regionName[i] = eMuSidebandName[i];
	}
	if (TnPEE || TnPMuMu) {
		numberRegions = 1;
		for (int i(0); i < numberRegions; i++) regionName[i] = TnPCRName[i];
	}

	if (signalEE) triggerName = "_signalEE";
	if (signalMuMu) triggerName = "_signalMuMu";
	if (eMuSideband) triggerName = "_eMuSideband";
	if (TnPEE) triggerName = "_TnPEE";
	if (TnPMuMu) triggerName = "_TnPMuMu";



	rho_histo = newTH1D("rho_histo"+triggerName, "rho_histo"+triggerName, "rho", 100, 0., 100.);
	nvtx_histo = newTH1D("nvtx_histo"+triggerName, "nvtx_histo"+triggerName, "nvtx", 100, 0., 100.);
	vtx_x_histo = newTH1D("vtx_x_histo"+triggerName, "vtx_x_histo"+triggerName, "vtx_x", 100, -0.3, 0.3);
	vtx_y_histo = newTH1D("vtx_y_histo"+triggerName, "vtx_y_histo"+triggerName, "vtx_y", 100, -0.3, 0.3);
	vtx_z_histo = newTH1D("vtx_z_histo"+triggerName, "vtx_z_histo"+triggerName, "vtx_z", 100, -20., 20.);

	// muon_dxy_histo = newTH1D("muon_dxy_histo", "muon_dxy_histo", "|dxy| muons", 100, 0., 0.5);

	// nLeadingLeptons_histo = newTH1D("nLeadingLeptons_histo", "nLeadingLeptons_histo", "nLeadingLepton", 10, 0, 10); 
	// nSubLeadingLeptons_histo = newTH1D("nSubLeadingLeptons_histo", "nSubLeadingLeptons_histo", "nSubLeadingLepton", 10, 0, 10); 



	const double mlmjBins[] = {150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500, 1650, 1800, 1950, 2100, 2250, 2400, 2550, 2700, 2850, 3000, 3150, 3300, 3450, 3600, 3750, 3900, 4150, 7000};
	// Int_t  nBins = sizeof(mlmjBins) / sizeof(Float_t) - 1;
	Int_t  nBins = 27;
	// cout << "nBins = " << nBins << endl;

	for (unsigned short l(0); l < numberRegions; l++) {
		for (unsigned short i(0); i < 8; i++) {
			pt_histo[i][l] = newTH1D("pt_"+objName[i]+regionName[l]+"_histo", "pt_"+objName[i]+regionName[l]+"_histo", "p_{T} [GeV] "+objName[i]+regionName[l], 200, 0., 1000.);
			eta_histo[i][l] = newTH1D("eta_"+objName[i]+regionName[l]+"_histo", "eta_"+objName[i]+regionName[l]+"_histo", "#eta "+objName[i]+regionName[l], 100, -2.5, 2.5);
			phi_histo[i][l] = newTH1D("phi_"+objName[i]+regionName[l]+"_histo", "phi_"+objName[i]+regionName[l]+"_histo", "#phi "+objName[i]+regionName[l], 100, -3.5, 3.5);		
		} 

		for (unsigned short n(0); n < nEtaMassRegions; n++) {
			// mass_dldj_histo[l][n] = newTH1D("mass_multiLeptonMultiJets"+regionName[l]+etaMassName[n]+"_histo", "mass_multiLeptonMultiJets"+regionName[l]+etaMassName[n]+"_histo", "m_{lljj}"+regionName[l], 300, 0., 6000.);
			mass_dldj_histo[l][n] = newTH1DvarBinSize("mass_multiLeptonMultiJets"+regionName[l]+etaMassName[n]+"_histo", "mass_multiLeptonMultiJets"+regionName[l]+etaMassName[n]+"_histo", "m_{lljj}"+regionName[l], nBins, mlmjBins);
			mass_dl_histo[l][n] = newTH1D("mass_dileptons"+regionName[l]+etaMassName[n]+"_histo", "mass_dileptons"+regionName[l]+etaMassName[n]+"_histo", "m_{ll}"+regionName[l], 100, 0., 2000.);
			mass_dj_histo[l][n] = newTH1D("mass_dijets"+regionName[l]+etaMassName[n]+"_histo", "mass_dijets"+regionName[l]+etaMassName[n]+"_histo", "m_{jj}"+regionName[l], 100, 0., 2000.);		
			mass_djLl_histo[l][n] = newTH1D("mass_dijetsLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "mass_dijetsLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "m_{jjl_{L}}"+regionName[l], 200, 0., 4000.);	
			mass_djSLl_histo[l][n] = newTH1D("mass_dijetsSubLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "mass_dijetsSubLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "m_{jjl_{SL}}"+regionName[l], 200, 0., 4000.);
		} 

		pt_dl_histo[l] = newTH1D("pt_dileptons"+regionName[l]+"_histo", "pt_dileptons_"+regionName[l]+"_histo", "diLepton p_{T} [GeV]", 100, 0., 1000.);
	}

	if ( !eMuSideband ) {
		for (unsigned short z(0); z < nZregions; z++) {
			for (unsigned short n(0); n < nEtaMassRegions; n++) {
				Zmass_histo[z][n] = newTH1D("Z"+Zname[z]+triggerName+"_mass"+etaMassName[n]+"_histo", "Z"+Zname[z]+triggerName+"_mass"+etaMassName[n]+"_histo", "m(Z) [GeV/c^{2}]", 200, 0, 200);
			}
			Zpt_histo[z] = newTH1D("Z"+Zname[z]+triggerName+"_pt_histo", "Z"+Zname[z]+triggerName+"_pt_histo", "pt(Z) [GeV/c^{2}]", 200, 50, 250); 
		}
	}

	for (unsigned short l(0); l < 1; l++) {
		for (unsigned short i(8); i < 16; i++) {
			pt_histo[i][l] = newTH1D("pt_"+objName[i]+triggerName+"_histo", "pt_"+objName[i]+triggerName+"_histo", "p_{T} [GeV] "+objName[i], 200, 0., 1000.);
			eta_histo[i][l] = newTH1D("eta_"+objName[i]+triggerName+"_histo", "eta_"+objName[i]+triggerName+"_histo", "#eta "+objName[i], 100, -2.5, 2.5);
			phi_histo[i][l] = newTH1D("phi_"+objName[i]+triggerName+"_histo", "phi_"+objName[i]+triggerName+"_histo", "#phi "+objName[i], 100, -3.5, 3.5);		
		}
	}

	if (saveHEEPvariables_) {
		for (unsigned short j(0); j < nEle; j++) {
			for (unsigned short m(0); m < nEtaRegions; m++) {
				etaSC_histo[j][m] = newTH1D("etaSC_"+eleName[j]+etaName[m]+"_histo", "etaSC_"+eleName[j]+etaName[m]+"_histo", "#eta_{SC} "+eleName[j], 100, -2.5, 2.5);
				isEcalDriven_histo[j][m] = newTH1D("isEcalDriven_"+eleName[j]+etaName[m]+"_histo", "isEcalDriven_"+eleName[j]+etaName[m]+"_histo", "isEcalDriven "+eleName[j], 2, 0., 2.);
				dEtaIn_histo[j][m] = newTH1D("dEtaIn_"+eleName[j]+etaName[m]+"_histo", "dEtaIn_"+eleName[j]+etaName[m]+"_histo", "dEtaIn "+eleName[j], 100, -0.2, 0.2);
				dPhiIn_histo[j][m] = newTH1D("dPhiIn_"+eleName[j]+etaName[m]+"_histo", "dPhiIn_"+eleName[j]+etaName[m]+"_histo", "dPhiIn "+eleName[j], 100, -0.2, 0.2);
				hOverE_histo[j][m] = newTH1D("hOverE_"+eleName[j]+etaName[m]+"_histo", "hOverE_"+eleName[j]+etaName[m]+"_histo", "H/E "+eleName[j], 40, 0., 2.);
				full5x5_r9_histo[j][m] = newTH1D("r9_"+eleName[j]+etaName[m]+"_histo", "r9_"+eleName[j]+etaName[m]+"_histo", "r9 "+eleName[j], 100, 0., 1.05);
				full5x5_sigmaIetaIeta_histo[j][m] = newTH1D("sigmaIetaIeta_"+eleName[j]+etaName[m]+"_histo", "sigmaIetaIeta_"+eleName[j]+etaName[m]+"_histo", "sigmaIetaIeta "+eleName[j], 100, 0., 0.05);
				full5x5_E2x5_Over_E5x5_histo[j][m] = newTH1D("e2x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e2x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e2x5/e5x5 "+eleName[j], 60, 0.5, 1.1);
				full5x5_E1x5_Over_E5x5_histo[j][m] = newTH1D("e1x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e1x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e1x5/e5x5 "+eleName[j], 100, 0., 1.1);
				EmHadDepth1Iso_histo[j][m] = newTH1D("EmHadDepth1Iso_"+eleName[j]+etaName[m]+"_histo", "EmHadDepth1Iso_"+eleName[j]+etaName[m]+"_histo", "EmHadDepth1Iso "+eleName[j], 200, 0., 200.);
				// ptTracksIso_histo[j][m];
				innerLayerLostHits_histo[j][m] = newTH1D("missingHits_"+eleName[j]+etaName[m]+"_histo", "missingHits_"+eleName[j]+etaName[m]+"_histo", "missingHits "+eleName[j], 10, 0., 10.);
				dxy_histo[j][m] = newTH1D("dxy_"+eleName[j]+etaName[m]+"_histo", "dxy_"+eleName[j]+etaName[m]+"_histo", "|dxy| "+eleName[j], 100, 0., 0.5);
				eOverP[j][m] = newTH1D("eOverP_"+eleName[j]+etaName[m]+"_histo", "eOverP_"+eleName[j]+etaName[m]+"_histo", "E/p "+eleName[j], 100, 0., 100);
			}
		}
	}

}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getElectrons(vector<eleStruct>& electrons){//, vector<eleStruct>& goodElectrons){
	unsigned short nTotElectrons(ele_e->size());
	for (unsigned short i(0); i < nTotElectrons; i++){
		eleStruct ele(ele_pt->at(i), 
						ele_eta->at(i), 						
						ele_phi->at(i), 
						ele_e->at(i), 
						ele_passHEEPId->at(i),
						0,
						ele_passMediumId->at(i), 
						ele_iso->at(i), 
						ele_dz->at(i),
						ele_d0->at(i),						
						ele_isMatchedToGen->at(i),						
						ele_charge->at(i),
						0,
						0,
						0,
						0,
						0,
						0,	
						0,
						0,
						0,
						0, 
						0,
						0,
						0, 
						0,
						0,
						0,
						0
					); 

		if (saveHEEPvariables_) {
			eleStruct ele(ele_pt->at(i), 
						ele_eta->at(i), 						
						ele_phi->at(i), 
						ele_e->at(i), 
						ele_passHEEPId->at(i),
						0,
						ele_passMediumId->at(i), 
						ele_iso->at(i), 
						ele_dz->at(i),
						ele_d0->at(i),						
						ele_isMatchedToGen->at(i),						
						ele_charge->at(i),
						ele_etaSC->at(i),
						ele_isEcalDriven->at(i),
						ele_dEtaIn->at(i),
						ele_dPhiIn->at(i),
						ele_hOverE->at(i),
						ele_full5x5_r9->at(i),	
						ele_full5x5_sigmaIetaIeta->at(i),
						ele_full5x5_E5x5->at(i),
						ele_full5x5_E1x5->at(i),
						ele_full5x5_E2x5->at(i), 
						ele_EmHadDepth1Iso->at(i),
						ele_ptTracksIso->at(i),
						ele_innerLayerLostHits->at(i), 
						ele_dxy->at(i),
						ele_eOverP->at(i),
						ele_ecalEnergy->at(i),
						ele_hcalOverEcal->at(i)
					); 
		}
		electrons.push_back(ele);
		// if (fabs(ele.v.Eta()) < 2.4 && ele.v.Pt() > 53 && ele.passHEEPId) goodElectrons.push_back(ele);
	} //End of loop over all the electrons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getLeadingElectrons(vector<eleStruct>& leadingElectrons){
	unsigned short nTotLeadingElectrons(leadingEle_iso->size());
	for (unsigned short i(0); i < nTotLeadingElectrons; i++){
		eleStruct ele(leadingLepton_pt->at(i),
						leadingLepton_eta->at(i),
						leadingLepton_phi->at(i),
						leadingLepton_e->at(i), 
						leadingEle_passHEEPId->at(i),
						0,
						leadingEle_passMediumId->at(i),
						leadingEle_iso->at(i),
						0,
						0,
						leadingEle_isMatchedToGen->at(i),
						leadingLepton_charge->at(i), 
						0,
						0,
						0,
						0,
						0,
						0,
						0,
						0,
						0,
						0,		
						0,
						0,
						0,
						0,
						0,
						0,
						0
					);

		if (saveHEEPvariables_) {
			eleStruct ele(leadingLepton_pt->at(i),
						leadingLepton_eta->at(i),
						leadingLepton_phi->at(i),
						leadingLepton_e->at(i), 
						leadingEle_passHEEPId->at(i),
						0,
						leadingEle_passMediumId->at(i),
						leadingEle_iso->at(i),
						leadingEle_dz->at(i),
						0,
						leadingEle_isMatchedToGen->at(i),
						leadingLepton_charge->at(i), 
						leadingEle_etaSC->at(i),
						leadingEle_isEcalDriven->at(i),
						leadingEle_dEtaIn->at(i),
						leadingEle_dPhiIn->at(i),
						leadingEle_hOverE->at(i),
						leadingEle_full5x5_r9->at(i),
						leadingEle_full5x5_sigmaIetaIeta->at(i),
						leadingEle_full5x5_E5x5->at(i),
						leadingEle_full5x5_E1x5->at(i),
						leadingEle_full5x5_E2x5->at(i),		
						leadingEle_EmHadDepth1Iso->at(i),
						leadingEle_ptTracksIso->at(i),
						leadingEle_innerLayerLostHits->at(i),
						leadingEle_dxy->at(i),
						leadingEle_eOverP->at(i),
						leadingEle_ecalEnergy->at(i),
						leadingEle_hcalOverEcal->at(i)
					);			
		}
		leadingElectrons.push_back(ele);
	} //End of loop over all the leading electrons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getSubLeadingElectrons(vector<eleStruct>& subLeadingElectrons){
	unsigned short nTotSubLeadingElectrons(subLeadingEle_iso->size());
	for (unsigned short i(0); i < nTotSubLeadingElectrons; i++){
		eleStruct ele(subLeadingLepton_pt->at(i),
						subLeadingLepton_eta->at(i),
						subLeadingLepton_phi->at(i),
						subLeadingLepton_e->at(i), 
						subLeadingEle_passHEEPId->at(i),
						0,
						subLeadingEle_passMediumId->at(i),
						subLeadingEle_iso->at(i),
						0,
						0,
						subLeadingEle_isMatchedToGen->at(i),
						subLeadingLepton_charge->at(i), 
						0,
						0,
						0,
						0,
						0,
						0,
						0,
						0,
						0,
						0,		
						0,
						0,
						0,
						0,
						0,
						0,
						0				
					);	

		if (saveHEEPvariables_) {
			eleStruct ele(subLeadingLepton_pt->at(i),
						subLeadingLepton_eta->at(i),
						subLeadingLepton_phi->at(i),
						subLeadingLepton_e->at(i), 
						subLeadingEle_passHEEPId->at(i),
						0,
						subLeadingEle_passMediumId->at(i),
						subLeadingEle_iso->at(i),
						subLeadingEle_dz->at(i),
						0,
						subLeadingEle_isMatchedToGen->at(i),
						subLeadingLepton_charge->at(i), 
						subLeadingEle_etaSC->at(i),
						subLeadingEle_isEcalDriven->at(i),
						subLeadingEle_dEtaIn->at(i),
						subLeadingEle_dPhiIn->at(i),
						subLeadingEle_hOverE->at(i),
						subLeadingEle_full5x5_r9->at(i),
						subLeadingEle_full5x5_sigmaIetaIeta->at(i),
						subLeadingEle_full5x5_E5x5->at(i),
						subLeadingEle_full5x5_E1x5->at(i),
						subLeadingEle_full5x5_E2x5->at(i),		
						subLeadingEle_EmHadDepth1Iso->at(i),
						subLeadingEle_ptTracksIso->at(i),
						subLeadingEle_innerLayerLostHits->at(i),
						subLeadingEle_dxy->at(i),
						subLeadingEle_eOverP->at(i),
						subLeadingEle_ecalEnergy->at(i),
						subLeadingEle_hcalOverEcal->at(i)				
					);
		} 		
		subLeadingElectrons.push_back(ele);
	} //End of loop over all the subLeading electrons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getMuons(vector<muonStruct>& muons, vector<muonStruct>& muonsRochCorr){ //, vector<muonStruct>& goodMuons, vector<muonStruct>& goodMuonsRochCorr){
	unsigned short nTotMuons(mu_e->size());
	for (unsigned short i(0); i < nTotMuons; i++) {
		muonStruct mu(mu_pt->at(i), 
				mu_eta->at(i), 
				mu_phi->at(i), 
				mu_e->at(i), 
				mu_iso->at(i),	
				mu_isHighPt->at(i),
				mu_isMatchedToGen->at(i),
				mu_charge->at(i),
				mu_dz->at(i),
				mu_dxy->at(i),
				mu_RochCor->at(i)
			);					
		muons.push_back(mu);
		// if (fabs(mu.v.Eta()) < 2.4 && mu.v.Pt() > 53 && mu.isHighPt) goodMuons.push_back(mu);

		mu.v.SetPtEtaPhiE(mu.v.Pt()*mu.rochCor, mu.v.Eta(), mu.v.Phi(), mu.v.E());
		muonsRochCorr.push_back(mu);
		// if (fabs(mu.v.Eta()) < 2.4 && mu.v.Pt() > 53 && mu.isHighPt) goodMuonsRochCorr.push_back(mu);
	}//End of loop over all the muons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getLeadingMuons(vector<muonStruct>& leadingMuons){//, vector<muonStruct>& leadingMuonsRochCorr){
	unsigned short nTotLeadingMuons(leadingMuon_isHighPt->size());
	for (unsigned short i(0); i < nTotLeadingMuons; i++) {
		muonStruct mu(leadingLepton_pt->at(i),
					leadingLepton_eta->at(i),
					leadingLepton_phi->at(i),
					leadingLepton_e->at(i),
					leadingMuon_iso->at(i),
					leadingMuon_isHighPt->at(i),
					leadingMuon_isMatchedToGen->at(i),
					leadingLepton_charge->at(i), 	
					leadingMuon_dz->at(i),
					leadingMuon_dxy->at(i),
					leadingMuon_RochCor->at(i)
				);					
		leadingMuons.push_back(mu);
	}//End of loop over all the leading muons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getSubLeadingMuons(vector<muonStruct>& subLeadingMuons){//, vector<muonStruct>& subLeadingMuonsRochCorr){
	unsigned short nTotSubLeadingMuons(subLeadingMuon_isHighPt->size());
	for (unsigned short i(0); i < nTotSubLeadingMuons; i++) {
		muonStruct mu(subLeadingLepton_pt->at(i),
					subLeadingLepton_eta->at(i),
					subLeadingLepton_phi->at(i),
					subLeadingLepton_e->at(i),
					subLeadingMuon_iso->at(i),
					subLeadingMuon_isHighPt->at(i),
					subLeadingMuon_isMatchedToGen->at(i),
					subLeadingLepton_charge->at(i), 	
					subLeadingMuon_dz->at(i),
					subLeadingMuon_dxy->at(i),
					subLeadingMuon_RochCor->at(i)
				);					
		subLeadingMuons.push_back(mu);
	}//End of loop over all the subleading muons
}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::getJets(vector<jetStruct>& jets){//, vector<jetStruct>& goodJets){
	unsigned short nTotJets(jet_e->size());
	for (unsigned short i(0); i < nTotJets; i++){
		jetStruct jet(jet_pt->at(i),
				jet_eta->at(i),
				jet_phi->at(i),
				jet_e->at(i),
				jet_isMatchedToGen->at(i),
				jet_isThight->at(i)
			);
		jets.push_back(jet);
		// if (fabs(jet.v.Eta()) < 2.4 && jet.v.Pt() > 40 && jet.isTight) goodJets.push_back(jet);
	} //End of loop over all the jets
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getLeadingJets(vector<jetStruct>& leadingJets){
	unsigned short nTotLeadingJets(leadingJet_pt->size());
	for (unsigned short i(0); i < nTotLeadingJets; i++){
		jetStruct jet(leadingJet_pt->at(i),
				leadingJet_eta->at(i),
				leadingJet_phi->at(i),
				leadingJet_e->at(i),
				leadingJet_isMatchedToGen->at(i),
				leadingJet_isThight->at(i)
			);
		leadingJets.push_back(jet);
	} //End of loop over all the leading jets
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getSubLeadingJets(vector<jetStruct>& subLeadingJets){
	unsigned short nTotSubLeadingJets(subLeadingJet_pt->size());
	for (unsigned short i(0); i < nTotSubLeadingJets; i++){
		jetStruct jet(subLeadingJet_pt->at(i),
				subLeadingJet_eta->at(i),
				subLeadingJet_phi->at(i),
				subLeadingJet_e->at(i),
				subLeadingJet_isMatchedToGen->at(i),
				subLeadingJet_isThight->at(i)
			);
		subLeadingJets.push_back(jet);
	} //End of loop over all the leading jets
}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::getMultiLeptonMultiJets(vector<multiLeptonMultiJetStruct>& multiLeptonMultiJets){
	unsigned short nTotmultiLeptonMultiJets(isEEJJ->size());
	for (unsigned short i(0); i < nTotmultiLeptonMultiJets; i++){
		multiLeptonMultiJetStruct mlmj(isEEJJ->at(i),
							isEETT->at(i),
							isMMJJ->at(i),
							isMMTT->at(i),
							isEMJJ->at(i),
							isSignalRegion->at(i),
							isLowMllCR->at(i),
							isLowMlljjCR->at(i),
							isBB->at(i),
							isEE->at(i),
							isEB->at(i),
							passPreselections->at(i),
							dRLeadLeptonLeadJet->at(i),
							dRLeadLeptonSubLeadJet->at(i),
							dRSubLeadLeptonLeadJet->at(i),
							dRSubLeadLeptonSubLeadJet->at(i),
							multiLeptonMultiJet_sumPt->at(i),
							multiLeptonMultiJet_invMass->at(i),
							diLepton_invMass->at(i),
							diLepton_pt->at(i), 
							diJet_invMass->at(i), 
							diJetLeadingLepton_invMass->at(i),
							diJetSubLeadingLepton_invMass->at(i)
						);
		multiLeptonMultiJets.push_back(mlmj);
	} //End of loop over all the multiLeptonMultiJets
}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::doEleDistributionsPlots(vector<eleStruct>& electrons, const int eleIdx, const int variableIdx, const int regionIdx) { 
	pt_histo[variableIdx][regionIdx]->Fill(electrons[eleIdx].v.Pt(),w);
	eta_histo[variableIdx][regionIdx]->Fill(electrons[eleIdx].v.Eta(),w);
	phi_histo[variableIdx][regionIdx]->Fill(electrons[eleIdx].v.Phi(),w); 
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::doMuonDistributionsPlots(vector<muonStruct>& muons, const int muIdx, const int variableIdx, const int regionIdx) { 
	pt_histo[variableIdx][regionIdx]->Fill(muons[muIdx].v.Pt(),w);
	eta_histo[variableIdx][regionIdx]->Fill(muons[muIdx].v.Eta(),w);
	phi_histo[variableIdx][regionIdx]->Fill(muons[muIdx].v.Phi(),w); 
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::doJetsDistributionsPlots(vector<jetStruct>& jets, const int jetIdx, const int variableIdx, const int regionIdx) { 
	pt_histo[variableIdx][regionIdx]->Fill(jets[jetIdx].v.Pt(),w);
	eta_histo[variableIdx][regionIdx]->Fill(jets[jetIdx].v.Eta(),w);
	phi_histo[variableIdx][regionIdx]->Fill(jets[jetIdx].v.Phi(),w); 
}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::doElePlots(vector<eleStruct>& electrons, const int eleIdx, const int variableIdx, const int etaIdx){
	etaSC_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].etaSC,w);
	isEcalDriven_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].isEcalDriven,w);
	dEtaIn_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].dEtaIn,w);
	dPhiIn_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].dPhiIn,w);
	hOverE_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].HoE,w);
	full5x5_r9_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].r9,w);
	full5x5_sigmaIetaIeta_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].sigmaIetaIeta,w);
	full5x5_E2x5_Over_E5x5_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].e2x5/electrons[eleIdx].e5x5,w);
	full5x5_E1x5_Over_E5x5_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].e1x5/electrons[eleIdx].e5x5,w);
	EmHadDepth1Iso_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].EmHadDepth1Iso,w);
	innerLayerLostHits_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].missingHits,w);
	dxy_histo[variableIdx][etaIdx]->Fill(electrons[eleIdx].dxy,w);
	eOverP[variableIdx][etaIdx]->Fill(electrons[eleIdx].eOverP,w);
}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::doMassPlots(vector<multiLeptonMultiJetStruct>& multiLeptonMultiJets, const int mlmjIdx, const int regionIdx, const int etaMassIdx){
	mass_dldj_histo[regionIdx][etaMassIdx]->Fill(multiLeptonMultiJets[mlmjIdx].multiLeptonMultiJet_invMass,w);
	mass_dl_histo[regionIdx][etaMassIdx]->Fill(multiLeptonMultiJets[mlmjIdx].diLepton_invMass,w);
	mass_dj_histo[regionIdx][etaMassIdx]->Fill(multiLeptonMultiJets[mlmjIdx].diJet_invMass,w);
	mass_djLl_histo[regionIdx][etaMassIdx]->Fill(multiLeptonMultiJets[mlmjIdx].diJetLeadingLepton_invMass,w);
	mass_djSLl_histo[regionIdx][etaMassIdx]->Fill(multiLeptonMultiJets[mlmjIdx].diJetSubLeadingLepton_invMass,w);

	if (etaMassIdx == 0) pt_dl_histo[regionIdx]->Fill(multiLeptonMultiJets[mlmjIdx].diLepton_pt,w);

}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::doZPlots(vector<multiLeptonMultiJetStruct>& multiLeptonMultiJets, const int mlmjIdx, const int variableIdx) { 	
	float ptZ = multiLeptonMultiJets[mlmjIdx].diLepton_pt;
	float mZ = multiLeptonMultiJets[mlmjIdx].diLepton_invMass;	
	bool isInEBEB = inEBEB(leadingLepton_eta->at(mlmjIdx), subLeadingLepton_eta->at(mlmjIdx));  
	bool isInEEEE = inEEEE(leadingLepton_eta->at(mlmjIdx), subLeadingLepton_eta->at(mlmjIdx));
	bool isInEBEE = inEBEE(leadingLepton_eta->at(mlmjIdx), subLeadingLepton_eta->at(mlmjIdx)); 			

	Zpt_histo[variableIdx]->Fill(ptZ,w);
	Zmass_histo[variableIdx][0]->Fill(mZ,w);
	if (isInEBEB) Zmass_histo[variableIdx][1]->Fill(mZ,w); 
	if (isInEEEE) Zmass_histo[variableIdx][2]->Fill(mZ,w);  
	if (isInEBEE) Zmass_histo[variableIdx][3]->Fill(mZ,w);  
	if (isInEEEE || isInEBEE) Zmass_histo[variableIdx][4]->Fill(mZ,w);  
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::doZPlots(TLorentzVector Llepton, TLorentzVector SLlepton, const int variableIdx) { 	
	float ptZ = (Llepton+SLlepton).Pt();
	float mZ = (Llepton+SLlepton).M();
	bool isInEBEB = inEBEB(Llepton.Eta(), SLlepton.Eta()); 
	bool isInEEEE = inEEEE(Llepton.Eta(), SLlepton.Eta());
	bool isInEBEE = inEBEE(Llepton.Eta(), SLlepton.Eta()); 

	Zpt_histo[variableIdx]->Fill(ptZ,w);
	Zmass_histo[variableIdx][0]->Fill(mZ,w);
	if (isInEBEB) Zmass_histo[variableIdx][1]->Fill(mZ,w); 
	if (isInEEEE) Zmass_histo[variableIdx][2]->Fill(mZ,w);  
	if (isInEBEE) Zmass_histo[variableIdx][3]->Fill(mZ,w);  
	if (isInEEEE || isInEBEE) Zmass_histo[variableIdx][4]->Fill(mZ,w);  
}
// ******************************************************************************************




// **************** 
void makePlots_multiLeptonMultiJet::Loop(){
	//--- Initialize the tree branches ---
	Init();
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	cout << "nentries = " << nentries << endl;


	// for (Long64_t i=809; i<813; i+=1) {
	for (Long64_t i=0; i<nentries; i+=1) {
		Long64_t ientry = LoadTree(i);
		if (ientry < 0) break;
		// cout << "evento " << ientry << endl;
		nEvents++;

		if(fChain->GetEntry(i) == 0){
			std::cerr << "Failed to read Tree entry " << i << "!\n";
			continue;
		}


	// --Collections
		vector<eleStruct> electrons, goodElectrons, leadingElectrons, subLeadingElectrons; 
		vector<muonStruct> muons, muonsRochCorr, goodMuons, goodMuonsRochCorr, leadingMuons, subLeadingMuons; 
		vector<jetStruct> jets, goodJets, leadingJets, subLeadingJets; 
		vector<multiLeptonMultiJetStruct> multiLeptonMultiJets;

		// cout << "get electrons" << endl;
		getElectrons(electrons);//, goodElectrons);

		// cout << "get muons" << endl;
		getMuons(muons, muonsRochCorr);//, goodMuons, goodMuonsRochCorr);

		// cout << "get jets" << endl;
		getJets(jets);//, goodJets);

		// cout << "get mlmj" << endl;
		getMultiLeptonMultiJets(multiLeptonMultiJets);

		// cout << "get leadingEle" << endl;
		getLeadingElectrons(leadingElectrons);
		// cout << "get subLeadingEle" << endl;
		getSubLeadingElectrons(subLeadingElectrons);

		// cout << "get leadingMuon" << endl;
		getLeadingMuons(leadingMuons);
		// cout << "get subLeadingMuon" << endl;
		getSubLeadingMuons(subLeadingMuons);

		// cout << "get leadingJet" << endl;
		getLeadingJets(leadingJets);
		// cout << "get subLeadingJet" << endl;
		getSubLeadingJets(subLeadingJets);


		unsigned short nElectrons = electrons.size();    
		// unsigned short nGoodElectrons = goodElectrons.size();
		unsigned short nMuons = muons.size();
		// unsigned short nGoodMuons = goodMuons.size();
		unsigned short nMuonsRochCor = muonsRochCorr.size();
		// unsigned short nGoodMuonsRochCor = goodMuonsRochCorr.size();		
		unsigned short nJets = jets.size();
		// unsigned short nGoodJets = goodJets.size();
		unsigned short nMultiLeptonMultiJets = multiLeptonMultiJets.size();
		unsigned short nLeadingEle = leadingElectrons.size();
		unsigned short nSubLeadingEle = subLeadingElectrons.size();
		unsigned short nLeadingMuons = leadingMuons.size();
		unsigned short nSubLeadingMuons = subLeadingMuons.size();
		unsigned short nLeadingJets = leadingJets.size();
		unsigned short nSubLeadingJets = subLeadingJets.size();

		if ( !(nMuonsRochCor == nMuons) ) cout << "nMuonsRochCor != nMuons" << endl;
		if ( !(nLeadingEle == nMultiLeptonMultiJets) ) cout << "nLeadingEle != nMultiLeptonMultiJets" << endl;
		if ( !(nSubLeadingEle == nMultiLeptonMultiJets) ) cout << "nSubLeadingEle != nMultiLeptonMultiJets" << endl;
		if ( !(nLeadingMuons == nMultiLeptonMultiJets) ) cout << "nLeadingMuons != nMultiLeptonMultiJets" << endl;
		if ( !(nSubLeadingMuons == nMultiLeptonMultiJets) ) cout << "nSubLeadingMuons != nMultiLeptonMultiJets" << endl;
		if ( !(nLeadingJets == nMultiLeptonMultiJets) ) cout << "nLeadingJets != nMultiLeptonMultiJets" << endl;
		if ( !(nSubLeadingJets == nMultiLeptonMultiJets) ) cout << "nSubLeadingJets != nMultiLeptonMultiJets" << endl;


		w = weight / nToDivideWeight;
		// cout << "weight = " << weight << ", nToDivideWeight = " << nToDivideWeight << " -> w = " << w << endl;
		// w = weight * nentries / nTotEvt;


	// --vtx
		rho_histo->Fill(rho,w);
		nvtx_histo->Fill(nvtx,w);

		for (unsigned short i(0); i < nvtx; i++){
			vtx_x_histo->Fill(vtx_x->at(i),w);
			vtx_y_histo->Fill(vtx_y->at(i),w);
			vtx_z_histo->Fill(vtx_z->at(i),w);
		}


	// --sort SingleObjects
		for (unsigned short j(0); j < nElectrons-1; j++){
			if ( electrons[j+1].v.Pt() > electrons[j].v.Pt() ) swap(electrons[j+1], electrons[j]);
		}

		for (unsigned short j(0); j < nMuons-1; j++){
			if ( muons[j+1].v.Pt() > muons[j].v.Pt() ) swap(muons[j+1], muons[j]);
		}

		for (unsigned short l(0); l < nJets-1; l++){
			if ( jets[l+1].v.Pt() > jets[l].v.Pt() ) swap(jets[l+1], jets[l]);
		}



	// --definizione tagli
		float ptMinLepton = 53;
		if (TnPEE || TnPMuMu) ptMinLepton = 28;

	// --electrons
		if (nElectrons > 1) {			// --electrons distributions for events with at least 2ele
			for (unsigned short j(0); j < nElectrons; j++){
				// doEleDistributionsPlots(electrons,j,8,0);
				// if (saveHEEPvariables_) {
				// 	int eleIdx = 4;
				// 	doElePlots(electrons,j,eleIdx,0);
				// 	if ( isInEB(electrons[j].v.Eta()) ) doElePlots(electrons,j,eleIdx,1);
				// 	if ( isInEE(electrons[j].v.Eta()) ) doElePlots(electrons,j,eleIdx,2);
				// }
				if (fabs(electrons[j].v.Eta()) < 2.4 && electrons[j].v.Pt() > ptMinLepton && electrons[j].passHEEPId) goodElectrons.push_back(electrons[j]);
			}

			for (unsigned short j(0); j < goodElectrons.size(); j++){
				doEleDistributionsPlots(electrons,j,9,0);
				if (saveHEEPvariables_) {
					int goodEleIdx = 5;
					doElePlots(goodElectrons,j,goodEleIdx,0);
					if ( isInEB(goodElectrons[j].v.Eta()) ) doElePlots(goodElectrons,j,goodEleIdx,1);
					if ( isInEE(goodElectrons[j].v.Eta()) ) doElePlots(goodElectrons,j,goodEleIdx,2);
				}
			}
		}
		// cout << "arriva a linea " << __LINE__ << endl;

	// --muons
		if (nMuons > 1) {    // --muons distributions for events with at least 2muons
			for (unsigned short i(0); i < nMuons; i++){ 
				doMuonDistributionsPlots(muons,i,10,0);
				// muon_dxy_histo->Fill( mu_dxy->at(i) );

				if (fabs(muons[i].v.Eta()) < 2.4 && muons[i].v.Pt() > ptMinLepton && muons[i].isHighPt) {
					goodMuons.push_back(muons[i]);
					goodMuonsRochCorr.push_back(muonsRochCorr[i]);
				}
			}

			for (unsigned short i(0); i < goodMuons.size(); i++){
				doMuonDistributionsPlots(goodMuons,i,11,0);
			}

			for (unsigned short i(0); i < nMuonsRochCor; i++){
				doMuonDistributionsPlots(muonsRochCorr,i,14,0);
			}

			for (unsigned short i(0); i < goodMuonsRochCorr.size(); i++){
				doMuonDistributionsPlots(goodMuonsRochCorr,i,15,0);
			}
		}
		// cout << "arriva a linea " << __LINE__ << endl;

	// --jets
		if (nJets > 1) {    // --jets distributions for events with at least 2jets
			for (unsigned short l(0); l < nJets; l++){
				// doJetsDistributionsPlots(jets,l,12,0);
				if (fabs(jets[l].v.Eta()) < 2.4 && jets[l].v.Pt() > 40 && jets[l].isTight) goodJets.push_back(jets[l]);
			}

			for (unsigned short l(0); l < goodJets.size(); l++){
				doJetsDistributionsPlots(jets,l,13,0);
			}
		}
		// cout << "arriva a linea " << __LINE__ << endl;	


		// for (unsigned short i(0); i < nElectrons; i++){
		// 	cout << "ele " << i << ": pt=" << electrons[i].v.Pt() << ", eta=" << electrons[i].v.Eta() << ", phi=" << electrons[i].v.Phi() << ", charge=" << electrons[i].charge << ", passHEEPId=" << electrons[i].passHEEPId << endl;						
		// }

		// for (unsigned short i(0); i < nMuons; i++){
		// 	cout << "muon " << i << ": pt=" << muons[i].v.Pt() << ", eta=" << muons[i].v.Eta() << ", phi=" << muons[i].v.Phi() << ", charge=" << muons[i].charge << endl;						
		// }

		// for (unsigned short i(0); i < nJets; i++){
		// 	cout << "jet " << i << ": pt=" << jets[i].v.Pt() << ", eta=" << jets[i].v.Eta() << ", phi=" << jets[i].v.Phi() << endl;						
		// }


		float ptMinLooseLlepton = 30;
		float ptMinLooseSLlepton = 28;
		if (signalEE) {
			ptMinLooseLlepton = 37;
			ptMinLooseSLlepton = 35;
		}
		if (signalMuMu) {
			ptMinLooseLlepton = 55;
			ptMinLooseSLlepton = 53;
		}


		if ( !eMuSideband ) {   
			// --Z using electrons
			if (nElectrons > 1) {
				unsigned int lIdx=0, slIdx=1; //sono ordinati in ordine decrescente in pt
				// cout << "leadingEle: pt=" << electrons[lIdx].v.Pt() << ", eta=" << electrons[lIdx].v.Eta() << ", phi=" << electrons[lIdx].v.Phi() << ", charge=" << electrons[lIdx].charge << ", passHEEPId=" << electrons[lIdx].passHEEPId << endl;
				// cout << "subLeadingEle: pt=" << electrons[slIdx].v.Pt() << ", eta=" << electrons[slIdx].v.Eta() << ", phi=" << electrons[slIdx].v.Phi() << ", charge=" << electrons[slIdx].charge << ", passHEEPId=" << electrons[slIdx].passHEEPId << endl;					
				if (electrons[lIdx].v.Pt() > ptMinLooseLlepton && electrons[slIdx].v.Pt() > ptMinLooseSLlepton) {
					nEventsWith2elePassingLoosePreselections++;
					// cout << "2elePassingLoosePreselections" << endl;
					if (electrons[lIdx].charge * electrons[slIdx].charge < 0) {
						nEventsWith2elePassingLoosePreselectionsAndCharge++;
						// cout << "2elePassingLoosePreselectionsAndCharge" << endl;	
						if (electrons[lIdx].passHEEPId && electrons[slIdx].passHEEPId) {  
							nEventsWith2elePassingLoosePreselectionsAndChargeAndHeepId++;
							TLorentzVector lEle, slEle;
							lEle.SetPtEtaPhiE(electrons[lIdx].v.Pt(), electrons[lIdx].v.Eta(), electrons[lIdx].v.Phi(), electrons[lIdx].v.E());
							slEle.SetPtEtaPhiE(electrons[slIdx].v.Pt(), electrons[slIdx].v.Eta(), electrons[slIdx].v.Phi(), electrons[slIdx].v.E());
							int singleEleIdx = 0;
							doZPlots(lEle, slEle, singleEleIdx);

							if (nJets > 1) {
								doZPlots(lEle, slEle, 4);
								nEventsFillingZtoEE++;
								// cout << "fillo ZtoEE " << endl;

								for (unsigned short j(0); j < nElectrons; j++){
									doEleDistributionsPlots(electrons,j,8,0);
									if (saveHEEPvariables_) {
										int eleIdx = 4;
										doElePlots(electrons,j,eleIdx,0);
										if ( isInEB(electrons[j].v.Eta()) ) doElePlots(electrons,j,eleIdx,1);
										if ( isInEE(electrons[j].v.Eta()) ) doElePlots(electrons,j,eleIdx,2);
									}
								}

								for (unsigned short l(0); l < nJets; l++) doJetsDistributionsPlots(jets,l,12,0);
								
							}
						}
					}
				}
			}

			// --Z using muons
			if (nMuons > 1) {
				unsigned int lIdx=0, slIdx=1; //sono ordinati in ordine decrescente in pt
				if (muons[lIdx].v.Pt() > ptMinLooseLlepton && muons[slIdx].v.Pt() > ptMinLooseSLlepton) {
					nEventsWith2muonsPassingLoosePreselections++;
					if (muons[lIdx].charge * muons[slIdx].charge < 0) {
						nEventsWith2muonsPassingLoosePreselectionsAndCharge++;	
						if (muons[lIdx].isHighPt && muons[slIdx].isHighPt) {  
							nEventsWith2muonsPassingLoosePreselectionsAndChargeAndHighPt++;
							TLorentzVector lMuon, slMuon;
							lMuon.SetPtEtaPhiE(muons[lIdx].v.Pt(), muons[lIdx].v.Eta(), muons[lIdx].v.Phi(), muons[lIdx].v.E());
							slMuon.SetPtEtaPhiE(muons[slIdx].v.Pt(), muons[slIdx].v.Eta(), muons[slIdx].v.Phi(), muons[slIdx].v.E());
							int singleMuonIdx = 1;
							doZPlots(lMuon, slMuon, singleMuonIdx);

							if (nJets > 1) doZPlots(lMuon, slMuon, 5);
						}
					}
				}
			}
		}
		// cout << "arriva a linea " << __LINE__ << endl;




	// --MultiLeptonMultiJet
		// cout << "nMultiLeptonMultiJets = " << nMultiLeptonMultiJets << endl;
		if (nMultiLeptonMultiJets < 1) continue;
		nEventsWithAtLeast1DLDJ++;

		//select MLMJ with right lepton pair
		vector<int> idxMLMJwithRightLeptonPair_vector;
		for (unsigned short j(0); j < nMultiLeptonMultiJets; j++){
			if ( ((signalEE || TnPEE) && multiLeptonMultiJets[j].isEEJJ) || ((signalMuMu || TnPMuMu) && multiLeptonMultiJets[j].isMMJJ) || (eMuSideband && multiLeptonMultiJets[j].isEMJJ) ) {
				// multiLeptonMultiJetsWithRightLeptonPair.push_back(multiLeptonMultiJets[j]);
				idxMLMJwithRightLeptonPair_vector.push_back(j);
			}
		}

		int nMLMJwithRightLeptonPair = idxMLMJwithRightLeptonPair_vector.size();	
		// cout << "nMLMJwithRightLeptonPair = " << nMLMJwithRightLeptonPair << endl;
		if (nMLMJwithRightLeptonPair < 1) continue;
		nEventsWithRightLeptonPair++;


		//choose mlmj candidate with sumPt >
		int idxDLDJ = idxMLMJwithRightLeptonPair_vector[0];
		// cout << "idxDLDJ=" << idxDLDJ << endl;

		if (nMLMJwithRightLeptonPair > 1) {
			for (unsigned short j(1); j < nMLMJwithRightLeptonPair; j++){		
				// cout << "sumPt mlmj[" << idxDLDJ << "] = " << multiLeptonMultiJets[idxDLDJ].multiLeptonMultiJet_sumPt << endl;
				// cout << "sumPt mlmj[" << idxMLMJwithRightLeptonPair_vector[j] << "] = " << multiLeptonMultiJets[idxMLMJwithRightLeptonPair_vector[j]].multiLeptonMultiJet_sumPt << endl;	
				if (multiLeptonMultiJets[idxMLMJwithRightLeptonPair_vector[j]].multiLeptonMultiJet_sumPt > multiLeptonMultiJets[idxDLDJ].multiLeptonMultiJet_sumPt) {
					idxDLDJ = idxMLMJwithRightLeptonPair_vector[j];
					// cout << "prendo idxDLDJ = " << idxMLMJwithRightLeptonPair_vector[j] << endl;
				}
			}
		} 
		// cout << "idxDLDJ = " << idxDLDJ << endl;


		if (signalEE || TnPEE) {
			if (leadingElectrons[idxDLDJ].v.Pt() > subLeadingElectrons[idxDLDJ].v.Pt()) nEventsWithLandSLleptonsRight++;
			if (leadingElectrons[idxDLDJ].v.Pt() <= subLeadingElectrons[idxDLDJ].v.Pt()) nEventsWithLandSLleptonsInverted++;
		}
		if (leadingJets[idxDLDJ].v.Pt() > subLeadingJets[idxDLDJ].v.Pt()) nEventsWithLandSLjetsRight++;
		if (leadingJets[idxDLDJ].v.Pt() <= subLeadingJets[idxDLDJ].v.Pt()) nEventsWithLandSLjetsInverted++;

		// cout << "leadingEle_MLMJ: pt=" << leadingElectrons[idxDLDJ].v.Pt() << ", eta=" << leadingElectrons[idxDLDJ].v.Eta() << ", phi=" << leadingElectrons[idxDLDJ].v.Phi() << ", charge=" << leadingElectrons[idxDLDJ].charge << ", passHEEPId=" << leadingElectrons[idxDLDJ].passHEEPId << endl;
		// cout << "subLeadingEle_MLMJ: pt=" << subLeadingElectrons[idxDLDJ].v.Pt() << ", eta=" << subLeadingElectrons[idxDLDJ].v.Eta() << ", phi=" << subLeadingElectrons[idxDLDJ].v.Phi() << ", charge=" << subLeadingElectrons[idxDLDJ].charge << ", passHEEPId=" << subLeadingElectrons[idxDLDJ].passHEEPId << endl;

		// cout << "leadingJet_MLMJ: pt=" << leadingJets[idxDLDJ].v.Pt() << ", eta=" << leadingJets[idxDLDJ].v.Eta() << ", phi=" << leadingJets[idxDLDJ].v.Phi() << endl;
		// cout << "subLeadingJet_MLMJ: pt=" << subLeadingJets[idxDLDJ].v.Pt() << ", eta=" << subLeadingJets[idxDLDJ].v.Eta() << ", phi=" << subLeadingJets[idxDLDJ].v.Phi() << endl;

		// if ( !(leadingElectrons[idxDLDJ].v.Pt() == electrons[0].v.Pt()) ) cout << "leadingEle_MLMJ pt = " << leadingElectrons[idxDLDJ].v.Pt() << ", ele0 pt = " << electrons[0].v.Pt() << endl;
		// if ( !(subLeadingElectrons[idxDLDJ].v.Pt() == electrons[1].v.Pt()) ) cout << "subLeadingEle_MLMJ pt = " << subLeadingElectrons[idxDLDJ].v.Pt() << ", ele1 pt = " << electrons[1].v.Pt() << endl;

		// if ( !(leadingJets[idxDLDJ].v.Pt() == jets[0].v.Pt()) ) cout << "leadingJet_MLMJ pt = " << leadingJets[idxDLDJ].v.Pt() << ", jet0 pt = " << jets[0].v.Pt() << endl;
		// if ( !(subLeadingJets[idxDLDJ].v.Pt() == jets[1].v.Pt()) ) cout << "subLeadingJet_MLMJ pt = " << subLeadingJets[idxDLDJ].v.Pt() << ", jet1 pt = " << jets[1].v.Pt() << endl;


		//Z->ee 
		if (TnPEE || (signalEE && leadingElectrons[idxDLDJ].v.Pt() > 37 && subLeadingElectrons[idxDLDJ].v.Pt() > 35) ) {  
			nEventsWithEEJJpassingLoosePreselections++; 
			// cout << "MLMJpassingLoosePreselections" << endl;
			// cout << "leadingEle_MLMJ: pt=" << leadingElectrons[idxDLDJ].v.Pt() << ", eta=" << leadingElectrons[idxDLDJ].v.Eta() << ", phi=" << leadingElectrons[idxDLDJ].v.Phi() << ", charge=" << leadingElectrons[idxDLDJ].charge << ", passHEEPId=" << leadingElectrons[idxDLDJ].passHEEPId << endl;
			// cout << "subLeadingEle_MLMJ: pt=" << subLeadingElectrons[idxDLDJ].v.Pt() << ", eta=" << subLeadingElectrons[idxDLDJ].v.Eta() << ", phi=" << subLeadingElectrons[idxDLDJ].v.Phi() << ", charge=" << subLeadingElectrons[idxDLDJ].charge << ", passHEEPId=" << subLeadingElectrons[idxDLDJ].passHEEPId << endl;

			if (leadingElectrons[idxDLDJ].charge * subLeadingElectrons[idxDLDJ].charge < 0 ) {	
				nEventsWithEEJJpassingLoosePreselectionsAndCharge++;		
				// cout << "MLMJpassingLoosePreselectionsAndCharge" << endl;
				if (leadingElectrons[idxDLDJ].passHEEPId && subLeadingElectrons[idxDLDJ].passHEEPId) {
					nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId++;
					int mlmjEleIdx = 2;
					doZPlots(multiLeptonMultiJets,idxDLDJ,mlmjEleIdx);
					nEventsFillingZtoEE_MLMJ++;
					// cout << "fillo ZtoEE_MLMJ " << endl;
				}
			}
		}
		
		//Z->mumu
		if (TnPMuMu || (signalMuMu && leadingMuons[idxDLDJ].v.Pt() > 55 && subLeadingMuons[idxDLDJ].v.Pt() > 53) ) {   
			nEventsWithMMJJpassingLoosePreselections++; 
			if (leadingMuons[idxDLDJ].charge * subLeadingMuons[idxDLDJ].charge < 0 ) { 
				nEventsWithMMJJpassingLoosePreselectionsAndCharge++;
				if (leadingMuons[idxDLDJ].isHighPt && subLeadingMuons[idxDLDJ].isHighPt) {
					nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt++;
					int mlmjMuonIdx = 3;
					doZPlots(multiLeptonMultiJets,idxDLDJ,mlmjMuonIdx);
				}
			}
		}




		int LleptonIdx=0, SLleptonIdx=1, LeleIdx=2, SLeleIdx=3, LmuIdx=4, SLmuIdx=5, LjetIdx=6, SLjetIdx=7;

		if (TnPEE || TnPMuMu) {	
			// --leadingEle
			if (multiLeptonMultiJets[idxDLDJ].isEEJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {			
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LleptonIdx,0); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LeleIdx,0);
			}

			// --subLeadingEle
			if (multiLeptonMultiJets[idxDLDJ].isEEJJ) {
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLleptonIdx,0); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLeleIdx,0);
			}

			// --leadingMuon
			if (multiLeptonMultiJets[idxDLDJ].isMMJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LleptonIdx,0); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LmuIdx,0);
			}

			// --subLeadingMuon
			if (multiLeptonMultiJets[idxDLDJ].isMMJJ) {
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLleptonIdx,0); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLmuIdx,0);
			}

			// --leadingJet & subLeadingJet
			doJetsDistributionsPlots(leadingJets,idxDLDJ,LjetIdx,0);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,SLjetIdx,0);

			// --MLMJ mass 
			bool isInEBEB = inEBEB(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));  
			bool isInEEEE = inEEEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));
			bool isInEBEE = inEBEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ)); 		

			doMassPlots(multiLeptonMultiJets,idxDLDJ,0,0);
			if (isInEBEB) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,1);
			if (isInEEEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,2);
			if (isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,3);
			if (isInEEEE || isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,4);
		}


		if (!(multiLeptonMultiJets[idxDLDJ].passPreselections)) continue;
		nEventsWithDLDJpassingPreselections++;


		// --leadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {
			nLeadingElePassingPreselections++;
			if (saveHEEPvariables_) {
				int leadingEleIdx = 0;
				doElePlots(leadingElectrons,idxDLDJ,leadingEleIdx,0);
				if ( isInEB(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,leadingEleIdx,1);
				if ( isInEE(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,leadingEleIdx,2);
			}
			// if ( leadingElectrons[idxDLDJ].isEcalDriven ) nLeadingElePassingIsEcalDriven++;
			// if ( fabs(leadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nLeadingElePassingdEtaIn++;
			// if ( fabs(leadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nLeadingElePassingdPhiIn++;
			// if ( (leadingElectrons[idxDLDJ].e2x5/leadingElectrons[idxDLDJ].e5x5) > 0.94 || (leadingElectrons[idxDLDJ].e1x5/leadingElectrons[idxDLDJ].e5x5) > 0.83 ) nLeadingElePassingE2x5OverE5x5++;
			// if ( leadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*leadingElectrons[idxDLDJ].v.Pt()) ) nLeadingElePassingEmHadDepth1Iso++;
			// if ( leadingElectrons[idxDLDJ].missingHits <=1 ) nLeadingElePassingMissingHits++;
			// if ( fabs(leadingElectrons[idxDLDJ].dxy)<0.02 ) nLeadingElePassingDxy++;
			// if ( leadingElectrons[idxDLDJ].passHEEPId ) nLeadingElePassingHeepId++;

			// if ( fabs(leadingElectrons[idxDLDJ].etaSC)<1.4442 ) {
			// 	nLeadingElePassingPreselectionsInEB++;
			// 	if ( leadingElectrons[idxDLDJ].isEcalDriven ) nLeadingEleInEBPassingIsEcalDriven++;
			// 	if ( fabs(leadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nLeadingEleInEBPassingdEtaIn++;
			// 	if ( fabs(leadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nLeadingEleInEBPassingdPhiIn++;
			// 	if ( (leadingElectrons[idxDLDJ].e2x5/leadingElectrons[idxDLDJ].e5x5) > 0.94 || (leadingElectrons[idxDLDJ].e1x5/leadingElectrons[idxDLDJ].e5x5) > 0.83 ) nLeadingEleInEBPassingE2x5OverE5x5++;
			// 	if ( leadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*leadingElectrons[idxDLDJ].v.Pt()) ) nLeadingEleInEBPassingEmHadDepth1Iso++;
			// 	if ( leadingElectrons[idxDLDJ].missingHits <=1 ) nLeadingEleInEBPassingMissingHits++;
			// 	if ( fabs(leadingElectrons[idxDLDJ].dxy)<0.02 ) nLeadingEleInEBPassingDxy++;
			// 	if ( leadingElectrons[idxDLDJ].passHEEPId ) nLeadingElePassingHeepIdInEB++;
			// }

			// if ( fabs(leadingElectrons[idxDLDJ].etaSC)>1.556 && fabs(leadingElectrons[idxDLDJ].etaSC)<2.5 ) {
			// 	nLeadingElePassingPreselectionsInEE++;
			// 	if ( leadingElectrons[idxDLDJ].isEcalDriven ) nLeadingEleInEEPassingIsEcalDriven++;
			// 	if ( fabs(leadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nLeadingEleInEEPassingdEtaIn++;
			// 	if ( fabs(leadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nLeadingEleInEEPassingdPhiIn++;
			// 	if ( (leadingElectrons[idxDLDJ].e2x5/leadingElectrons[idxDLDJ].e5x5) > 0.94 || (leadingElectrons[idxDLDJ].e1x5/leadingElectrons[idxDLDJ].e5x5) > 0.83 ) nLeadingEleInEEPassingE2x5OverE5x5++;
			// 	if ( leadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*leadingElectrons[idxDLDJ].v.Pt()) ) nLeadingEleInEEPassingEmHadDepth1Iso++;
			// 	if ( leadingElectrons[idxDLDJ].missingHits <=1 ) nLeadingEleInEEPassingMissingHits++;
			// 	if ( fabs(leadingElectrons[idxDLDJ].dxy)<0.02 ) nLeadingEleInEEPassingDxy++;
			// 	if ( leadingElectrons[idxDLDJ].passHEEPId ) nLeadingElePassingHeepIdInEE++;
			// }
		}


		// --subLeadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ) {
			nSubLeadingElePassingPreselections++;
			if (saveHEEPvariables_) {
				int subLeadingEleIdx = 1;
				doElePlots(subLeadingElectrons,idxDLDJ,subLeadingEleIdx,0);
				if ( isInEB(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,subLeadingEleIdx,1);
				if ( isInEE(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,subLeadingEleIdx,2);
			}
			// if ( subLeadingElectrons[idxDLDJ].isEcalDriven ) nSubLeadingElePassingIsEcalDriven++;
			// if ( fabs(subLeadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nSubLeadingElePassingdEtaIn++;
			// if ( fabs(subLeadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nSubLeadingElePassingdPhiIn++;
			// if ( (subLeadingElectrons[idxDLDJ].e2x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.94 || (subLeadingElectrons[idxDLDJ].e1x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.83 ) nSubLeadingElePassingE2x5OverE5x5++;
			// if ( subLeadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*subLeadingElectrons[idxDLDJ].v.Pt()) ) nSubLeadingElePassingEmHadDepth1Iso++;
			// if ( subLeadingElectrons[idxDLDJ].missingHits <=1 ) nSubLeadingElePassingMissingHits++;
			// if ( fabs(subLeadingElectrons[idxDLDJ].dxy)<0.02 ) nSubLeadingElePassingDxy++;
			// if ( subLeadingElectrons[idxDLDJ].passHEEPId ) nSubLeadingElePassingHeepId++;


			// if ( fabs(subLeadingElectrons[idxDLDJ].etaSC)<1.4442 ) {
			// 	nSubLeadingElePassingPreselectionsInEB++;
			// 	if ( subLeadingElectrons[idxDLDJ].isEcalDriven ) nSubLeadingEleInEBPassingIsEcalDriven++;
			// 	if ( fabs(subLeadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nSubLeadingEleInEBPassingdEtaIn++;
			// 	if ( fabs(subLeadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nSubLeadingEleInEBPassingdPhiIn++;
			// 	if ( (subLeadingElectrons[idxDLDJ].e2x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.94 || (subLeadingElectrons[idxDLDJ].e1x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.83 ) nSubLeadingEleInEBPassingE2x5OverE5x5++;
			// 	if ( subLeadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*subLeadingElectrons[idxDLDJ].v.Pt()) ) nSubLeadingEleInEBPassingEmHadDepth1Iso++;
			// 	if ( subLeadingElectrons[idxDLDJ].missingHits <=1 ) nSubLeadingEleInEBPassingMissingHits++;
			// 	if ( fabs(subLeadingElectrons[idxDLDJ].dxy)<0.02 ) nSubLeadingEleInEBPassingDxy++;
			// 	if ( subLeadingElectrons[idxDLDJ].passHEEPId ) nSubLeadingElePassingHeepIdInEB++;
			// }

			// if ( fabs(subLeadingElectrons[idxDLDJ].etaSC)>1.556 && fabs(subLeadingElectrons[idxDLDJ].etaSC)<2.5 ) {
			// 	nSubLeadingElePassingPreselectionsInEE++;
			// 	if ( subLeadingElectrons[idxDLDJ].isEcalDriven ) nSubLeadingEleInEEPassingIsEcalDriven++;
			// 	if ( fabs(subLeadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nSubLeadingEleInEEPassingdEtaIn++;
			// 	if ( fabs(subLeadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nSubLeadingEleInEEPassingdPhiIn++;
			// 	if ( (subLeadingElectrons[idxDLDJ].e2x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.94 || (subLeadingElectrons[idxDLDJ].e1x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.83 ) nSubLeadingEleInEEPassingE2x5OverE5x5++;
			// 	if ( subLeadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*subLeadingElectrons[idxDLDJ].v.Pt()) ) nSubLeadingEleInEEPassingEmHadDepth1Iso++;
			// 	if ( subLeadingElectrons[idxDLDJ].missingHits <=1 ) nSubLeadingEleInEEPassingMissingHits++;
			// 	if ( fabs(subLeadingElectrons[idxDLDJ].dxy)<0.02 ) nSubLeadingEleInEEPassingDxy++;
			// 	if ( subLeadingElectrons[idxDLDJ].passHEEPId ) nSubLeadingElePassingHeepIdInEE++;
			// }
		}


		bool passSelections = 0;

		if (leadingJets[idxDLDJ].isTight && subLeadingJets[idxDLDJ].isTight) {
			if (multiLeptonMultiJets[idxDLDJ].isEEJJ && leadingElectrons[idxDLDJ].passHEEPId && subLeadingElectrons[idxDLDJ].passHEEPId) passSelections = 1;
			if (multiLeptonMultiJets[idxDLDJ].isMMJJ && leadingMuons[idxDLDJ].isHighPt && subLeadingMuons[idxDLDJ].isHighPt) passSelections = 1;
			if (multiLeptonMultiJets[idxDLDJ].isEMJJ && leadingElectrons[idxDLDJ].passHEEPId && leadingMuons[idxDLDJ].isHighPt) passSelections = 1;			
		}

		if (!passSelections) continue;
		nEventsWithDLDJpassingSelections++;
		

		//definition of signal/control regions
		bool signalRegion=0, flavourSidebandCR=0, lowMlljjCR=0, lowMllCR=0;

		if (signalEE || signalMuMu) {
			if (multiLeptonMultiJets[idxDLDJ].isSignalRegion) {
				nEventsWithDLDJpassingSelectionsInSignalRegion++;
				signalRegion=true;
			}
			if (multiLeptonMultiJets[idxDLDJ].isLowMlljjCR) {
				nEventsWithDLDJpassingSelectionsInLowMlljjCR++;
				lowMlljjCR=true;
			}
			if (multiLeptonMultiJets[idxDLDJ].isLowMllCR) {
				nEventsWithDLDJpassingSelectionsInLowMllCR++;
				lowMllCR=true;
			}
		}

		if (eMuSideband && multiLeptonMultiJets[idxDLDJ].isSignalRegion) {
			nEventsWithDLDJpassingSelectionsInFlavourSidebandCR++;
			flavourSidebandCR=true;  
		}


		bool isInEBEB = inEBEB(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));  
		bool isInEEEE = inEEEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));
		bool isInEBEE = inEBEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ)); 		

		
		if (signalRegion) {		
			doMassPlots(multiLeptonMultiJets,idxDLDJ,0,0);
			if (isInEBEB) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,1);
			if (isInEEEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,2);
			if (isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,3);
			if (isInEEEE || isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,4);
		}
		if (lowMlljjCR) {
			doMassPlots(multiLeptonMultiJets,idxDLDJ,1,0);
			if (isInEBEB) doMassPlots(multiLeptonMultiJets,idxDLDJ,1,1);
			if (isInEEEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,1,2);
			if (isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,1,3);
			if (isInEEEE || isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,1,4);
		}
		if (lowMllCR) {
			doMassPlots(multiLeptonMultiJets,idxDLDJ,2,0);
			if (isInEBEB) doMassPlots(multiLeptonMultiJets,idxDLDJ,2,1);
			if (isInEEEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,2,2);
			if (isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,2,3);
			if (isInEEEE || isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,2,4);
		}

		if (flavourSidebandCR) {
			doMassPlots(multiLeptonMultiJets,idxDLDJ,0,0);
			if (isInEBEB) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,1);
			if (isInEEEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,2);
			if (isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,3);
			if (isInEEEE || isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,0,4);
		}



		int nLeadingLeptons=0, nSubLeadingLeptons=0;


		// --leadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {
			nLeadingLeptons++;  
			nEventsWithLeadingElePassingSelections++;

			if (signalRegion) {	
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LleptonIdx,0); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LeleIdx,0);
			}
			if (lowMlljjCR) {
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LleptonIdx,1); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LeleIdx,1);				
			}
			if (lowMllCR) {
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LleptonIdx,2); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LeleIdx,2);					
			}
			if (flavourSidebandCR) {
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LleptonIdx,0); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,LeleIdx,0);					
			}

			if (saveHEEPvariables_) {
				doElePlots(leadingElectrons,idxDLDJ,LeleIdx,0);
				if ( isInEB(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,LeleIdx,1);
				if ( isInEE(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,LeleIdx,2);	   
			}
		}


		// --subLeadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ) {
			nSubLeadingLeptons++;
			nEventsWithSubLeadingElePassingSelections++;

			if (signalRegion) {	
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLleptonIdx,0); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLeleIdx,0);
			}
			if (lowMlljjCR) {
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLleptonIdx,1); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLeleIdx,1);				
			}
			if (lowMllCR) {
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLleptonIdx,2); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLeleIdx,2);					
			}
			if (flavourSidebandCR) {
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLleptonIdx,0); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,SLeleIdx,0);					
			}

			if (saveHEEPvariables_) {
				doElePlots(subLeadingElectrons,idxDLDJ,SLeleIdx,0);
				if ( isInEB(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,SLeleIdx,1);
				if ( isInEE(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,SLeleIdx,2); 
			}
		}


		// --leadingMuon
		if (multiLeptonMultiJets[idxDLDJ].isMMJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {
			nLeadingLeptons++;
			nEventsWithLeadingMuonPassingSelections++;

			if (signalRegion) {	
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LleptonIdx,0); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LmuIdx,0);
			}
			if (lowMlljjCR) {
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LleptonIdx,1); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LmuIdx,1);				
			}
			if (lowMllCR) {
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LleptonIdx,2); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LmuIdx,2);					
			}
			if (flavourSidebandCR) {
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LleptonIdx,0); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,LmuIdx,0);					
			}
		}				


		// --subLeadingMuon
		if (multiLeptonMultiJets[idxDLDJ].isMMJJ) {
			nSubLeadingLeptons++;
			nEventsWithSubLeadingMuonPassingSelections++;

			if (signalRegion) {	
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLleptonIdx,0); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLmuIdx,0);
			}
			if (lowMlljjCR) {
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLleptonIdx,1); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLmuIdx,1);				
			}
			if (lowMllCR) {
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLleptonIdx,2); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLmuIdx,2);					
			}
			if (flavourSidebandCR) {
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLleptonIdx,0); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,SLmuIdx,0);					
			}	
		}


		// --leadingJet & subLeadingJet
		nEventsWithLeadingJetPassingSelections++;
		nEventsWithSubLeadingJetPassingSelections++;

		if (signalRegion) {	
			doJetsDistributionsPlots(leadingJets,idxDLDJ,LjetIdx,0);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,SLjetIdx,0);
		}
		if (lowMlljjCR) {
			doJetsDistributionsPlots(leadingJets,idxDLDJ,LjetIdx,1);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,SLjetIdx,1);		
		}
		if (lowMllCR) {
			doJetsDistributionsPlots(leadingJets,idxDLDJ,LjetIdx,2);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,SLjetIdx,2);				
		}
		if (flavourSidebandCR) {
			doJetsDistributionsPlots(leadingJets,idxDLDJ,LjetIdx,0);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,SLjetIdx,0);		
		}

		// nLeadingLeptons_histo -> Fill(nLeadingLeptons);
		// nSubLeadingLeptons_histo -> Fill(nSubLeadingLeptons);


	} //fine loop su eventi

	cout << "nEventsFillingZtoEE = " << nEventsFillingZtoEE << endl;
	cout << "nEventsFillingZtoEE_MLMJ = " << nEventsFillingZtoEE_MLMJ << endl;

	cout << "nEventsWith2elePassingLoosePreselections = " << nEventsWith2elePassingLoosePreselections << endl;
	cout << "nEventsWith2elePassingLoosePreselectionsAndCharge = " << nEventsWith2elePassingLoosePreselectionsAndCharge << endl;
	cout << "nEventsWith2elePassingLoosePreselectionsAndChargeAndHeepId = " << nEventsWith2elePassingLoosePreselectionsAndChargeAndHeepId << endl;

	cout << "nEventsWithEEJJpassingLoosePreselections = " << nEventsWithEEJJpassingLoosePreselections << endl;
	cout << "nEventsWithEEJJpassingLoosePreselectionsAndCharge = " << nEventsWithEEJJpassingLoosePreselectionsAndCharge << endl;
	cout << "nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId = " << nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId << endl;

	if (nEventsWithLandSLleptonsInverted != 0) cout << "nEventsWithLandSLleptonsInverted = " << nEventsWithLandSLleptonsInverted << endl;

}
// ******************************************************************************************




// **************** 
void makePlots_multiLeptonMultiJet::saveHistosAndOutputFile(TString& ouputdir){
	// setTDRStyle();
	// writeExtraText = true; 
	// extraText = "Preliminary"; 
	// lumi_13TeV = "0.6 fb^{-1}"; 

	gSystem->Exec(Form("mkdir -p %s",string(outputdir).data()));

	unsigned short numbOfHistograms = listOfHistograms.size();


	cout << "triggerName = " << triggerName << endl;
	TFile f(outputdir+"/"+"distributions_histos"+triggerName+".root","recreate");
	for (unsigned short l(0); l < numbOfHistograms; l++){
		string histoName = listOfHistograms[l]->GetName();
		if ( histoName.find(triggerName) != string::npos ) listOfHistograms[l]->Write(); 
	}


	cout << "numberRegions = " << numberRegions << endl;
	for (unsigned short n(0); n < numberRegions; n++){
		TFile f(outputdir+"/"+"distributions_histos"+regionName[n]+".root","recreate");
		for (unsigned short l(0); l < numbOfHistograms; l++){
			string histoName = listOfHistograms[l]->GetName();
			if ( histoName.find(regionName[n]) != string::npos ) listOfHistograms[l]->Write(); 
		}			
	}


	string RegionNameFile = "";
	if (signalEE) RegionNameFile = "signalEE";
	if (signalMuMu) RegionNameFile = "signalMuMu";
	if (eMuSideband) RegionNameFile = "flavourSideband";
	if (TnPEE) RegionNameFile = "TnPEE";
	if (TnPMuMu) RegionNameFile = "TnPMuMu";

	ofstream fOutput;
	fOutput.open(outputdir+"/nEvents"+RegionNameFile+".dat");

	fOutput << "nEvents = " << nEvents << " \n";
	fOutput << "nEventsWithAtLeast1DLDJ = " << nEventsWithAtLeast1DLDJ << " \n";
	fOutput << "nEventsWithRightLeptonPair = " << nEventsWithRightLeptonPair << " \n";

	fOutput << "nEventsWithLandSLleptonsRight = " << nEventsWithLandSLleptonsRight << " \n";
	fOutput << "nEventsWithLandSLleptonsInverted = " << nEventsWithLandSLleptonsInverted << " \n";
	fOutput << "nEventsWithLandSLjetsRight = " << nEventsWithLandSLjetsRight << " \n";
	fOutput << "nEventsWithLandSLjetsInverted = " << nEventsWithLandSLjetsInverted << " \n";

	fOutput << "nEventsWithDLDJpassingPreselections = " << nEventsWithDLDJpassingPreselections << " \n";
	fOutput << "nEventsWithDLDJpassingSelections = " << nEventsWithDLDJpassingSelections << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInSignalRegion = " << nEventsWithDLDJpassingSelectionsInSignalRegion << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInLowMlljjCR = " << nEventsWithDLDJpassingSelectionsInLowMlljjCR << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInLowMllCR = " << nEventsWithDLDJpassingSelectionsInLowMllCR << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInFlavourSidebandCR = " << nEventsWithDLDJpassingSelectionsInFlavourSidebandCR << " \n" << " \n";
	fOutput << "nEventsWithLeadingElePassingSelections = " << nEventsWithLeadingElePassingSelections << " \n";
	fOutput << "nEventsWithSubLeadingElePassingSelections = " << nEventsWithSubLeadingElePassingSelections << " \n";
	fOutput << "nEventsWithLeadingMuonPassingSelections = " << nEventsWithLeadingMuonPassingSelections << " \n";
	fOutput << "nEventsWithSubLeadingMuonPassingSelections = " << nEventsWithSubLeadingMuonPassingSelections << " \n";
	fOutput << "nEventsWithLeadingJetPassingSelections = " << nEventsWithLeadingJetPassingSelections << " \n";
	fOutput << "nEventsWithSubLeadingJetPassingSelections = " << nEventsWithSubLeadingJetPassingSelections << " \n" << " \n";

	fOutput << "nEventsWithEEJJpassingLoosePreselections = " << nEventsWithEEJJpassingLoosePreselections << " \n";
	fOutput << "nEventsWithEEJJpassingLoosePreselectionsAndCharge = " << nEventsWithEEJJpassingLoosePreselectionsAndCharge << " \n";
	fOutput << "nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId = " << nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId << " \n";
	fOutput << "nEventsWithMMJJpassingLoosePreselections = " << nEventsWithMMJJpassingLoosePreselections << " \n";
	fOutput << "nEventsWithMMJJpassingLoosePreselectionsAndCharge = " << nEventsWithMMJJpassingLoosePreselectionsAndCharge << " \n";
	fOutput << "nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt = " << nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt << " \n" << " \n";


	fOutput << "nLeadingElePassingPreselections = " << nLeadingElePassingPreselections << " \n";
	fOutput << "nLeadingElePassingIsEcalDriven = " << nLeadingElePassingIsEcalDriven << " \n";
	fOutput << "nLeadingElePassingdEtaIn = " << nLeadingElePassingdEtaIn << " \n";
	fOutput << "nLeadingElePassingdPhiIn = " << nLeadingElePassingdPhiIn << " \n";
	fOutput << "nLeadingElePassingE2x5OverE5x5 = " << nLeadingElePassingE2x5OverE5x5 << " \n";
	fOutput << "nLeadingElePassingEmHadDepth1Iso = " << nLeadingElePassingEmHadDepth1Iso << " \n";
	fOutput << "nLeadingElePassingMissingHits = " << nLeadingElePassingMissingHits << " \n";
	fOutput << "nLeadingElePassingDxy = " << nLeadingElePassingDxy << " \n";
	fOutput << "nLeadingElePassingHeepId = " << nLeadingElePassingHeepId << " \n" << " \n";

	fOutput << "nLeadingElePassingPreselectionsInEB = " << nLeadingElePassingPreselectionsInEB << " \n";
	fOutput << "nLeadingEleInEBPassingIsEcalDriven = " << nLeadingEleInEBPassingIsEcalDriven << " \n";
	fOutput << "nLeadingEleInEBPassingdEtaIn = " << nLeadingEleInEBPassingdEtaIn << " \n";
	fOutput << "nLeadingEleInEBPassingdPhiIn = " << nLeadingEleInEBPassingdPhiIn << " \n";
	fOutput << "nLeadingEleInEBPassingE2x5OverE5x5 = " << nLeadingEleInEBPassingE2x5OverE5x5 << " \n";
	fOutput << "nLeadingEleInEBPassingEmHadDepth1Iso = " << nLeadingEleInEBPassingEmHadDepth1Iso << " \n";
	fOutput << "nLeadingEleInEBPassingMissingHits = " << nLeadingEleInEBPassingMissingHits << " \n";
	fOutput << "nLeadingEleInEBPassingDxy = " << nLeadingEleInEBPassingDxy << " \n";
	fOutput << "nLeadingElePassingHeepIdInEB = " << nLeadingElePassingHeepIdInEB << " \n" << " \n";

	fOutput << "nLeadingElePassingPreselectionsInEE = " << nLeadingElePassingPreselectionsInEE << " \n";
	fOutput << "nLeadingEleInEEPassingIsEcalDriven = " << nLeadingEleInEEPassingIsEcalDriven << " \n";
	fOutput << "nLeadingEleInEEPassingdEtaIn = " << nLeadingEleInEEPassingdEtaIn << " \n";
	fOutput << "nLeadingEleInEEPassingdPhiIn = " << nLeadingEleInEEPassingdPhiIn << " \n";
	fOutput << "nLeadingEleInEEPassingE2x5OverE5x5 = " << nLeadingEleInEEPassingE2x5OverE5x5 << " \n";
	fOutput << "nLeadingEleInEEPassingEmHadDepth1Iso = " << nLeadingEleInEEPassingEmHadDepth1Iso << " \n";
	fOutput << "nLeadingEleInEEPassingMissingHits = " << nLeadingEleInEEPassingMissingHits << " \n";
	fOutput << "nLeadingEleInEEPassingDxy = " << nLeadingEleInEEPassingDxy << " \n";
	fOutput << "nLeadingElePassingHeepIdInEE = " << nLeadingElePassingHeepIdInEE << " \n" << " \n";


	fOutput << "nSubLeadingElePassingPreselections = " << nSubLeadingElePassingPreselections << " \n";
	fOutput << "nSubLeadingElePassingIsEcalDriven = " << nSubLeadingElePassingIsEcalDriven << " \n";
	fOutput << "nSubLeadingElePassingdEtaIn = " << nSubLeadingElePassingdEtaIn << " \n";
	fOutput << "nSubLeadingElePassingdPhiIn = " << nSubLeadingElePassingdPhiIn << " \n";
	fOutput << "nSubLeadingElePassingE2x5OverE5x5 = " << nSubLeadingElePassingE2x5OverE5x5 << " \n";
	fOutput << "nSubLeadingElePassingEmHadDepth1Iso = " << nSubLeadingElePassingEmHadDepth1Iso << " \n";
	fOutput << "nSubLeadingElePassingMissingHits = " << nSubLeadingElePassingMissingHits << " \n";
	fOutput << "nSubLeadingElePassingDxy = " << nSubLeadingElePassingDxy << " \n";
	fOutput << "nSubLeadingElePassingHeepId = " << nSubLeadingElePassingHeepId << " \n" << " \n";

	fOutput << "nSubLeadingElePassingPreselectionsInEB = " << nSubLeadingElePassingPreselectionsInEB << " \n";
	fOutput << "nSubLeadingEleInEBPassingIsEcalDriven = " << nSubLeadingEleInEBPassingIsEcalDriven << " \n";
	fOutput << "nSubLeadingEleInEBPassingdEtaIn = " << nSubLeadingEleInEBPassingdEtaIn << " \n";
	fOutput << "nSubLeadingEleInEBPassingdPhiIn = " << nSubLeadingEleInEBPassingdPhiIn << " \n";
	fOutput << "nSubLeadingEleInEBPassingE2x5OverE5x5 = " << nSubLeadingEleInEBPassingE2x5OverE5x5 << " \n";
	fOutput << "nSubLeadingEleInEBPassingEmHadDepth1Iso = " << nSubLeadingEleInEBPassingEmHadDepth1Iso << " \n";
	fOutput << "nSubLeadingEleInEBPassingMissingHits = " << nSubLeadingEleInEBPassingMissingHits << " \n";
	fOutput << "nSubLeadingEleInEBPassingDxy = " << nSubLeadingEleInEBPassingDxy << " \n";
	fOutput << "nSubLeadingElePassingHeepIdInEB = " << nSubLeadingElePassingHeepIdInEB << " \n" << " \n";

	fOutput << "nSubLeadingElePassingPreselectionsInEE = " << nSubLeadingElePassingPreselectionsInEE << " \n";
	fOutput << "nSubLeadingEleInEEPassingIsEcalDriven = " << nSubLeadingEleInEEPassingIsEcalDriven << " \n";
	fOutput << "nSubLeadingEleInEEPassingdEtaIn = " << nSubLeadingEleInEEPassingdEtaIn << " \n";
	fOutput << "nSubLeadingEleInEEPassingdPhiIn = " << nSubLeadingEleInEEPassingdPhiIn << " \n";
	fOutput << "nSubLeadingEleInEEPassingE2x5OverE5x5 = " << nSubLeadingEleInEEPassingE2x5OverE5x5 << " \n";
	fOutput << "nSubLeadingEleInEEPassingEmHadDepth1Iso = " << nSubLeadingEleInEEPassingEmHadDepth1Iso << " \n";
	fOutput << "nSubLeadingEleInEEPassingMissingHits = " << nSubLeadingEleInEEPassingMissingHits << " \n";
	fOutput << "nSubLeadingEleInEEPassingDxy = " << nSubLeadingEleInEEPassingDxy << " \n";
	fOutput << "nSubLeadingElePassingHeepIdInEE = " << nSubLeadingElePassingHeepIdInEE << " \n" << " \n";

	// fOutput << " = " <<  << " \n";	
	fOutput.close();

}
// ******************************************************************************************



makePlots_multiLeptonMultiJet::makePlots_multiLeptonMultiJet(TString filename_, TString outputdir_, bool MC_, bool signalEE_, bool signalMuMu_, bool eMuSideband_, bool TnPEE_, bool TnPMuMu_, int nToDivideWeight_):
	filename(filename_), outputdir(outputdir_), MC(MC_), signalEE(signalEE_), signalMuMu(signalMuMu_), eMuSideband(eMuSideband_), TnPEE(TnPEE_), TnPMuMu(TnPMuMu_), nToDivideWeight(nToDivideWeight_)
{
	fChain = new TChain("", "");

	TFile *f = TFile::Open(filename);
	if(!f){
		cerr << "Failed to open file " << filename << ".\n";
	} else {
		cout << "Reading input file " << filename << "\n";
		TString treePath = filename + "/analysisTree/event";
		if (fChain) fChain->Add(treePath);
	}

	SetHistos();
	Loop();
	saveHistosAndOutputFile(outputdir);
}



void runMakePlots_multiLeptonMultiJet() {

	// string inputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/miniTrees-Moriond17/MLMJwithSingleObjects/selectedTrigger/";   
	// string outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger/";
	// string outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger_L-SL_swapped/";

	string inputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/miniTrees-Moriond17/MLMJwithSingleObjects/selectedTrigger/rightPt_L-SL/";   
	// string outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger/rightPt_L-SL/";
	string outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/Distributions-Moriond17/MLMJwithSingleObjects/selectedTrigger/rightPt_L-SL/plotsForEleZ/";

	bool signalEE = false;  //true for DoubleEG (if not TnP)   
	bool signalMuMu = false;  //true for SingleMuon (if not TnP)
	bool eMuSideband = false; //true for MuonEG

	bool TnPEE = true; //true for DoubleEG 
	bool TnPMuMu = false; //true for SingleMuon 

	string triggerName;

	if (signalEE) triggerName = "signalEE/";
	if (signalMuMu) triggerName = "signalMuMu/";
	if (eMuSideband) triggerName = "eMuSideband/";
	if (TnPEE) triggerName = "TnPEE/";
	if (TnPMuMu) triggerName = "TnPMuMu/";

	inputDir = inputDir+triggerName;
	outputDir = outputDir+triggerName;


	// makePlots_multiLeptonMultiJet("/afs/cern.ch/user/g/gnegro/work/NuAnalysis-Moriond17/CMSSW_8_0_26_patch1/src/dafne/output.root", "prova", false, signalEE, false, eMuSideband, TnPEE, false, 1);

	// makePlots_multiLeptonMultiJet(inputDir+"DoubleEG/output_DoubleEG_miniTree.root", outputDir+"DoubleEG", false, signalEE, false, false, TnPEE, false, 1);
	// makePlots_multiLeptonMultiJet(inputDir+"DoubleEG/output_DoubleEG_miniTree.root", "prova", false, signalEE, false, false, TnPEE, false, 1);
	// makePlots_multiLeptonMultiJet("outputDoubleEG.root", "provaDoubleEG", false, signalEE, false, false, TnPEE, false, 1);
	// makePlots_multiLeptonMultiJet("outputDY.root", "provaDY", false, signalEE, false, false, TnPEE, false, 1);

	// makePlots_multiLeptonMultiJet(inputDir+"SingleEle/output_SingleEle_miniTree.root", outputDir+"SingleEle", false, signalEE, false, false, TnPEE, false, 1);

	// makePlots_multiLeptonMultiJet(inputDir+"SingleMuon/output_SingleMuon_miniTree.root", outputDir+"SingleMuon", false, false, signalMuMu, false, false, TnPMuMu, 1);

	// makePlots_multiLeptonMultiJet(inputDir+"MuonEG/output_MuonEG_miniTree.root", outputDir+"MuonEG", false, false, false, true, false, false, 1);


	makePlots_multiLeptonMultiJet(inputDir+"TTJets/output_TTJets_miniTree.root", outputDir+"TTJets", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 1);

	makePlots_multiLeptonMultiJet(inputDir+"WJets/output_WJetsToLNu_miniTree.root", outputDir+"WJets", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 1);

	makePlots_multiLeptonMultiJet(inputDir+"WZ/output_WZ_miniTree.root", outputDir+"WZ", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 1);

	makePlots_multiLeptonMultiJet(inputDir+"ZZ/output_ZZ_miniTree.root", outputDir+"ZZ", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 1);

	makePlots_multiLeptonMultiJet(inputDir+"WW/output_WW_miniTree.root", outputDir+"WW", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 1);


	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY50/output_DYJetsToLL_Pt-50To100_amcatnloFXFX-v3_dafne-Moriond17_miniTree.root", outputDir+"DYJetsPtBinned/DY50-1", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 2);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY50/output_DYJetsToLL_Pt-50To100_amcatnloFXFX-ext3-v1_dafne-Moriond17_miniTree.root", outputDir+"DYJetsPtBinned/DY50-2", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 2);

	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY100/output_DYJetsToLL_Pt-100To250_amcatnloFXFX-v2_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY100-1", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY100/output_DYJetsToLL_Pt-100To250_amcatnloFXFX-ext1-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY100-2", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY100/output_DYJetsToLL_Pt-100To250_amcatnloFXFX-ext2-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY100-3", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);

	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY250/output_DYJetsToLL_Pt-250To400_amcatnloFXFX-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY250-1", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY250/output_DYJetsToLL_Pt-250To400_amcatnloFXFX-ext1-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY250-2", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY250/output_DYJetsToLL_Pt-250To400_amcatnloFXFX-ext2-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY250-3", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);

	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY400/output_DYJetsToLL_Pt-400To650_amcatnloFXFX-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY400-1", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY400/output_DYJetsToLL_Pt-400To650_amcatnloFXFX-ext1-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY400-2", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// // makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY400/output_DYJetsToLL_Pt-400To650_amcatnloFXFX-ext2-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY400-3", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);

	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY650/output_DYJetsToLL_Pt-650ToInf_amcatnloFXFX-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY650-1", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY650/output_DYJetsToLL_Pt-650ToInf_amcatnloFXFX-ext1-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY650-2", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);
	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/DY650/output_DYJetsToLL_Pt-650ToInf_amcatnloFXFX-ext2-v1_dafne-Moriond17_miniTree.root", outputDir+"/DYJetsPtBinned/DY650-3", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 3);


	// makePlots_multiLeptonMultiJet(inputDir+"DYinclusive-amcatnlo/output_DYJetsToLL-amcatnloFXFX-inclusive_miniTree.root", outputDir+"/DYinclusive-amcatnlo", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 1);
	// makePlots_multiLeptonMultiJet(inputDir+"DYinclusive-madgraph/output_DYJetsToLL-madgraphMLM_inclusive_miniTree.root", outputDir+"/DYinclusive-madgraph", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMuMu, 1);

}



#ifndef __CLING__
int main() {
  runMakePlots_multiLeptonMultiJet();
  return 0;
}
#endif