#include "dafne/Validation/interface/makePlots_multiLeptonMultiJet.h"
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
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
	ele_HEEPBitMapValues = 0;
    ele_passTightId = 0;
    ele_passMediumId = 0;
    ele_passLooseId = 0;
    ele_passVetoId = 0;
    ele_passMVATightId = 0;
    ele_passMVAMediumId = 0;
	ele_idmva = 0;
	ele_iso = 0;
	ele_dz = 0;
	ele_d0 = 0;
	ele_isMatchedToGen = 0;
	ele_charge = 0;
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
	mu_e = 0;
	mu_pt = 0;
	mu_eta = 0;
	mu_phi = 0;
	mu_iso = 0;
	mu_PFiso = 0;
	mu_isTight = 0;
	mu_isMedium = 0;
	mu_isLoose = 0;
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
	jet_bdiscriminant = 0;
	jet_partonFlavour = 0;
	jet_hadronFlavour = 0;
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
	leadingEle_HEEPBitMapValues = 0;
    leadingEle_passTightId = 0;
    leadingEle_passMediumId = 0;
    leadingEle_passLooseId = 0;
    leadingEle_passVetoId = 0;
    leadingEle_passMVATightId = 0;
    leadingEle_passMVAMediumId = 0;
    leadingEle_idmva = 0;
    leadingEle_iso = 0;
    leadingEle_isMatchedToGen = 0;
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
	subLeadingEle_passHEEPId = 0;
	subLeadingEle_HEEPBitMapValues = 0;
   	subLeadingEle_passTightId = 0;
    subLeadingEle_passMediumId = 0;
    subLeadingEle_passLooseId = 0;
    subLeadingEle_passVetoId = 0;
    subLeadingEle_passMVATightId = 0;
    subLeadingEle_passMVAMediumId = 0;
    subLeadingEle_idmva = 0;
   	subLeadingEle_iso = 0;
   	subLeadingEle_isMatchedToGen = 0;
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
	leadingMuon_iso = 0;
	leadingMuon_PFiso = 0;
	leadingMuon_isHighPt = 0;
	leadingMuon_isTight = 0;
    leadingMuon_isMedium = 0;
    leadingMuon_isLoose = 0;
    leadingMuon_isMatchedToGen = 0;
    leadingMuon_dz = 0;
    leadingMuon_dxy = 0;
    leadingMuon_RochCor = 0;
    subLeadingMuon_iso = 0;
    subLeadingMuon_PFiso = 0;
   	subLeadingMuon_isHighPt = 0;
   	subLeadingMuon_isTight = 0;
    subLeadingMuon_isMedium = 0;
    subLeadingMuon_isLoose = 0;
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
	fChain->SetBranchAddress("passEEJJhlt", &passEEJJhlt, &b_passEEJJhlt);
	fChain->SetBranchAddress("passMMJJhlt", &passMMJJhlt, &b_passMMJJhlt);
	fChain->SetBranchAddress("passEMJJhlt", &passEMJJhlt, &b_passEMJJhlt);
	fChain->SetBranchAddress("passTandPEEhlt", &passTandPEEhlt, &b_passTandPEEhlt);
	fChain->SetBranchAddress("passTandPMMhlt", &passTandPMMhlt, &b_passTandPMMhlt);
	fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
	fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
	fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
	fChain->SetBranchAddress("ele_e", &ele_e, &b_ele_e);
	fChain->SetBranchAddress("ele_pt", &ele_pt, &b_ele_pt);
	fChain->SetBranchAddress("ele_eta", &ele_eta, &b_ele_eta);
	fChain->SetBranchAddress("ele_phi", &ele_phi, &b_ele_phi);
	fChain->SetBranchAddress("ele_passHEEPId", &ele_passHEEPId, &b_ele_passHEEPId);
    fChain->SetBranchAddress("ele_HEEPBitMapValues", &ele_HEEPBitMapValues, &b_ele_HEEPBitMapValues);
    fChain->SetBranchAddress("ele_passTightId", &ele_passTightId, &b_ele_passTightId);
    fChain->SetBranchAddress("ele_passMediumId", &ele_passMediumId, &b_ele_passMediumId);
    fChain->SetBranchAddress("ele_passLooseId", &ele_passLooseId, &b_ele_passLooseId);
    fChain->SetBranchAddress("ele_passVetoId", &ele_passVetoId, &b_ele_passVetoId);
    fChain->SetBranchAddress("ele_passMVATightId", &ele_passMVATightId, &b_ele_passMVATightId);
    fChain->SetBranchAddress("ele_passMVAMediumId", &ele_passMVAMediumId, &b_ele_passMVAMediumId);
	fChain->SetBranchAddress("ele_idmva", &ele_idmva, &b_ele_idmva);
	fChain->SetBranchAddress("ele_iso", &ele_iso, &b_ele_iso);
	fChain->SetBranchAddress("ele_dz", &ele_dz, &b_ele_dz);
	fChain->SetBranchAddress("ele_d0", &ele_d0, &b_ele_d0);
	fChain->SetBranchAddress("ele_isMatchedToGen", &ele_isMatchedToGen, &b_ele_isMatchedToGen);
	fChain->SetBranchAddress("ele_charge", &ele_charge, &b_ele_charge);
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
	fChain->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
	fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
	fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
	fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
	fChain->SetBranchAddress("mu_iso", &mu_iso, &b_mu_iso);
	fChain->SetBranchAddress("mu_PFiso", &mu_PFiso, &b_mu_PFiso);
	fChain->SetBranchAddress("mu_isTight", &mu_isTight, &b_mu_isTight);
	fChain->SetBranchAddress("mu_isMedium", &mu_isMedium, &b_mu_isMedium);
	fChain->SetBranchAddress("mu_isLoose", &mu_isLoose, &b_mu_isLoose);
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
	fChain->SetBranchAddress("jet_bdiscriminant", &jet_bdiscriminant, &b_jet_bdiscriminant);
	fChain->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
	fChain->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
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
	fChain->SetBranchAddress("leadingEle_HEEPBitMapValues", &leadingEle_HEEPBitMapValues, &b_leadingEle_HEEPBitMapValues);
    fChain->SetBranchAddress("leadingEle_passTightId", &leadingEle_passTightId, &b_leadingEle_passTightId);
    fChain->SetBranchAddress("leadingEle_passMediumId", &leadingEle_passMediumId, &b_leadingEle_passMediumId);
    fChain->SetBranchAddress("leadingEle_passLooseId", &leadingEle_passLooseId, &b_leadingEle_passLooseId);
    fChain->SetBranchAddress("leadingEle_passVetoId", &leadingEle_passVetoId, &b_leadingEle_passVetoId);
    fChain->SetBranchAddress("leadingEle_passMVATightId", &leadingEle_passMVATightId, &b_leadingEle_passMVATightId);
    fChain->SetBranchAddress("leadingEle_passMVAMediumId", &leadingEle_passMVAMediumId, &b_leadingEle_passMVAMediumId);
    fChain->SetBranchAddress("leadingEle_idmva", &leadingEle_idmva, &b_leadingEle_idmva);
    fChain->SetBranchAddress("leadingEle_iso", &leadingEle_iso, &b_leadingEle_iso);
    fChain->SetBranchAddress("leadingEle_isMatchedToGen", &leadingEle_isMatchedToGen, &b_leadingEle_isMatchedToGen);
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
	fChain->SetBranchAddress("subLeadingEle_passHEEPId", &subLeadingEle_passHEEPId, &b_subLeadingEle_passHEEPId);
	fChain->SetBranchAddress("subLeadingEle_HEEPBitMapValues", &subLeadingEle_HEEPBitMapValues, &b_subLeadingEle_HEEPBitMapValues);
    fChain->SetBranchAddress("subLeadingEle_passTightId", &subLeadingEle_passTightId, &b_subLeadingEle_passTightId);
    fChain->SetBranchAddress("subLeadingEle_passMediumId", &subLeadingEle_passMediumId, &b_subLeadingEle_passMediumId);
    fChain->SetBranchAddress("subLeadingEle_passLooseId", &subLeadingEle_passLooseId, &b_subLeadingEle_passLooseId);
    fChain->SetBranchAddress("subLeadingEle_passVetoId", &subLeadingEle_passVetoId, &b_subLeadingEle_passVetoId);
    fChain->SetBranchAddress("subLeadingEle_passMVATightId", &subLeadingEle_passMVATightId, &b_subLeadingEle_passMVATightId);
    fChain->SetBranchAddress("subLeadingEle_passMVAMediumId", &subLeadingEle_passMVAMediumId, &b_subLeadingEle_passMVAMediumId);
    fChain->SetBranchAddress("subLeadingEle_idmva", &subLeadingEle_idmva, &b_subLeadingEle_idmva);
    fChain->SetBranchAddress("subLeadingEle_iso", &subLeadingEle_iso, &b_subLeadingEle_iso);
    fChain->SetBranchAddress("subLeadingEle_isMatchedToGen", &subLeadingEle_isMatchedToGen, &b_subLeadingEle_isMatchedToGen);
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
    fChain->SetBranchAddress("leadingMuon_iso", &leadingMuon_iso, &b_leadingMuon_iso);
    fChain->SetBranchAddress("leadingMuon_PFiso", &leadingMuon_PFiso, &b_leadingMuon_PFiso);
	fChain->SetBranchAddress("leadingMuon_isHighPt", &leadingMuon_isHighPt, &b_leadingMuon_isHighPt);
 	fChain->SetBranchAddress("leadingMuon_isTight", &leadingMuon_isTight, &b_leadingMuon_isTight);
    fChain->SetBranchAddress("leadingMuon_isMedium", &leadingMuon_isMedium, &b_leadingMuon_isMedium);
    fChain->SetBranchAddress("leadingMuon_isLoose", &leadingMuon_isLoose, &b_leadingMuon_isLoose);
    fChain->SetBranchAddress("leadingMuon_isMatchedToGen", &leadingMuon_isMatchedToGen, &b_leadingMuon_isMatchedToGen);
    fChain->SetBranchAddress("leadingMuon_dz", &leadingMuon_dz, &b_leadingMuon_dz);
    fChain->SetBranchAddress("leadingMuon_dxy", &leadingMuon_dxy, &b_leadingMuon_dxy);
    fChain->SetBranchAddress("leadingMuon_RochCor", &leadingMuon_RochCor, &b_leadingMuon_RochCor);
    fChain->SetBranchAddress("subLeadingMuon_iso", &subLeadingMuon_iso, &b_subLeadingMuon_iso);
    fChain->SetBranchAddress("subLeadingMuon_PFiso", &subLeadingMuon_PFiso, &b_subLeadingMuon_PFiso);
   	fChain->SetBranchAddress("subLeadingMuon_isHighPt", &subLeadingMuon_isHighPt, &b_subLeadingMuon_isHighPt);
   	fChain->SetBranchAddress("subLeadingMuon_isTight", &subLeadingMuon_isTight, &b_subLeadingMuon_isTight);
   	fChain->SetBranchAddress("subLeadingMuon_isMedium", &subLeadingMuon_isMedium, &b_subLeadingMuon_isMedium);
    fChain->SetBranchAddress("subLeadingMuon_isLoose", &subLeadingMuon_isLoose, &b_subLeadingMuon_isLoose);
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

	rho_histo = newTH1D("rho_histo", "rho_histo", "rho", 100, 0., 100.);
	nvtx_histo = newTH1D("nvtx_histo", "nvtx_histo", "nvtx", 100, 0., 100.);
	vtx_x_histo = newTH1D("vtx_x_histo", "vtx_x_histo", "vtx_x", 100, -0.3, 0.3);
	vtx_y_histo = newTH1D("vtx_y_histo", "vtx_y_histo", "vtx_y", 100, -0.3, 0.3);
	vtx_z_histo = newTH1D("vtx_z_histo", "vtx_z_histo", "vtx_z", 100, -20., 20.);

	muon_dxy_histo = newTH1D("muon_dxy_histo", "muon_dxy_histo", "|dxy| muons", 100, 0., 0.5);

	nLeadingLeptons_histo = newTH1D("nLeadingLeptons_histo", "nLeadingLeptons_histo", "nLeadingLepton", 10, 0, 10); 
	nSubLeadingLeptons_histo = newTH1D("nSubLeadingLeptons_histo", "nSubLeadingLeptons_histo", "nSubLeadingLepton", 10, 0, 10); 


	for (unsigned short i(0); i < 6; i++) {
		pt_histo[i][0] = newTH1D("pt_"+objName[i]+"_histo", "pt_"+objName[i]+"_histo", "p_{T} [GeV] "+objName[i], 200, 0., 1000.);
		eta_histo[i][0] = newTH1D("eta_"+objName[i]+"_histo", "eta_"+objName[i]+"_histo", "#eta "+objName[i], 100, -2.5, 2.5);
		phi_histo[i][0] = newTH1D("phi_"+objName[i]+"_histo", "phi_"+objName[i]+"_histo", "#phi "+objName[i], 100, -3.5, 3.5);		
	}

	for (unsigned short i(14); i < 16; i++) {
		pt_histo[i][0] = newTH1D("pt_"+objName[i]+"_histo", "pt_"+objName[i]+"_histo", "p_{T} [GeV] "+objName[i], 200, 0., 1000.);
		eta_histo[i][0] = newTH1D("eta_"+objName[i]+"_histo", "eta_"+objName[i]+"_histo", "#eta "+objName[i], 100, -2.5, 2.5);
		phi_histo[i][0] = newTH1D("phi_"+objName[i]+"_histo", "phi_"+objName[i]+"_histo", "#phi "+objName[i], 100, -3.5, 3.5);		
	}


	for (unsigned short l(0); l < 4; l++) {

		for (unsigned short i(6); i < 14; i++) {
			pt_histo[i][l] = newTH1D("pt_"+objName[i]+regionName[l]+"_histo", "pt_"+objName[i]+regionName[l]+"_histo", "p_{T} [GeV] "+objName[i]+regionName[l], 200, 0., 1000.);
			eta_histo[i][l] = newTH1D("eta_"+objName[i]+regionName[l]+"_histo", "eta_"+objName[i]+regionName[l]+"_histo", "#eta "+objName[i]+regionName[l], 100, -2.5, 2.5);
			phi_histo[i][l] = newTH1D("phi_"+objName[i]+regionName[l]+"_histo", "phi_"+objName[i]+regionName[l]+"_histo", "#phi "+objName[i]+regionName[l], 100, -3.5, 3.5);		
		} 

		for (unsigned short n(0); n < 5; n++) {
			mass_dldj_histo[l][n] = newTH1D("mass_multiLeptonMultiJets"+regionName[l]+etaMassName[n]+"_histo", "mass_multiLeptonMultiJets"+regionName[l]+etaMassName[n]+"_histo", "m_{lljj}"+regionName[l], 300, 0., 6000.);
			mass_dl_histo[l][n] = newTH1D("mass_dileptons"+regionName[l]+etaMassName[n]+"_histo", "mass_dileptons"+regionName[l]+etaMassName[n]+"_histo", "m_{ll}"+regionName[l], 100, 0., 2000.);
			mass_dj_histo[l][n] = newTH1D("mass_dijets"+regionName[l]+etaMassName[n]+"_histo", "mass_dijets"+regionName[l]+etaMassName[n]+"_histo", "m_{jj}"+regionName[l], 100, 0., 2000.);		
			mass_djLl_histo[l][n] = newTH1D("mass_dijetsLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "mass_dijetsLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "m_{jjl_{L}}"+regionName[l], 200, 0., 4000.);	
			mass_djSLl_histo[l][n] = newTH1D("mass_dijetsSubLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "mass_dijetsSubLeadingLepton"+regionName[l]+etaMassName[n]+"_histo", "m_{jjl_{SL}}"+regionName[l], 200, 0., 4000.);
		}

		pt_dl_histo[l] = newTH1D("pt_dileptons_"+regionName[l]+"_histo", "pt_dileptons_"+regionName[l]+"_histo", "diLepton p_{T} [GeV]", 100, 0., 1000.);

	}


	for (unsigned short j(0); j < 6; j++) {
		for (unsigned short m(0); m < 3; m++) {
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

	for (unsigned short n(0); n < 5; n++) {
		for (unsigned short z(0); z < 4; z++) {
			Zmass_histo[z][n] = newTH1D("Z"+Zname[z]+"_mass"+etaMassName[n]+"_histo", "Z"+Zname[z]+"_mass"+etaMassName[n]+"_histo", "m(Z"+Zname[z]+") [GeV/c^{2}]", 200, 0, 200); 
		}
	}



}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::getElectrons(vector<eleStruct>& electrons, vector<eleStruct>& goodElectrons){
	unsigned short nTotElectrons(ele_e->size());
	for (unsigned short i(0); i < nTotElectrons; i++){
		eleStruct ele(ele_pt->at(i), 
						ele_eta->at(i), 						
						ele_phi->at(i), 
						ele_e->at(i), 
						ele_passHEEPId->at(i),
						ele_HEEPBitMapValues->at(i),
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
		electrons.push_back(ele);
		if (fabs(ele.v.Eta()) < 2.4 && ele.v.Pt() > 53 && ele.passHEEPId) goodElectrons.push_back(ele);
	} //End of loop over all the electrons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getLeadingElectrons(vector<eleStruct>& leadingElectrons){
	unsigned short nTotLeadingElectrons(leadingEle_etaSC->size());
	for (unsigned short i(0); i < nTotLeadingElectrons; i++){
		eleStruct ele(leadingLepton_pt->at(i),
						leadingLepton_eta->at(i),
						leadingLepton_phi->at(i),
						leadingLepton_e->at(i), 
						leadingEle_passHEEPId->at(i),
						leadingEle_HEEPBitMapValues->at(i),
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
		leadingElectrons.push_back(ele);
	} //End of loop over all the leading electrons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getSubLeadingElectrons(vector<eleStruct>& subLeadingElectrons){
	unsigned short nTotSubLeadingElectrons(subLeadingEle_etaSC->size());
	for (unsigned short i(0); i < nTotSubLeadingElectrons; i++){
		eleStruct ele(subLeadingLepton_pt->at(i),
						subLeadingLepton_eta->at(i),
						subLeadingLepton_phi->at(i),
						subLeadingLepton_e->at(i), 
						subLeadingEle_passHEEPId->at(i),
						subLeadingEle_HEEPBitMapValues->at(i),
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
		subLeadingElectrons.push_back(ele);
	} //End of loop over all the subLeading electrons
}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::getMuons(vector<muonStruct>& muons, vector<muonStruct>& goodMuons, vector<muonStruct>& muonsRochCorr, vector<muonStruct>& goodMuonsRochCorr){
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
		if (fabs(mu.v.Eta()) < 2.4 && mu.v.Pt() > 53 && mu.isHighPt) goodMuons.push_back(mu);

		mu.v.SetPtEtaPhiE(mu.v.Pt()*mu.rochCor, mu.v.Eta(), mu.v.Phi(), mu.v.E());
		muonsRochCorr.push_back(mu);
		if (fabs(mu.v.Eta()) < 2.4 && mu.v.Pt() > 53 && mu.isHighPt) goodMuonsRochCorr.push_back(mu);
	}//End of loop over all the muons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getLeadingMuons(vector<muonStruct>& leadingMuons, vector<muonStruct>& leadingMuonsRochCorr){
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

		mu.v.SetPtEtaPhiE(mu.v.Pt()*mu.rochCor, mu.v.Eta(), mu.v.Phi(), mu.v.E());
		leadingMuonsRochCorr.push_back(mu);
	}//End of loop over all the leading muons
}
// ******************************************************************************************


// **************** 
void makePlots_multiLeptonMultiJet::getSubLeadingMuons(vector<muonStruct>& subLeadingMuons, vector<muonStruct>& subLeadingMuonsRochCorr){
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
		
		mu.v.SetPtEtaPhiE(mu.v.Pt()*mu.rochCor, mu.v.Eta(), mu.v.Phi(), mu.v.E());
		subLeadingMuonsRochCorr.push_back(mu);
	}//End of loop over all the subleading muons
}
// ******************************************************************************************



// **************** 
void makePlots_multiLeptonMultiJet::getJets(vector<jetStruct>& jets, vector<jetStruct>& goodJets){
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
		if (fabs(jet.v.Eta()) < 2.4 && jet.v.Pt() > 40 && jet.isTight) goodJets.push_back(jet);
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
void makePlots_multiLeptonMultiJet::doZmassPlots(vector<multiLeptonMultiJetStruct>& multiLeptonMultiJets, const int mlmjIdx, const int variableIdx, const int etaIdx) { 	
	float mZ = multiLeptonMultiJets[mlmjIdx].diLepton_invMass;
	Zmass_histo[variableIdx][etaIdx]->Fill(mZ,w);
}
// ******************************************************************************************


void makePlots_multiLeptonMultiJet::doZmassPlots(TLorentzVector Llepton, TLorentzVector SLlepton, const int variableIdx) { 	
	float mZ = (Llepton+SLlepton).M();
	bool isInEBEB = inEBEB(Llepton.Eta(), SLlepton.Eta()); 
	bool isInEEEE = inEEEE(Llepton.Eta(), SLlepton.Eta());
	bool isInEBEE = inEBEE(Llepton.Eta(), SLlepton.Eta()); 

	Zmass_histo[variableIdx][0]->Fill(mZ,w);
	if (isInEBEB) Zmass_histo[variableIdx][1]->Fill(mZ,w);  
	if (isInEEEE) Zmass_histo[variableIdx][2]->Fill(mZ,w);  
	if (isInEBEE) Zmass_histo[variableIdx][3]->Fill(mZ,w);  
	if (isInEEEE || isInEBEE) Zmass_histo[variableIdx][4]->Fill(mZ,w);  
}






// **************** 
void makePlots_multiLeptonMultiJet::Loop(){
	//--- Initialize the tree branches ---
	Init();
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	cout << "nentries = " << nentries << endl;

	// for (Long64_t i=0; i<300; i+=1) {
	for (Long64_t i=0; i<nentries; i+=1) {
		Long64_t ientry = LoadTree(i);
		if (ientry < 0) break;
		// cout << "evento " << ientry << endl;
		nEvents++;

		if(fChain->GetEntry(i) == 0){
			std::cerr << "Failed to read Tree entry " << i << "!\n";
			continue;
		}


	// --Trigger
		bool passTrigger = 0;

		if (passEEJJhlt) nEventsPassingEEJJhlt++;
		if (passMMJJhlt) nEventsPassingMMJJhlt++;
		if (passEMJJhlt) nEventsPassingEMJJhlt++;
		if (passTandPEEhlt) nEventsPassingTandPEEhlt++;
		if (passTandPMMhlt) nEventsPassingTandPMMhlt++;

		if (signalEE) passTrigger = passEEJJhlt;   
		if (signalMuMu) passTrigger = passMMJJhlt;
		if (eMuSideband) passTrigger = passEMJJhlt;
		if (TnPEE) passTrigger = passTandPEEhlt;
		if (TnPMM) passTrigger = passTandPMMhlt;

		if (!passTrigger) continue;  //se evento non passa il trigger passo a quello successivo
		nEventsPassingTrigger++;



	// --Collections
		vector<eleStruct> electrons, goodElectrons, leadingElectrons, subLeadingElectrons;
		vector<muonStruct> muons, muonsRochCorr, goodMuons, goodMuonsRochCorr, leadingMuons, leadingMuonsRochCorr, subLeadingMuons, subLeadingMuonsRochCorr;
		vector<jetStruct> jets, goodJets, leadingJets, subLeadingJets;
		vector<multiLeptonMultiJetStruct> multiLeptonMultiJets;

		// cout << "get electrons" << endl;
		getElectrons(electrons, goodElectrons);

		// cout << "get muons" << endl;
		getMuons(muons, goodMuons, muonsRochCorr, goodMuonsRochCorr);

		// cout << "get jets" << endl;
		getJets(jets, goodJets);

		// cout << "get mlmj" << endl;
		getMultiLeptonMultiJets(multiLeptonMultiJets);

		// cout << "get leadingEle" << endl;
		getLeadingElectrons(leadingElectrons);
		// cout << "get subLeadingEle" << endl;
		getSubLeadingElectrons(subLeadingElectrons);

		// cout << "get leadingMuon" << endl;
		getLeadingMuons(leadingMuons, leadingMuonsRochCorr);
		// cout << "get subLeadingMuon" << endl;
		getSubLeadingMuons(subLeadingMuons, subLeadingMuonsRochCorr);

		// cout << "get leadingJet" << endl;
		getLeadingJets(leadingJets);
		// cout << "get subLeadingJet" << endl;
		getSubLeadingJets(subLeadingJets);

		unsigned short nElectrons = electrons.size();    
		unsigned short nGoodElectrons = goodElectrons.size();
		unsigned short nMuons = muons.size();
		unsigned short nGoodMuons = goodMuons.size();
		unsigned short nMuonsRochCor = muonsRochCorr.size();
		unsigned short nGoodMuonsRochCor = goodMuonsRochCorr.size();		
		unsigned short nJets = jets.size();
		unsigned short nGoodJets = goodJets.size();
		unsigned short nMultiLeptonMultiJets = multiLeptonMultiJets.size();
		unsigned short nLeadingEle = leadingElectrons.size();
		unsigned short nSubLeadingEle = subLeadingElectrons.size();
		unsigned short nLeadingMuons = leadingMuons.size();
		unsigned short nSubLeadingMuons = subLeadingMuons.size();
		unsigned short nLeadingJets = leadingJets.size();
		unsigned short nSubLeadingJets = subLeadingJets.size();

		if ( !(nLeadingEle == nMultiLeptonMultiJets) ) cout << "nLeadingEle != nMultiLeptonMultiJets" << endl;
		if ( !(nSubLeadingEle == nMultiLeptonMultiJets) ) cout << "nSubLeadingEle != nMultiLeptonMultiJets" << endl;
		if ( !(nLeadingMuons == nMultiLeptonMultiJets) ) cout << "nLeadingMuons != nMultiLeptonMultiJets" << endl;
		if ( !(nSubLeadingMuons == nMultiLeptonMultiJets) ) cout << "nSubLeadingMuons != nMultiLeptonMultiJets" << endl;
		if ( !(nLeadingJets == nMultiLeptonMultiJets) ) cout << "nLeadingJets != nMultiLeptonMultiJets" << endl;
		if ( !(nSubLeadingJets == nMultiLeptonMultiJets) ) cout << "nSubLeadingJets != nMultiLeptonMultiJets" << endl;
		if ( !(nMuonsRochCor == nMuons) ) cout << "nMuonsRochCor != nMuons" << endl;

		// if ( !(nGoodMuonsRochCor == nGoodMuons) ) cout << "nGoodMuonsRochCor != nGoodMuons" << endl;
		// puo essere diverso perche nell'if uso pt cambiato

		nOfElectrons += nElectrons;
		nOfGoodElectrons += nGoodElectrons;
		nOfMuons += nMuons;
		nOfGoodMuons += nGoodMuons;
		nOfGoodMuonsWithRochCor += nGoodMuonsRochCor;
		nOfJets += nJets;
		nOfGoodJets += nGoodJets;



		w = weight;
		// cout << "weight = " << w << endl;


	// --vtx
		rho_histo->Fill(rho,w);
		nvtx_histo->Fill(nvtx,w);

		for (unsigned short i(0); i < nvtx; i++){
			vtx_x_histo->Fill(vtx_x->at(i),w);
			vtx_y_histo->Fill(vtx_y->at(i),w);
			vtx_z_histo->Fill(vtx_z->at(i),w);
		}


	// --electrons
		for (unsigned short j(0); j < nElectrons; j++){
			doEleDistributionsPlots(electrons,j,0,0);
			doElePlots(electrons,j,0,0);
			if ( isInEB(electrons[j].v.Eta()) ) doElePlots(electrons,j,0,1);
			if ( isInEE(electrons[j].v.Eta()) ) doElePlots(electrons,j,0,2);
		}

		for (unsigned short j(0); j < nGoodElectrons; j++){
			doEleDistributionsPlots(goodElectrons,j,1,0);
			doElePlots(goodElectrons,j,1,0);
			if ( isInEB(goodElectrons[j].v.Eta()) ) doElePlots(goodElectrons,j,1,1);
			if ( isInEE(goodElectrons[j].v.Eta()) ) doElePlots(goodElectrons,j,1,2);
		}


	// --muons
		for (unsigned short i(0); i < nMuons; i++){
			doMuonDistributionsPlots(muons,i,2,0);
			muon_dxy_histo->Fill( mu_dxy->at(i) );
		}

		for (unsigned short i(0); i < nGoodMuons; i++){
			doMuonDistributionsPlots(goodMuons,i,3,0);
		}

		for (unsigned short i(0); i < nMuonsRochCor; i++){
			doMuonDistributionsPlots(muonsRochCorr,i,14,0);
		}

		for (unsigned short i(0); i < nGoodMuonsRochCor; i++){
			doMuonDistributionsPlots(goodMuonsRochCorr,i,15,0);
		}


	// --jets
		for (unsigned short l(0); l < nJets; l++){
			doJetsDistributionsPlots(jets,l,4,0);
		}

		for (unsigned short l(0); l < nGoodJets; l++){
			doJetsDistributionsPlots(goodJets,l,5,0);
		}


	// --Zmass using electrons
		if (TnPEE) {   
			if (nElectrons < 2) continue;
			unsigned int lIdx=0, slIdx=1; //sono ordinati in ordine decrescente in pt

			if (electrons[lIdx].v.Pt() < 35 || electrons[lIdx].v.Pt() < 35) continue; 
			nEventsWith2elePassingLoosePreselections++;

			if (electrons[lIdx].charge * electrons[slIdx].charge < 0) {
				nEventsWith2elePassingLoosePreselectionsAndCharge++;	

				if (electrons[lIdx].passHEEPId && electrons[slIdx].passHEEPId) {  
					nEventsWith2elePassingLoosePreselectionsAndChargeAndHeepId++;

					TLorentzVector lEle, slEle;
					lEle.SetPtEtaPhiE(electrons[lIdx].v.Pt(), electrons[lIdx].v.Eta(), electrons[lIdx].v.Phi(), electrons[lIdx].v.E());
					slEle.SetPtEtaPhiE(electrons[slIdx].v.Pt(), electrons[slIdx].v.Eta(), electrons[slIdx].v.Phi(), electrons[slIdx].v.E());

					doZmassPlots(lEle, slEle, 0);
				}
			}
		}

	// --Zmass using muons
		if (TnPMM) {   
			if (nMuons < 2) continue;
			unsigned int lIdx=0, slIdx=1; //sono ordinati in ordine decrescente in pt

			if (muons[lIdx].v.Pt() < 35 || muons[lIdx].v.Pt() < 35) continue; 
			nEventsWith2muonsPassingLoosePreselections++;

			if (muons[lIdx].charge * muons[slIdx].charge < 0) {
				nEventsWith2muonsPassingLoosePreselectionsAndCharge++;	

				if (muons[lIdx].isHighPt && muons[slIdx].isHighPt) {  
					nEventsWith2muonsPassingLoosePreselectionsAndChargeAndHighPt++;

					TLorentzVector lMuon, slMuon;
					lMuon.SetPtEtaPhiE(muons[lIdx].v.Pt(), muons[lIdx].v.Eta(), muons[lIdx].v.Phi(), muons[lIdx].v.E());
					slMuon.SetPtEtaPhiE(muons[slIdx].v.Pt(), muons[slIdx].v.Eta(), muons[slIdx].v.Phi(), muons[slIdx].v.E());

					doZmassPlots(lMuon, slMuon, 1);
				}
			}
		}




	// --MultiLeptonMultiJet
		//choose mlmj candidate
		if (nMultiLeptonMultiJets < 1) continue;
		nEventsWithAtLeast1DLDJ++;

		int idxDLDJ = 0;

		if (nMultiLeptonMultiJets > 1) {
			for (unsigned short j(1); j < nMultiLeptonMultiJets; j++){
				// cout << "sumPt mlmj[" << idxDLDJ << "] = " << multiLeptonMultiJets[idxDLDJ].multiLeptonMultiJet_sumPt << endl;
				// cout << "sumPt mlmj[" << j << "] = " << multiLeptonMultiJets[j].multiLeptonMultiJet_sumPt << endl;				
				if (multiLeptonMultiJets[j].multiLeptonMultiJet_sumPt > multiLeptonMultiJets[idxDLDJ].multiLeptonMultiJet_sumPt) {
					idxDLDJ = j;
					// cout << "prendo idxDLDJ = " << j << endl;
				}
			}
		} 
		// cout << "idxDLDJ = " << idxDLDJ << endl;


		//select right lepton pair
		if ( (signalEE || TnPEE) && (!(multiLeptonMultiJets[idxDLDJ].isEEJJ)) ) continue;
		if ( (signalMuMu || TnPMM) && (!(multiLeptonMultiJets[idxDLDJ].isMMJJ)) ) continue;
		if ( eMuSideband && !(multiLeptonMultiJets[idxDLDJ].isEMJJ) ) continue;
		nEventsWithRightLeptonPair++;


		bool isInEBEB = inEBEB(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));  //usare eta jets per djets mass ?
		bool isInEEEE = inEEEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));
		bool isInEBEE = inEBEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ)); 			


	// --Zmass
		if (TnPEE) {   //Z->ee 
			nEventsWithEEJJpassingLoosePreselections++; //=nEventsWithRightLeptonPair
			if (leadingElectrons[idxDLDJ].charge * subLeadingElectrons[idxDLDJ].charge < 0 ) {	
				nEventsWithEEJJpassingLoosePreselectionsAndCharge++;		
				if (leadingElectrons[idxDLDJ].passHEEPId && subLeadingElectrons[idxDLDJ].passHEEPId) {
					nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId++;
					doZmassPlots(multiLeptonMultiJets,idxDLDJ,2,0);
					if (isInEBEB) doZmassPlots(multiLeptonMultiJets,idxDLDJ,2,1);
					if (isInEEEE) doZmassPlots(multiLeptonMultiJets,idxDLDJ,2,2);
					if (isInEBEE) doZmassPlots(multiLeptonMultiJets,idxDLDJ,2,3);
					if (isInEEEE || isInEBEE) doZmassPlots(multiLeptonMultiJets,idxDLDJ,2,4);
				}
			}
		}

		if (TnPMM) {   //Z->mumu
			nEventsWithMMJJpassingLoosePreselections++; //=nEventsWithRightLeptonPair
			if (leadingMuons[idxDLDJ].charge * subLeadingMuons[idxDLDJ].charge < 0 ) { 
				nEventsWithMMJJpassingLoosePreselectionsAndCharge++;
				if (leadingMuons[idxDLDJ].isHighPt && subLeadingMuons[idxDLDJ].isHighPt) {
					nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt++;
					doZmassPlots(multiLeptonMultiJets,idxDLDJ,3,0);
					if (isInEBEB) doZmassPlots(multiLeptonMultiJets,idxDLDJ,3,1);
					if (isInEEEE) doZmassPlots(multiLeptonMultiJets,idxDLDJ,3,2);
					if (isInEBEE) doZmassPlots(multiLeptonMultiJets,idxDLDJ,3,3);
					if (isInEEEE || isInEBEE) doZmassPlots(multiLeptonMultiJets,idxDLDJ,3,4);
				}
			}
		}
	// --


		if (!(multiLeptonMultiJets[idxDLDJ].passPreselections)) continue;
		nEventsWithDLDJpassingPreselections++;


		// --leadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {
			nLeadingElePassingPreselections++;
			doElePlots(leadingElectrons,idxDLDJ,2,0);
			if ( isInEB(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,2,1);
			if ( isInEE(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,2,2);

			if ( leadingElectrons[idxDLDJ].isEcalDriven ) nLeadingElePassingIsEcalDriven++;
			if ( fabs(leadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nLeadingElePassingdEtaIn++;
			if ( fabs(leadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nLeadingElePassingdPhiIn++;
			if ( (leadingElectrons[idxDLDJ].e2x5/leadingElectrons[idxDLDJ].e5x5) > 0.94 || (leadingElectrons[idxDLDJ].e1x5/leadingElectrons[idxDLDJ].e5x5) > 0.83 ) nLeadingElePassingE2x5OverE5x5++;
			if ( leadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*leadingElectrons[idxDLDJ].v.Pt()) ) nLeadingElePassingEmHadDepth1Iso++;
			if ( leadingElectrons[idxDLDJ].missingHits <=1 ) nLeadingElePassingMissingHits++;
			if ( fabs(leadingElectrons[idxDLDJ].dxy)<0.02 ) nLeadingElePassingDxy++;
			if ( leadingElectrons[idxDLDJ].passHEEPId ) nLeadingElePassingHeepId++;

			if ( fabs(leadingElectrons[idxDLDJ].etaSC)<1.4442 ) {
				nLeadingElePassingPreselectionsInEB++;
				if ( leadingElectrons[idxDLDJ].isEcalDriven ) nLeadingEleInEBPassingIsEcalDriven++;
				if ( fabs(leadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nLeadingEleInEBPassingdEtaIn++;
				if ( fabs(leadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nLeadingEleInEBPassingdPhiIn++;
				if ( (leadingElectrons[idxDLDJ].e2x5/leadingElectrons[idxDLDJ].e5x5) > 0.94 || (leadingElectrons[idxDLDJ].e1x5/leadingElectrons[idxDLDJ].e5x5) > 0.83 ) nLeadingEleInEBPassingE2x5OverE5x5++;
				if ( leadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*leadingElectrons[idxDLDJ].v.Pt()) ) nLeadingEleInEBPassingEmHadDepth1Iso++;
				if ( leadingElectrons[idxDLDJ].missingHits <=1 ) nLeadingEleInEBPassingMissingHits++;
				if ( fabs(leadingElectrons[idxDLDJ].dxy)<0.02 ) nLeadingEleInEBPassingDxy++;
				if ( leadingElectrons[idxDLDJ].passHEEPId ) nLeadingElePassingHeepIdInEB++;
			}

			if ( fabs(leadingElectrons[idxDLDJ].etaSC)>1.556 && fabs(leadingElectrons[idxDLDJ].etaSC)<2.5 ) {
				nLeadingElePassingPreselectionsInEE++;
				if ( leadingElectrons[idxDLDJ].isEcalDriven ) nLeadingEleInEEPassingIsEcalDriven++;
				if ( fabs(leadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nLeadingEleInEEPassingdEtaIn++;
				if ( fabs(leadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nLeadingEleInEEPassingdPhiIn++;
				if ( (leadingElectrons[idxDLDJ].e2x5/leadingElectrons[idxDLDJ].e5x5) > 0.94 || (leadingElectrons[idxDLDJ].e1x5/leadingElectrons[idxDLDJ].e5x5) > 0.83 ) nLeadingEleInEEPassingE2x5OverE5x5++;
				if ( leadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*leadingElectrons[idxDLDJ].v.Pt()) ) nLeadingEleInEEPassingEmHadDepth1Iso++;
				if ( leadingElectrons[idxDLDJ].missingHits <=1 ) nLeadingEleInEEPassingMissingHits++;
				if ( fabs(leadingElectrons[idxDLDJ].dxy)<0.02 ) nLeadingEleInEEPassingDxy++;
				if ( leadingElectrons[idxDLDJ].passHEEPId ) nLeadingElePassingHeepIdInEE++;
			}
		}


		// --subLeadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ) {
			nSubLeadingElePassingPreselections++;
			doElePlots(subLeadingElectrons,idxDLDJ,3,0);
			if ( isInEB(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,3,1);
			if ( isInEE(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,3,2);

			if ( subLeadingElectrons[idxDLDJ].isEcalDriven ) nSubLeadingElePassingIsEcalDriven++;
			if ( fabs(subLeadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nSubLeadingElePassingdEtaIn++;
			if ( fabs(subLeadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nSubLeadingElePassingdPhiIn++;
			if ( (subLeadingElectrons[idxDLDJ].e2x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.94 || (subLeadingElectrons[idxDLDJ].e1x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.83 ) nSubLeadingElePassingE2x5OverE5x5++;
			if ( subLeadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*subLeadingElectrons[idxDLDJ].v.Pt()) ) nSubLeadingElePassingEmHadDepth1Iso++;
			if ( subLeadingElectrons[idxDLDJ].missingHits <=1 ) nSubLeadingElePassingMissingHits++;
			if ( fabs(subLeadingElectrons[idxDLDJ].dxy)<0.02 ) nSubLeadingElePassingDxy++;
			if ( subLeadingElectrons[idxDLDJ].passHEEPId ) nSubLeadingElePassingHeepId++;


			if ( fabs(subLeadingElectrons[idxDLDJ].etaSC)<1.4442 ) {
				nSubLeadingElePassingPreselectionsInEB++;
				if ( subLeadingElectrons[idxDLDJ].isEcalDriven ) nSubLeadingEleInEBPassingIsEcalDriven++;
				if ( fabs(subLeadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nSubLeadingEleInEBPassingdEtaIn++;
				if ( fabs(subLeadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nSubLeadingEleInEBPassingdPhiIn++;
				if ( (subLeadingElectrons[idxDLDJ].e2x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.94 || (subLeadingElectrons[idxDLDJ].e1x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.83 ) nSubLeadingEleInEBPassingE2x5OverE5x5++;
				if ( subLeadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*subLeadingElectrons[idxDLDJ].v.Pt()) ) nSubLeadingEleInEBPassingEmHadDepth1Iso++;
				if ( subLeadingElectrons[idxDLDJ].missingHits <=1 ) nSubLeadingEleInEBPassingMissingHits++;
				if ( fabs(subLeadingElectrons[idxDLDJ].dxy)<0.02 ) nSubLeadingEleInEBPassingDxy++;
				if ( subLeadingElectrons[idxDLDJ].passHEEPId ) nSubLeadingElePassingHeepIdInEB++;
			}

			if ( fabs(subLeadingElectrons[idxDLDJ].etaSC)>1.556 && fabs(subLeadingElectrons[idxDLDJ].etaSC)<2.5 ) {
				nSubLeadingElePassingPreselectionsInEE++;
				if ( subLeadingElectrons[idxDLDJ].isEcalDriven ) nSubLeadingEleInEEPassingIsEcalDriven++;
				if ( fabs(subLeadingElectrons[idxDLDJ].dEtaIn)<0.004 ) nSubLeadingEleInEEPassingdEtaIn++;
				if ( fabs(subLeadingElectrons[idxDLDJ].dPhiIn)<0.006 ) nSubLeadingEleInEEPassingdPhiIn++;
				if ( (subLeadingElectrons[idxDLDJ].e2x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.94 || (subLeadingElectrons[idxDLDJ].e1x5/subLeadingElectrons[idxDLDJ].e5x5) > 0.83 ) nSubLeadingEleInEEPassingE2x5OverE5x5++;
				if ( subLeadingElectrons[idxDLDJ].EmHadDepth1Iso<(0.28*rho + 2 + 0.03*subLeadingElectrons[idxDLDJ].v.Pt()) ) nSubLeadingEleInEEPassingEmHadDepth1Iso++;
				if ( subLeadingElectrons[idxDLDJ].missingHits <=1 ) nSubLeadingEleInEEPassingMissingHits++;
				if ( fabs(subLeadingElectrons[idxDLDJ].dxy)<0.02 ) nSubLeadingEleInEEPassingDxy++;
				if ( subLeadingElectrons[idxDLDJ].passHEEPId ) nSubLeadingElePassingHeepIdInEE++;
			}
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
			doMassPlots(multiLeptonMultiJets,idxDLDJ,3,0);
			if (isInEBEB) doMassPlots(multiLeptonMultiJets,idxDLDJ,3,1);
			if (isInEEEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,3,2);
			if (isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,3,3);
			if (isInEEEE || isInEBEE) doMassPlots(multiLeptonMultiJets,idxDLDJ,3,4);
		}



		int nLeadingLeptons=0, nSubLeadingLeptons=0;

		// --leadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {
			nLeadingLeptons++;  
			nEventsWithLeadingElePassingSelections++;

			if (signalRegion) {	
	 			doEleDistributionsPlots(leadingElectrons,idxDLDJ,6,0); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,8,0);
			}
			if (lowMlljjCR) {
	 			doEleDistributionsPlots(leadingElectrons,idxDLDJ,6,1); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,8,1);				
			}
			if (lowMllCR) {
	 			doEleDistributionsPlots(leadingElectrons,idxDLDJ,6,2); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,8,2);					
			}
			if (flavourSidebandCR) {
	 			doEleDistributionsPlots(leadingElectrons,idxDLDJ,6,3); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,8,3);					
			}

			doElePlots(leadingElectrons,idxDLDJ,4,0);
			if ( isInEB(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,4,1);
			if ( isInEE(leadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(leadingElectrons,idxDLDJ,4,2);	   
		}


		// --subLeadingEle
		if (multiLeptonMultiJets[idxDLDJ].isEEJJ) {
			nSubLeadingLeptons++;
			nEventsWithSubLeadingElePassingSelections++;

			if (signalRegion) {	
	 			doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,7,0); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,9,0);
			}
			if (lowMlljjCR) {
	 			doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,7,1); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,9,1);				
			}
			if (lowMllCR) {
	 			doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,7,2); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,9,2);					
			}
			if (flavourSidebandCR) {
	 			doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,7,3); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,9,3);					
			}

			doElePlots(subLeadingElectrons,idxDLDJ,5,0);
			if ( isInEB(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,5,1);
			if ( isInEE(subLeadingElectrons[idxDLDJ].v.Eta()) ) doElePlots(subLeadingElectrons,idxDLDJ,5,2); 
		}


		// --leadingMuon
		if (multiLeptonMultiJets[idxDLDJ].isMMJJ || multiLeptonMultiJets[idxDLDJ].isEMJJ) {
			nLeadingLeptons++;
			nEventsWithLeadingMuonPassingSelections++;

			if (signalRegion) {	
	 			doMuonDistributionsPlots(leadingMuons,idxDLDJ,6,0); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,10,0);
			}
			if (lowMlljjCR) {
	 			doMuonDistributionsPlots(leadingMuons,idxDLDJ,6,1); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,10,1);				
			}
			if (lowMllCR) {
	 			doMuonDistributionsPlots(leadingMuons,idxDLDJ,6,2); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,10,2);					
			}
			if (flavourSidebandCR) {
	 			doMuonDistributionsPlots(leadingMuons,idxDLDJ,6,3); 
				doMuonDistributionsPlots(leadingMuons,idxDLDJ,10,3);					
			}
			// devo usare i leadingMuonsRochCorr ??
		}				


		// --subLeadingMuon
		if (multiLeptonMultiJets[idxDLDJ].isMMJJ) {
			nSubLeadingLeptons++;
			nEventsWithSubLeadingMuonPassingSelections++;

			if (signalRegion) {	
	 			doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,7,0); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,11,0);
			}
			if (lowMlljjCR) {
	 			doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,7,1); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,11,1);				
			}
			if (lowMllCR) {
	 			doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,7,2); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,11,2);					
			}
			if (flavourSidebandCR) {
	 			doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,7,3); 
				doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,11,3);					
			}	
			// devo usare i subLeadingMuonsRochCorr	??
		}


		// --leadingJet & subLeadingJet
		nEventsWithLeadingJetPassingSelections++;
		nEventsWithSubLeadingJetPassingSelections++;

		if (signalRegion) {	
			doJetsDistributionsPlots(leadingJets,idxDLDJ,12,0);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,13,0);
		}
		if (lowMlljjCR) {
			doJetsDistributionsPlots(leadingJets,idxDLDJ,12,1);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,13,1);		
		}
		if (lowMllCR) {
			doJetsDistributionsPlots(leadingJets,idxDLDJ,12,2);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,13,2);				
		}
		if (flavourSidebandCR) {
			doJetsDistributionsPlots(leadingJets,idxDLDJ,12,3);		
			doJetsDistributionsPlots(subLeadingJets,idxDLDJ,13,3);				
		}


		nLeadingLeptons_histo -> Fill(nLeadingLeptons);
		nSubLeadingLeptons_histo -> Fill(nSubLeadingLeptons);


	} //fine loop su eventi

	
}
// ******************************************************************************************




// **************** 
void makePlots_multiLeptonMultiJet::saveHistosAndOutputFile(TString& ouputdir){
	// setTDRStyle();
	// writeExtraText = true; 
	// extraText = "Preliminary"; 
	// lumi_13TeV = "0.6 fb^{-1}"; 

	unsigned short numbOfHistograms = listOfHistograms.size();

	if (TnPEE || TnPMM) {
		TFile f(outputdir+suff+"/"+"Zmass_histos"+suff+".root","recreate");
		for (unsigned short l(0); l < numbOfHistograms; l++){
			string histoName = listOfHistograms[l]->GetName();
			if (histoName.find("Z") != string::npos) listOfHistograms[l]->Write(); 
		}

	} else {
		TFile f(outputdir+suff+"/"+"distributions_histos"+suff+".root","recreate");
		for (unsigned short l(0); l < numbOfHistograms; l++){
			string histoName = listOfHistograms[l]->GetName();
			if ( !(histoName.find("Z") != string::npos) ) listOfHistograms[l]->Write(); 

		// TCanvas* canv = new TCanvas();
		// gPad->SetRightMargin(0.05);
		// // CMS_lumi(canv, 4, 33);
		// listOfHistograms[l]->Draw();
		// canv->SaveAs(outputdir+suff+"/"+histoName+".png");
		// delete canv;
		// // }
		}
	}


	ofstream fOutput;
	if (TnPEE || TnPMM) fOutput.open(outputdir+suff+"/nEventsTnP"+suff+".dat");
	else fOutput.open(outputdir+suff+"/nEvents"+suff+".dat");
	fOutput << "nEvents = " << nEvents << " \n";
	fOutput << "nEventsPassingEEJJhlt = " << nEventsPassingEEJJhlt << " \n";
	fOutput << "nEventsPassingMMJJhlt = " << nEventsPassingMMJJhlt << " \n";
	fOutput << "nEventsPassingEMJJhlt = " << nEventsPassingEMJJhlt << " \n";
	fOutput << "nEventsPassingTandPEEhlt = " << nEventsPassingTandPEEhlt << " \n";
	fOutput << "nEventsPassingTandPMMhlt = " << nEventsPassingTandPMMhlt << " \n";
	fOutput << "nEventsPassingTrigger = " << nEventsPassingTrigger << " \n" << " \n";

	fOutput << "nElectrons = " << nOfElectrons << ", nGoodElectrons = " << nOfGoodElectrons << " \n" ;
	fOutput << "nMuons = " << nOfMuons << ", nGoodMuons = " << nOfGoodMuons << ", nGoodMuonsWithRochCor = " << nOfGoodMuonsWithRochCor << " \n" ;
	fOutput << "nJets = " << nOfJets << ", nGoodJets = " << nOfGoodJets << " \n" << " \n";

	fOutput << "nEventsWithAtLeast1DLDJ = " << nEventsWithAtLeast1DLDJ << " \n";
	fOutput << "nEventsWithRightLeptonPair = " << nEventsWithRightLeptonPair << " \n";
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

	fOutput << "nEventsWith2elePassingLoosePreselections = " << nEventsWith2elePassingLoosePreselections << " \n";
	fOutput << "nEventsWith2elePassingLoosePreselectionsAndCharge = " << nEventsWith2elePassingLoosePreselectionsAndCharge << " \n";
	fOutput << "nEventsWith2elePassingLoosePreselectionsAndChargeAndheepId = " << nEventsWith2elePassingLoosePreselectionsAndChargeAndHeepId << " \n";
	fOutput << "nEventsWith2muonsPassingLoosePreselections = " << nEventsWith2muonsPassingLoosePreselections << " \n";
	fOutput << "nEventsWith2muonsPassingLoosePreselectionsAndCharge = " << nEventsWith2muonsPassingLoosePreselectionsAndCharge << " \n";
	fOutput << "nEventsWith2muonsPassingLoosePreselectionsAndChargeAndHighPt = " << nEventsWith2muonsPassingLoosePreselectionsAndChargeAndHighPt << " \n" << " \n";


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



makePlots_multiLeptonMultiJet::makePlots_multiLeptonMultiJet(TString filename_, TString outputdir_, bool MC_, bool signalEE_, bool signalMuMu_, bool eMuSideband_, bool TnPEE_, bool TnPMM_):
	filename(filename_), outputdir(outputdir_), MC(MC_), signalEE(signalEE_), signalMuMu(signalMuMu_), eMuSideband(eMuSideband_), TnPEE(TnPEE_), TnPMM(TnPMM_)
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

	string inputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/miniTrees-Moriond17/";   
	string outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/Distributions-Moriond17/DoubleEG";

	bool signalEE = true;  //true for DoubleEG (if not TnP)   
	bool signalMuMu = false;  //true for SingleMuon (if not TnP)
	bool eMuSideband = false; //true for MuonEG

	bool TnPEE = false; //true for DoubleEG 
	bool TnPMM = false; //true for SingleMuon 

	if (signalMuMu || TnPMM) outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/Distributions-Moriond17/SingleMuon";
	if (eMuSideband) outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_26_patch1/src/dafne/Distributions-Moriond17/MuonEG";


	// makePlots_multiLeptonMultiJet(inputDir+"DYJetsPtBinned/output_DYJetsPtBinned_miniTree.root", outputDir+"/DYJetsPtBinned", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots_multiLeptonMultiJet(inputDir+"TTJets/output_TTJets_miniTree.root", outputDir+"/TTJets", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots_multiLeptonMultiJet(inputDir+"WJets/output_WJetsToLNu_miniTree.root", outputDir+"/WJets", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots_multiLeptonMultiJet(inputDir+"WZ/output_WZ_miniTree.root", outputDir+"/WZ", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots_multiLeptonMultiJet(inputDir+"ZZ/output_ZZ_miniTree.root", outputDir+"/ZZ", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots_multiLeptonMultiJet(inputDir+"WW/output_WW_miniTree.root", outputDir+"/WW", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots_multiLeptonMultiJet(inputDir+"SingleTop/output_SingleTop_miniTree.root", outputDir+"/SingleTop", true, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);


	makePlots_multiLeptonMultiJet(inputDir+"DoubleEG/output_DoubleEG_miniTree.root", outputDir, false, signalEE, false, false, TnPEE, false);

	// makePlots_multiLeptonMultiJet(inputDir+"SingleMuon/output_SingleMuon_miniTree.root", outputDir, false, false, signalMuMu, false, false, TnPMM);

	// makePlots_multiLeptonMultiJet(inputDir+"MuonEG/output_MuonEG_miniTree.root", outputDir, false, false, false, true, false, false);

}



#ifndef __CLING__
int main() {
  runMakePlots_multiLeptonMultiJet();
  return 0;
}
#endif