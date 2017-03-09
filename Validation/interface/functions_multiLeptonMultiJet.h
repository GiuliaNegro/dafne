#ifndef _functions_multiLeptonMultiJet_h_
#define _functions_multiLeptonMultiJet_h_
#include <iostream>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <TLorentzVector.h>

using namespace std;



struct eleStruct{
	eleStruct();
	eleStruct(float pt_, float eta_, float phi_, float en_, bool passHEEPId_, unsigned int HEEPBitMapValues_, 
		bool passMediumId_, float iso_, float dzvtx_, float d0vtx_, int isMatchedToGen_, int charge_, 
		float etaSC_, bool isEcalDriven_, float dEtaIn_, float dPhiIn_, float HoE_, float r9_, 
		float sigmaIetaIeta_, float e5x5_, float e1x5_, float e2x5_, float EmHadDepth1Iso_, float ptTracksIso_, 
		int missingHits_, float dxy_, float eOverP_, float ecalEnergy_, float hcalOverEcal_) {
		v.SetPtEtaPhiE(pt_, eta_, phi_, en_);
		passHEEPId = passHEEPId_;
		HEEPBitMapValues = HEEPBitMapValues_;
		passMediumId = passMediumId_;
		iso = iso_;
		dzvtx = dzvtx_;
		d0vtx = d0vtx_;
		isMatchedToGen = isMatchedToGen_;
		charge = charge_;
		etaSC = etaSC_;
		isEcalDriven = isEcalDriven_;
		dEtaIn = dEtaIn_;
		dPhiIn = dPhiIn_;
		HoE = HoE_;		
		r9 = r9_;
		sigmaIetaIeta = sigmaIetaIeta_;
		e5x5 = e5x5_;
		e1x5 = e1x5_;
		e2x5 = e2x5_;
		EmHadDepth1Iso = EmHadDepth1Iso_;
		ptTracksIso = ptTracksIso_;
		missingHits = missingHits_;
		dxy = dxy_;
		eOverP = eOverP_;
		ecalEnergy = ecalEnergy_;
		hcalOverEcal = hcalOverEcal_;
	}
	TLorentzVector v;
	float iso, dzvtx, d0vtx, etaSC, dEtaIn, dPhiIn, HoE, r9, sigmaIetaIeta, e5x5, e1x5, e2x5, EmHadDepth1Iso, 
	ptTracksIso, dxy, eOverP, ecalEnergy, hcalOverEcal; 
	bool passHEEPId, passMediumId, isEcalDriven;
	int isMatchedToGen, charge, missingHits;
	unsigned int HEEPBitMapValues;
};


struct muonStruct{
	muonStruct();
	muonStruct(float pt_, float eta_, float phi_, float en_, float iso_, bool isHighPt_, int isMatchedToGen_, 
		int charge_, float dzvtx_, float dxy_, float rochCor_) {
		v.SetPtEtaPhiE(pt_, eta_, phi_, en_);
		iso = iso_;
		isHighPt = isHighPt_;
		isMatchedToGen = isMatchedToGen_;
		charge = charge_;
		dzvtx = dzvtx_;
		dxy = dxy_;
		rochCor = rochCor_;
	}
	TLorentzVector v;
	float iso, dzvtx, dxy, rochCor; 
	bool isHighPt;
	int isMatchedToGen, charge;
};


struct jetStruct{
	jetStruct();
	jetStruct(double pt_, double eta_, double phi_, double en_, int isMatchedToGen_, bool isTight_) {
		v.SetPtEtaPhiE(pt_, eta_, phi_,en_);
		isMatchedToGen = isMatchedToGen_;
		isTight = isTight_;
	}
	TLorentzVector v;
	int isMatchedToGen;
	bool isTight;
};



struct multiLeptonMultiJetStruct{
	multiLeptonMultiJetStruct();
	multiLeptonMultiJetStruct(bool isEEJJ_, bool isEETT_, bool isMMJJ_, bool isMMTT_, bool isEMJJ_, 
		bool isSignalRegion_, bool isLowMllCR_, bool isLowMlljjCR_, bool isBB_, bool isEE_, bool isEB_, 
		bool passPreselections_, float dRLeadLeptonLeadJet_, float dRLeadLeptonSubLeadJet_, 
		float dRSubLeadLeptonLeadJet_, float dRSubLeadLeptonSubLeadJet_, float multiLeptonMultiJet_sumPt_, 
		float multiLeptonMultiJet_invMass_, float diLepton_invMass_, float diLepton_pt_, float diJet_invMass_, 
		float diJetLeadingLepton_invMass_, float diJetSubLeadingLepton_invMass_){
		isEEJJ = isEEJJ_;
		isEETT = isEETT_;
		isMMJJ = isMMJJ_;
		isMMTT = isMMTT_;
		isEMJJ = isEMJJ_;
		isSignalRegion = isSignalRegion_;
		isLowMllCR = isLowMllCR_;
		isLowMlljjCR = isLowMlljjCR_;
		isBB = isBB_;
		isEE = isEE_;
		isEB = isEB_;
		passPreselections = passPreselections_;
		dRLeadLeptonLeadJet = dRLeadLeptonLeadJet_;
		dRLeadLeptonSubLeadJet = dRLeadLeptonSubLeadJet_;
		dRSubLeadLeptonLeadJet = dRSubLeadLeptonLeadJet_;
		dRSubLeadLeptonSubLeadJet = dRSubLeadLeptonSubLeadJet_;
		multiLeptonMultiJet_sumPt = multiLeptonMultiJet_sumPt_;
		multiLeptonMultiJet_invMass = multiLeptonMultiJet_invMass_;
		diLepton_invMass = diLepton_invMass_;
		diLepton_pt = diLepton_pt_;
		diJet_invMass = diJet_invMass_;
		diJetLeadingLepton_invMass = diJetLeadingLepton_invMass_;
		diJetSubLeadingLepton_invMass = diJetSubLeadingLepton_invMass_;		
	}
	bool isEEJJ, isEETT, isMMJJ, isMMTT, isEMJJ, isSignalRegion, isLowMllCR, isLowMlljjCR, isBB, isEE, isEB, 
	passPreselections; 
	float dRLeadLeptonLeadJet, dRLeadLeptonSubLeadJet, dRSubLeadLeptonLeadJet, dRSubLeadLeptonSubLeadJet, 
	multiLeptonMultiJet_sumPt, multiLeptonMultiJet_invMass, diLepton_invMass, diLepton_pt, diJet_invMass, 
	diJetLeadingLepton_invMass, diJetSubLeadingLepton_invMass;
};



#endif