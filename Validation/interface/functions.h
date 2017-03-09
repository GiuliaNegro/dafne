#ifndef _functions_h_
#define _functions_h_
#include <iostream>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <TLorentzVector.h>

using namespace std;



bool isInEB(const float& eta){
	bool isInEB = fabs(eta) < 1.4442;
	return isInEB;
}


bool isInEE(const float& eta){
	bool isInEE = fabs(eta) > 1.566  &&  fabs(eta) < 2.5;
	return isInEE;
}



bool inEBEB(const float& eta1, const float& eta2){
	bool inEBEB = isInEB(eta1) && isInEB(eta2);
	return inEBEB;
}

bool inEEEE(const float& eta1, const float& eta2){
	bool inEEEE = isInEE(eta1) && isInEE(eta2);
	return inEEEE;
}

bool inEBEE(const float& eta1, const float& eta2){
	bool inEBEE = ((isInEB(eta1) && isInEE(eta2)) || (isInEE(eta1) && isInEB(eta2)));
	return inEBEE;
}



double phi0to2pi(double phi){
	double pi = 3.141592653589793238;
	while (phi >= 2.*pi) phi -= 2.*pi;
	while (phi < 0.) phi += 2.*pi;
	return phi;
}


double deltaPhi(TLorentzVector v1, TLorentzVector v2){
	// build the delta Phi angle between the two vectors
	// double pi = 3.141592653589793238;
	double dPhi = v1.Phi() - v2.Phi();
	return phi0to2pi(dPhi);
} 


double deltaR(TLorentzVector v1, TLorentzVector v2){
	double dEta = v1.Eta() - v2.Eta();
	double dPhi = deltaPhi(v1, v2);
	return sqrt(dEta * dEta + dPhi * dPhi);
}


float eff_sigma(std::vector<float> & v)
{
		size_t n = v.size();
		if (n < 2) return 0;
		std::sort(v.begin(), v.end());
		int s = floor(0.68269 * n);
		float d_min = v[s] - v[0];
		for (size_t i = s; i < n; ++i) {
				float d = v[i] - v[i - s];
				if (d < d_min) d_min = d;
		}
		return d_min / 2.;
}


#endif