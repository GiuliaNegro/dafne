import FWCore.ParameterSet.Config as cms
import json
import math
import os
from ROOT import TH2F, TFile


#Run2016B (5.929/fb)
#Run2016C (2.646/fb)
#Run2016D (4.353/fb)
#Run2016E (4.117/fb)
#Run2016F (3.186/fb)
#Run2016G (7.721/fb)
#Run2016H_v2 (8.636/fb)
#Run2016H_v3 (0.221/fb)

w1 = (5.929+2.646+4.353+4.117+3.186)/(5.929+2.646+4.353+4.117+3.186+7.721+8.636+0.221)
w2 = 1.0 - w1
print "w1 = ", w1
print "w2 = ", w2

Additional_ID_Systematics = 0.01    
Additional_ISO_Systematics = 0.01 
Additional_TR_Systematics = 0.01


f1 = open(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_IdBF.json'), 'r')
f2 = open(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_IdGH.json'), 'r')

g1 = open(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_IsoBF.json'), 'r')
g2 = open(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_IsoGH.json'), 'r')

h1 = open(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_TrigBF.json'), 'r')
h2 = open(os.path.join(os.environ['CMSSW_BASE'],'src/dafne/data/EfficienciesAndSF_TrigGH.json'), 'r')

ID_results1 = json.load(f1)
ID_results2 = json.load(f2)

Iso_results1 = json.load(g1)
Iso_results2 = json.load(g2)

Trig_results1 = json.load(h1)
Trig_results2 = json.load(h2)

SF_ID_C = [[0 for x in range(7)] for y in range(4)]   #x=pt, y=eta
SF_ID_E = [[0 for x in range(7)] for y in range(4)] 
SF_ID_C1 = [[0 for x in range(7)] for y in range(4)] 
SF_ID_E1 = [[0 for x in range(7)] for y in range(4)] 
SF_ID_C2 = [[0 for x in range(7)] for y in range(4)] 
SF_ID_E2 = [[0 for x in range(7)] for y in range(4)] 

SF_ISO_C = [[0 for x in range(7)] for y in range(4)] 
SF_ISO_E = [[0 for x in range(7)] for y in range(4)] 
SF_ISO_C1 = [[0 for x in range(7)] for y in range(4)] 
SF_ISO_E1 = [[0 for x in range(7)] for y in range(4)] 
SF_ISO_C2 = [[0 for x in range(7)] for y in range(4)] 
SF_ISO_E2 = [[0 for x in range(7)] for y in range(4)] 

SF_TR_C = [[0 for x in range(8)] for y in range(4)]  
SF_TR_E = [[0 for x in range(8)] for y in range(4)] 
SF_TR_C1 = [[0 for x in range(8)] for y in range(4)] 
SF_TR_E1 = [[0 for x in range(8)] for y in range(4)] 
SF_TR_C2 = [[0 for x in range(8)] for y in range(4)] 
SF_TR_E2 = [[0 for x in range(8)] for y in range(4)] 



print "ID"

i=0
for etaKey, values in sorted(ID_results1["MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta"]["abseta_pair_ne_ratio"].iteritems()) :
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		# print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
		SF_ID_C1[i][j] = result["value"]
		# SF_ID_E1[i][j] = math.sqrt(pow(result["error"],2) + pow(Additional_ID_Systematics,2))  
		SF_ID_E1[i][j] = result["error"]
		j += 1
	i += 1

i=0
for etaKey, values in sorted(ID_results2["MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta"]["abseta_pair_ne_ratio"].iteritems()) :
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		# print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
		SF_ID_C2[i][j] = result["value"]
		# SF_ID_E2[i][j] = math.sqrt(pow(result["error"],2) + pow(Additional_ID_Systematics,2))
		SF_ID_E2[i][j] = result["error"]
		j += 1
	i += 1	

i=0
for etaKey, values in sorted(ID_results2["MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta"]["abseta_pair_ne_ratio"].iteritems()) :
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		SF_ID_C[i][j] = SF_ID_C1[i][j]*w1 + SF_ID_C2[i][j]*w2
		SF_ID_E[i][j] = math.sqrt(pow(SF_ID_E1[i][j]*w1,2) + pow(SF_ID_E2[i][j]*w2,2))
		print etaKey, ptKey, "SF_ID_C = ", SF_ID_C[i][j], ", SF_ID_E = ", SF_ID_E[i][j]
		j += 1
	i += 1	



print "ISO"

i=0
for etaKey, values in sorted(Iso_results1["tkLooseISO_highptID_newpt_eta"]["abseta_pair_ne_ratio"].iteritems()) :
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		# print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
		SF_ISO_C1[i][j] = result["value"]
		# SF_ISO_E1[i][j] = math.sqrt(pow(result["error"],2) + pow(Additional_ISO_Systematics,2))
		SF_ISO_E1[i][j] = result["error"]
		j += 1
	i += 1

i=0
for etaKey, values in sorted(Iso_results2["tkLooseISO_highptID_newpt_eta"]["abseta_pair_ne_ratio"].iteritems()) :
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		# print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
		SF_ISO_C2[i][j] = result["value"]
		# SF_ISO_E2[i][j] = math.sqrt(pow(result["error"],2) + pow(Additional_ISO_Systematics,2))
		SF_ISO_E2[i][j] = result["error"]
		j += 1
	i += 1

i=0
for etaKey, values in sorted(Iso_results2["tkLooseISO_highptID_newpt_eta"]["abseta_pair_ne_ratio"].iteritems()) :
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		SF_ISO_C[i][j] = SF_ISO_C1[i][j]*w1 + SF_ISO_C2[i][j]*w2
		SF_ISO_E[i][j] = math.sqrt(pow(SF_ISO_E1[i][j]*w1,2) + pow(SF_ISO_E2[i][j]*w2,2))
		print etaKey, ptKey, "SF_ISO_C = ", SF_ISO_C[i][j], ", SF_ISO_E = ", SF_ISO_E[i][j]
		j += 1
	i += 1	



print "Trigger"

i=0
for etaKey, values in sorted(Trig_results1["Mu50_OR_TkMu50_PtEtaBins"]["abseta_pt_ratio"].iteritems()) :   
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		# print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
		SF_TR_C1[i][j] = result["value"]
		# SF_TR_E1[i][j] = math.sqrt(pow(result["error"],2) + pow(Additional_TR_Systematics,2))   
		SF_TR_E1[i][j] = result["error"]
		j += 1
	i += 1

i=0
for etaKey, values in sorted(Trig_results2["Mu50_OR_TkMu50_PtEtaBins"]["abseta_pt_ratio"].iteritems()) :   
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		# print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
		SF_TR_C2[i][j] = result["value"]
		# SF_TR_E2[i][j] = math.sqrt(pow(result["error"],2) + pow(Additional_TR_Systematics,2))     
		SF_TR_E2[i][j] = result["error"]
		j += 1
	i += 1

i=0
for etaKey, values in sorted(Trig_results2["Mu50_OR_TkMu50_PtEtaBins"]["abseta_pt_ratio"].iteritems()) :    
	j=0
	for ptKey, result in sorted(values.iteritems()) :
		SF_TR_C[i][j] = SF_TR_C1[i][j]*w1 + SF_TR_C2[i][j]*w2
		SF_TR_E[i][j] = math.sqrt(pow(SF_TR_E1[i][j]*w1,2) + pow(SF_TR_E2[i][j]*w2,2))
		print etaKey, ptKey, "SF_TR_C = ", SF_TR_C[i][j], ", SF_TR_E = ", SF_TR_E[i][j]
		j += 1
	i += 1	



