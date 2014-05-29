/*                          *
 *  Pythia_eh_Correlation.C * 
 *                          */
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#define PI 3.141592654

struct particle_t
{
	int event_no;
	
	int   id;
	
	float pt;
    float pz;
    float phi;
    float eta;
    float y;
	
	int mother_id;
};

//calculate pt bins
int calc_asso_pt_bin(double asso_pt) {
	int asso_pt_bin = -1;
	if(asso_pt > 0.2 && asso_pt < 0.3) asso_pt_bin = 0;
	else if(asso_pt > 0.3 && asso_pt < 0.4) asso_pt_bin = 1;
	else if(asso_pt > 0.4 && asso_pt < 0.5) asso_pt_bin = 2;
	else if(asso_pt > 0.5 && asso_pt < 0.6) asso_pt_bin = 3;
	else if(asso_pt > 0.6 && asso_pt < 0.8) asso_pt_bin = 4;
	else if(asso_pt > 0.8 && asso_pt < 1.0) asso_pt_bin = 5;
	else if(asso_pt > 1.0 && asso_pt < 1.5) asso_pt_bin = 6;
	else if(asso_pt > 1.5 && asso_pt < 2.0) asso_pt_bin = 7;
	else if(asso_pt > 2.0 && asso_pt < 2.5) asso_pt_bin = 8;
	else if(asso_pt > 2.5) asso_pt_bin = 9;
	
	return asso_pt_bin;
}

int calc_prim_pt_bin(double prim_pt) {
	
	int prim_pt_bin = -1;
	if(prim_pt<=1.0) prim_pt_bin = 0;
	else if(prim_pt>=1.0&&prim_pt<=2.0) prim_pt_bin = 1;
	else if(prim_pt>=2.0&&prim_pt<=3.0) prim_pt_bin = 2;
	else if(prim_pt>=3.0&&prim_pt<=4.0) prim_pt_bin = 3;
	else if(prim_pt>=4.0&&prim_pt<=5.0) prim_pt_bin = 4;
	else if(prim_pt>=5.0&&prim_pt<=6.0) prim_pt_bin = 5;
	else if(prim_pt>=6.0&&prim_pt<=7.0) prim_pt_bin = 6;
	else if(prim_pt>=7.0&&prim_pt<=8.0) prim_pt_bin = 7;
	else if(prim_pt>=8.0&&prim_pt<=9.0) prim_pt_bin = 8;
	else if(prim_pt>=9.0) prim_pt_bin = 9;
	
	return prim_pt_bin;
	
}


bool isInAcceptance(particle_t decay) {
	return (decay.eta < 1.0 && decay.eta > -1.0);
}

float dPhi(float phi1, float phi2) {
	float deltaphi = phi1 - phi2;
	if(deltaphi < -PI/2.0) deltaphi += 2.0*PI;
	if(deltaphi > 3.0*PI/2.0) deltaphi -= 2.0*PI;
	
	return deltaphi;
	
}

void Pythia_NPE_jet() {
	
	string inputfiles;
	cout << "Enter Input Files: " << endl;
	cin >> inputfiles;
	
	TChain *chain = new TChain("tree");
	int nfile = 0;
	nfile += chain->Add(inputfiles.data());
	cout << "Added " << nfile << " to chain." << endl;
	
	string ofname = "/Users/jaydunkelberger/PYTHIA/pythiasimulations/eh_corr/data/pythia_NPE_jet_output.root";
	TFile *outfile = new TFile(ofname.data(), "RECREATE");
	
	//Histograms
	TH1D *hJetSum[10];
	TH1D *hNPEFraction[10];
	
	for (int i=0; i<10; i++) {
		hJetSum[i] = new TH1D(Form("hJetSum_%d", i), Form("Sum of energy in jet for %d", i), 400, 0.0, 40.0);
		hNPEFraction[i] = new TH1D(Form("hNPEFraction_%d", i), Form("NPE/jettotal %d", i), 100, 0.0, 1.0);
	}
	
	particle_t current_particle;
	vector<particle_t> event_electrons;
	vector<particle_t> event_assoparticles;
	
	chain->SetBranchAddress("asso_particle", &current_particle);
	
	int nentries =  chain->GetEntries();
	int prev_event = -1;
	bool first_event = true;
	int NUpsilons = 0;
	
	for (int ichain=0; ichain<nentries; ichain++) {
		if (ichain % 10000 == 0) cout << "Analyzing Entry: " << ichain << endl;
		
		chain->GetEntry(ichain);
		
		int current_event = current_particle.event_no;
		
		if (current_event != prev_event && !first_event) { //end of event, build correlation
			for (vector<particle_t>::iterator e_iter = event_electrons.begin(); e_iter != event_electrons.end(); e_iter++) {
				//check electron pt/acceptance
				double jet_sum = 0.0;
				double e_E = 0.0;
				double htot_E = 0.0;
				if (e_iter->eta < .7 && e_iter->eta > -.7 && e_iter->pt > 3.0 && e_iter->pt < 10.0) {
					e_E = sqrt(e_iter->pt*e_iter->pt + e_iter->pz*e_iter->pz);
					for (vector<particle_t>::iterator asso_iter = event_assoparticles.begin(); asso_iter != event_assoparticles.end(); asso_iter++) {
						if (asso_iter->eta < 1.0 && asso_iter->eta > -1.0) {
							if (asso_iter->pt < .2) continue;
							double dEta = asso_iter->eta - e_iter->eta;
							double dPhi = asso_iter->phi - e_iter->phi;
							if (sqrt(dEta*dEta + dPhi*dPhi) < .2) { //jet cone condition
								htot_E += sqrt(asso_iter->pt*asso_iter->pt + asso_iter->pz*asso_iter->pz); //just using as a proxy for E now will change later
							}
						}
					}
					//at the end of each seed e sum up the jet
					for (int i=0; i<10; i++) {
						hJetSum[i]->Fill(htot_E + e_E);
						hNPEFraction[i]->Fill(e_E/(htot_E + e_E));
					}
				}
			}
			event_electrons.clear();
			event_assoparticles.clear();
		}
		else {
			if (current_particle.id == fabs(11)) { //is an electron
				event_electrons.push_back(current_particle);
			}
			else {
				event_assoparticles.push_back(current_particle);
			}
			first_event = false;
		}
		
		
		prev_event = current_event;
		
	}
	
	outfile->Write();
	return;
}