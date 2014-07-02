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
	
	float E;
	
	float pt;
    float pz;
    float phi;
    float eta;
    float y;
	
	int jet_id;
	int mother_id;
};

struct jet_t {
	int NParticles;
	
	double E;
	double E_T;
	double axis_phi;
	double axis_eta;
	double e_E;
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

float deltaPhi(float phi1, float phi2) {
	float deltaphi = phi1 - phi2;
	if (deltaphi > PI) {
		deltaphi -= 2.0*PI;
	}
	else if (deltaphi < -1.0*PI) {
		deltaphi += 2.0*PI;
	}
	
	return deltaphi;
	
}

void fixed_cone_e_jet(jet_t& jet, particle_t event_electron, vector<particle_t> event_assoparticles, double conesize) {
	
	jet.NParticles = 0;
	jet.E = 0.;
	jet.E_T = 0.;
	jet.axis_phi = 0.;
	jet.axis_eta = 0.;
	jet.e_E = 0.;
	
	double jet_sum = 0.0;
	double e_pt = 0.0;
	double htot_E = 0.0;
	double htot_pt = 0.0;
	
	double Ephi = 0.0;
	double Eeta = 0.0;
	
	if (event_electron.eta < .7 && event_electron.eta > -.7 && event_electron.pt > 3.0 && event_electron.pt < 10.0) {
		e_pt = event_electron.pt;
		for (vector<particle_t>::iterator asso_iter = event_assoparticles.begin(); asso_iter != event_assoparticles.end(); asso_iter++) {
			if (asso_iter->eta < 1.0 && asso_iter->eta > -1.0) {
				if (asso_iter->pt < .2) continue;
				double dEta = asso_iter->eta - event_electron.eta;
				double dPhi = deltaPhi(asso_iter->phi, event_electron.phi);
				if (sqrt(dEta*dEta + dPhi*dPhi) < conesize) {
					htot_pt += asso_iter->pt;
					htot_E += asso_iter->E;
					//cout << "asso_iter->pt = " << asso_iter->pt << endl;
					//cout << "asso_iter->E = " << asso_iter->E << endl;
					jet.NParticles++;
					Ephi += deltaPhi(event_electron.phi , asso_iter->phi)*asso_iter->pt; //phi coordinate defined relative to e
					Eeta += asso_iter->eta*asso_iter->pt;
					}
				}
			}
		}
	jet.NParticles++;
	htot_E += event_electron.E;
	htot_pt += event_electron.pt;
	Eeta += event_electron.eta*event_electron.pt;
	//cout << "htot_E = " << htot_E << endl;
	//cout << "htot_pt = " << htot_pt << endl;
	
	//calculate phi coordinate
	double phi_of_jet = event_electron.phi - Ephi/htot_pt;
	if (phi_of_jet < 0.0) {
		phi_of_jet += 2.0*PI;
	} 
	else if (phi_of_jet > 2.0*PI) {
		phi_of_jet -= 2.0*PI;
	}

	
	//make jet
	jet.E = htot_E;
	jet.E_T = htot_pt;
	jet.axis_phi = phi_of_jet;
	jet.axis_eta = Eeta/htot_pt;
	jet.e_E = e_pt;
}

void Pythia_NPE_jet() {
	
	string inputfiles;
	cout << "Enter Input Files: " << endl;
	cin >> inputfiles;
	
	TChain *chain = new TChain("tree");
	int nfile = 0;
	nfile += chain->Add(inputfiles.data());
	cout << "Added " << nfile << " to chain." << endl;
	
	string ofname = "/Users/jaydunkelberger/PYTHIA/pythiasimulations/e_jet/data/pythia_NPE_jet_output.root";
	TFile *outfile = new TFile(ofname.data(), "RECREATE");
	
	//Histograms
	TH1D *hJetSum_conesize[10];
	TH1D *hJet_ET_conesize[10];
	TH1D *hNPEFraction_conesize[10];
	TH1D *hJetSum_ptcut[10];
	TH1D *hNPEFraction_ptcut[10];
	TH1D *hNParticles_conesize[10];
	
	for (int i=0; i<5; i++) {
		hJetSum_conesize[i] = new TH1D(Form("hJetSum_conesize_%d", i), Form("Sum of energy in jet for cone size %d", i), 400, 0.0, 40.0);
		hJet_ET_conesize[i] = new TH1D(Form("hJet_ET_conesize_%d", i), Form("Sum of E_T in jet for cone size %d", i), 400, 0.0, 40.0);
		hNPEFraction_conesize[i] = new TH1D(Form("hNPEFraction_conesize_%d", i), Form("NPE/jettotal for cone size %d", i), 100, 0.0, 1.0);
		
		//hJetSum_ptcut[i] = new TH1D(Form("hJetSum_ptcut_%d", i), Form("Sum of energy in jet for pt cut %d", i), 400, 0.0, 40.0);
		//hNPEFraction_ptcut[i] = new TH1D(Form("hNPEFraction_ptcut_%d", i), Form("NPE/jettotal for pt cut %d", i), 100, 0.0, 1.0);
		
		hNParticles_conesize[i] = new TH1D(Form("hNParticles_conesize_%d", i), Form("NParticles in jet for conesize %d", i), 20, 0, 20);
	}
	
	TH1D *hjet_phi = new TH1D("hjet_phi", "phi centroid of jet", 100, 0, 2.0*PI);
	TH1D *hjet_eta = new TH1D("hjet_eta", "eta centroid of jet", 100, -1.0, 1.0);
	
	particle_t current_particle;
	vector<particle_t> event_electrons;
	vector<particle_t> event_assoparticles;
	
	chain->SetBranchAddress("asso_particle", &current_particle);
	
	int nentries =  chain->GetEntries();
	int prev_event = -1;
	bool first_event = true;
	int NUpsilons = 0;
	
	jet_t current_jet;
	
	for (int ichain=0; ichain<nentries; ichain++) {
		if (ichain % 10000 == 0) cout << "Analyzing Entry: " << ichain << endl;
		
		chain->GetEntry(ichain);
		
		int current_event = current_particle.event_no;
		
		if (current_event != prev_event && !first_event) { //end of event, build correlation
			
			for (vector<particle_t>::iterator e_iter = event_electrons.begin(); e_iter != event_electrons.end(); e_iter++) {
	
				if (e_iter->pt < 3 || e_iter->pt > 10 || e_iter->eta > .7 || e_iter->eta < -.7) continue; 
				fixed_cone_e_jet(current_jet, (*e_iter), event_assoparticles, 0.3);

				//cout << "current_jet.E = " << current_jet.E << endl;
				//cout << "current_jet.E_T = " << current_jet.E_T << endl;
				//cout << "current_jet.e_E = " << current_jet.e_E << endl;

				if(current_jet.E > 0.0) { //valid jet fill histograms
					for (int i=0; i<5; i++) {
						hJetSum_conesize[i]->Fill(current_jet.E);
						hJet_ET_conesize[i]->Fill(current_jet.E_T);
						hNPEFraction_conesize[i]->Fill(current_jet.e_E/current_jet.E);
						hNParticles_conesize[i]->Fill(current_jet.NParticles);
					}
					hjet_phi->Fill(current_jet.axis_phi);
					hjet_eta->Fill(current_jet.axis_eta);
				}
			}
			
			for (vector<particle_t>::iterator e_iter = event_electrons.begin(); e_iter != event_electrons.end(); e_iter++) {
				//check electron pt/acceptance
				double jet_sum = 0.0;
				double e_pt = 0.0;
				double htot_E[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
				double htot_E_ptcut[5] =  {0.0, 0.0, 0.0, 0.0, 0.0};
				/*
				if (e_iter->eta < .7 && e_iter->eta > -.7 && e_iter->pt > 3.0 && e_iter->pt < 10.0) {
					e_pt = e_iter->pt;
					for (vector<particle_t>::iterator asso_iter = event_assoparticles.begin(); asso_iter != event_assoparticles.end(); asso_iter++) {
						if (asso_iter->eta < 1.0 && asso_iter->eta > -1.0) {
							if (asso_iter->pt < .2) continue;
							double dEta = asso_iter->eta - e_iter->eta;
							double dPhi = asso_iter->phi - e_iter->phi;
							double minimum_jet_cone_size = 0.2;
							for (int jet_cone_inc_number = 0; jet_cone_inc_number<5; jet_cone_inc_number++) {
								if (sqrt(dEta*dEta + dPhi*dPhi) < (minimum_jet_cone_size + .1*jet_cone_inc_number)) {
									htot_E[jet_cone_inc_number] += asso_iter->pt;
								}
							}
							for (int pt_low_cutoff_index = 0; pt_low_cutoff_index < 5; pt_low_cutoff_index++) {
								if (asso_iter->pt > (.2 + .1*pt_low_cutoff_index)) {
									htot_E_ptcut[pt_low_cutoff_index] += asso_iter->pt; 
								}
							}
						}
					}
					
				}
				 */
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
