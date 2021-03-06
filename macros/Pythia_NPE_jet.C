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

struct final_state_particle_t
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
	
	float mother_E;
	float mother_pt;
	float mother_pz;
    float mother_phi;
    float mother_eta;
    float mother_y;
	
	int HFQ_id;
	
	float HFQ_E;
	float HFQ_pt;
	float HFQ_pz;
    float HFQ_phi;
    float HFQ_eta;
    float HFQ_y;
};

struct jet_t {
	int NParticles;
	
	float E;
	float E_T;
	float axis_phi;
	float axis_eta;
	float e_E;
};

//calculate pt bins
int calc_asso_pt_bin(float asso_pt) {
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

int calc_prim_pt_bin(float prim_pt) {
	
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


bool isInAcceptance(final_state_particle_t decay) {
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

float rotatePhi(float angle, float rotation) {
	float new_angle = angle + rotation;
	if (new_angle > PI) {
		new_angle -= 2.0*PI;
	}
	else if (new_angle < -PI) {
		new_angle += 2.0*PI;
	}
	return new_angle;
}

void fixed_cone_e_jet(jet_t& jet, final_state_particle_t event_electron, vector<final_state_particle_t> event_assoparticles, float conesize) {
	
	jet.NParticles = 0;
	jet.E = 0.;
	jet.E_T = 0.;
	jet.axis_phi = 0.;
	jet.axis_eta = 0.;
	jet.e_E = 0.;
	
	float jet_sum = 0.0;
	float e_pt = 0.0;
	float htot_E = 0.0;
	float htot_pt = 0.0;
	
	float Ephi = 0.0;
	float Eeta = 0.0;
	
	e_pt = event_electron.pt;
	for (vector<final_state_particle_t>::iterator asso_iter = event_assoparticles.begin(); asso_iter != event_assoparticles.end(); asso_iter++) {
		if (asso_iter->eta < 1.0 && asso_iter->eta > -1.0) {
			if (asso_iter->pt < .2) continue;
			float dEta = asso_iter->eta - event_electron.eta;
			float dPhi = deltaPhi(asso_iter->phi, event_electron.phi);
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

	jet.NParticles++;
	htot_E += event_electron.E;
	htot_pt += event_electron.pt;
	Eeta += event_electron.eta*event_electron.pt;
	//cout << "htot_E = " << htot_E << endl;
	//cout << "htot_pt = " << htot_pt << endl;
	
	//calculate phi coordinate
	float phi_of_jet = event_electron.phi - Ephi/htot_pt;
	if (phi_of_jet < -PI) {
		phi_of_jet += 2.0*PI;
	} 
	else if (phi_of_jet > PI) {
		phi_of_jet -= 2.0*PI;
	}

	
	//make jet
	jet.E = htot_E;
	jet.E_T = htot_pt;
	jet.axis_phi = phi_of_jet;
	jet.axis_eta = Eeta/htot_pt;
	jet.e_E = e_pt;
}
/*
void kt_e_jet(jet_t& jet, final_state_particle_t event_electron, vector<final_state_particle_t> event_assoparticles, float power) {
	
	jet.NParticles = 0;
	jet.E = 0.;
	jet.E_T = 0.;
	jet.axis_phi = 0.;
	jet.axis_eta = 0.;
	jet.e_E = 0.;
	
	float jet_sum = 0.0;
	float e_pt = 0.0;
	float htot_E = 0.0;
	float htot_pt = 0.0;
	
	float Ephi = 0.0;
	float Eeta = 0.0;
	
	e_pt = event_electron.pt;
	for (vector<final_state_particle_t>::iterator asso_iter = event_assoparticles.begin(); asso_iter != event_assoparticles.end(); asso_iter++) {
		if (asso_iter->eta < 1.0 && asso_iter->eta > -1.0) {
			if (asso_iter->pt < .2) continue;
			float dEta = asso_iter->eta - event_electron.eta;
			float dPhi = deltaPhi(asso_iter->phi, event_electron.phi);
			float kt1_power = pow(asso_iter->pt, 2.0*power);
			float kt2_power = pow(e_pt, 2.0*power);
			float kt_power = (kt1_power > kt2_power) ? kt2_power : kt1_power;
			
			//cout << "kt1_power: " << kt1_power << endl;
			//cout << "kt2_power: " << kt2_power << endl;
			//cout << "kt_power: " << kt_power << endl;
			
			float distance = kt_power*sqrt(dEta*dEta + dPhi*dPhi);
			
			float dBeam = pow(asso_iter->pt, 2.0*power);
			
			if (distance < dBeam) {
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
	jet.NParticles++;
	htot_E += event_electron.E;
	htot_pt += event_electron.pt;
	Eeta += event_electron.eta*event_electron.pt;
	//cout << "htot_E = " << htot_E << endl;
	//cout << "htot_pt = " << htot_pt << endl;
	
	//calculate phi coordinate
	float phi_of_jet = event_electron.phi - Ephi/htot_pt;
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
*/
void Pythia_NPE_jet() {
	
	string inputfiles;
	cout << "Enter Input Files: " << endl;
	cin >> inputfiles;
	
	TChain *chain = new TChain("tree");
	int nfile = 0;
	nfile += chain->Add(inputfiles.data());
	cout << "Added " << nfile << " to chain." << endl;
	
	string ofname = "/Users/jaydunkelberger/PYTHIA/pythiasimulations/e_jet/data/pythia_NPE_awayjet_weight_25M_tight_eta_output.root";
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
	
	TH1D *hjet_phi = new TH1D("hjet_phi", "phi centroid of jet", 100, -PI, PI);
	TH1D *hjet_eta = new TH1D("hjet_eta", "eta centroid of jet", 100, -1.0, 1.0);
	
	TH1D *hawayside_Etot = new TH1D("hawayside_Etot", "Total energy in away side cone", 400, 0.0, 40.0);
	TH1D *hdPhi_mother = new TH1D("hdPhi_mother", "dPhi of electron to parent particle", 80, -PI, PI);
	TH1D *hdPhi_HFQ = new TH1D("hdPhi_HFQ", "dPhi of electron to initial heavy flavor quark", 80, -PI, PI);
	
	TH2D *hJet_HFQ_Ecorr = new TH2D("hJet_HFQ_Ecorr", "Correlation of jet energy to initial HFQ energy", 400, 0.0, 40.0, 400, 0.0, 40.0);
	TH2D *hHFQ_charm_Ecorr = new TH2D("hHFQ_charm_Ecorr", "Correlation of jet energy to initial charm HFQ energy", 400, 0.0, 40.0, 400, 0.0, 40.0);
	TH2D *hHFQ_bottom_Ecorr = new TH2D("hHFQ_bottom_Ecorr", "Correlation of jet energy to initial bottom HFQ energy", 400, 0.0, 40.0, 400, 0.0, 40.0);
	TH1D *hHFQ_Pt_all_jets = new TH1D("hHFQ_Pt_all_jets", "HFQ Pt for high energy jet", 400, 0.0, 40.0);
	TH1D *hHFQ_Pt_high_away_side = new TH1D("hHFQ_Pt_high_away_side", "HFQ Pt for back to back high energy jets", 400, 0.0, 40.0);
	TH1D *hHFQ_Pt_ee_pair = new TH1D("hHFQ_Pt_ee_pair", "HFQ Pt for pair of high pt electrons", 400, 0.0, 40.0);
	TH1D *hHFQ_weighted_Pt_ee_pair = new TH1D("hHFQ_weighted_Pt_ee_pair", "Weighted HFQ Pt for pair of high pt electrons", 400, 0.0, 40.0);

	final_state_particle_t current_particle;
	vector<final_state_particle_t> event_electrons;
	vector<final_state_particle_t> event_assoparticles;
	
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
	
			//one special loop for 2 electron events
			if (event_electrons.size() == 2) {
				bool valid_prev_NPE = false;
				final_state_particle_t prev_NPE;
				for (vector<final_state_particle_t>::iterator e_iter = event_electrons.begin(); e_iter != event_electrons.end(); e_iter++) {
					if (e_iter->pt < 3 || e_iter->pt > 10 || e_iter->eta > .7 || e_iter->eta < -.7) continue;
					
					if (valid_prev_NPE) {
						double HFQ_PT = e_iter->HFQ_pt;
						if (HFQ_PT < .1) continue;
						hHFQ_Pt_ee_pair->Fill(HFQ_PT);	
						hHFQ_weighted_Pt_ee_pair->Fill(HFQ_PT, 1.0/(HFQ_PT*HFQ_PT*HFQ_PT*HFQ_PT*HFQ_PT*HFQ_PT));
						/*if (fabs(e_iter->HFQ_id) == 4) {
							hHFQ_charm_Ecorr->Fill(e_iter->HFQ_pt, e_iter->pt);	
						}
						else if (fabs(e_iter->HFQ_id) == 5) {
							hHFQ_bottom_Ecorr->Fill(e_iter->HFQ_pt, e_iter->pt);	
						}*/
					}
					valid_prev_NPE = true;
				}
			}
			
			for (vector<final_state_particle_t>::iterator e_iter = event_electrons.begin(); e_iter != event_electrons.end(); e_iter++) {
	
				if (e_iter->eta < .7 && e_iter->eta > -.7) {
					if (fabs(e_iter->HFQ_id) == 4) {
						hHFQ_charm_Ecorr->Fill(e_iter->HFQ_pt, e_iter->pt);	
					}
					else if (fabs(e_iter->HFQ_id) == 5) {
						hHFQ_bottom_Ecorr->Fill(e_iter->HFQ_pt, e_iter->pt);	
					}
				}
					
				if (e_iter->pt < 3 || e_iter->pt > 10 || e_iter->eta > .2 || e_iter->eta < -.2) continue; 
				
				hdPhi_mother->Fill(deltaPhi(e_iter->phi, e_iter->mother_phi));
				hdPhi_HFQ->Fill(deltaPhi(e_iter->phi, e_iter->HFQ_phi));
				
				fixed_cone_e_jet(current_jet, (*e_iter), event_assoparticles, 0.3);
				//kt_e_jet(current_jet, (*e_iter), event_assoparticles, 1.0);
				//cout << "current_jet.E = " << current_jet.E << endl;
				//cout << "current_jet.E_T = " << current_jet.E_T << endl;
				//cout << "current_jet.e_E = " << current_jet.e_E << endl;

				if(current_jet.E > 0.0/* && current_jet.NParticles >= 2*/) { //valid jet fill histograms, at least one clustered particle
					for (int i=0; i<5; i++) {
						hJetSum_conesize[i]->Fill(current_jet.E);
						hJet_ET_conesize[i]->Fill(current_jet.E_T);
						hNPEFraction_conesize[i]->Fill(current_jet.e_E/current_jet.E);
						hNParticles_conesize[i]->Fill(current_jet.NParticles);
						hJet_HFQ_Ecorr->Fill(e_iter->HFQ_pt, current_jet.E_T);
					}
					hjet_phi->Fill(current_jet.axis_phi);
					hjet_eta->Fill(current_jet.axis_eta);
					
					//loop through all FS particles and sum up those in away side jet
					float away_jet_axis_phi = rotatePhi(current_jet.axis_phi, -PI);
					float away_jet_axis_eta = -current_jet.axis_eta;
					float away_jet_Etot = 0.0;
					for (vector<final_state_particle_t>::iterator fs_iter = event_assoparticles.begin(); fs_iter != event_assoparticles.end(); fs_iter++) {
						if (fs_iter->eta < 1.0 && fs_iter->eta > -1.0 && fs_iter->pt > .2) { //acceptance cuts
							float dEta = fs_iter->eta - away_jet_axis_eta;
							float dPhi = deltaPhi(fs_iter->phi, away_jet_axis_phi);
							if (fabs(dEta) < 1.0 && fabs(dPhi) < 1.0) { //within the away side
								away_jet_Etot += fs_iter->pt;
							}
						}
					}
	
					hawayside_Etot->Fill(away_jet_Etot);
					hHFQ_Pt_all_jets->Fill(e_iter->HFQ_pt);
					if (away_jet_Etot > .6*current_jet.E_T) hHFQ_Pt_high_away_side->Fill(e_iter->HFQ_pt); 
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
