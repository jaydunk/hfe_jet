//=============================================================================
// NOTES:
// 1:tau does decay to electron(10^-15s), but not counted in the branch ratio, thus not recorded
// muon decays to electrons, but meantime is (10^-6s)
//=============================================================================
// V1.0 of eK_BJpsi: study e-kaon and B-jpsi correlations
// adopted from pmainHF2e, see the statement belo.
// largely modified to simplife it, e.g. no Ntuples being made.
// Jan. 2011 Wenqin Xu
//
// V1.1
//         if(id>9900000) id-=9900000; //deal with directly gluon pair
//         added to Flavor determination
// Jan. 2011 Wenqin 
//=============================================================================
// Statement of pmainHF2e
// V1.0 adopted from Thomas directly. See his statement below 
// Feb. 2010 Wenqin Xu
// 
// V2.0 added seeking decaying chain 
// and changed Ntuple structure from "an event per entry" to "an electron per entry" 
// Mar 9, 2010 Wenqin Xu
// 
// V2.1 added seeking B->bottom chain and histograms for the y of this chain
// and changed the flavor calculation to deal with rear excited states, e.g 10521
// Mar 15, 2010 Wenqin Xu
//
// V2.2 added seeking recoiled bottom-> original bottom produced by initial interaction. 
// The criteria is that non-recoil bottom should not have bottom mother
// Mar 16, 2010 Wenqin Xu
//
// V2.3 for the D->e branch, changed the orig to be Bottom if D is from B after generations
// The old way didn't deal with B->intermediate generations->D->e case.
// Mar 18, 2010 Wenqin Xu
//
// V2.4 added isAncestorb to electrons Ntuples, read from event::isAncestor function
// and used them to check my result (This isAncestorb could be misused)
// added a histogram to record the original b returned by event.iTopCopyId
// replaced id-n*10000 with id-n*pow(10.,static_cast<int>(log10(id))) 
// to deal with even higher states, i.e. Jpsi(2s)
// Mar 19, 2010 Wenqin Xu
//
// V2.5 added omega (pdgcode=223) to the not wanted electron parents, because isAncestorb could have a problme
// of identifying wrong b ancestor. 
// Mar 29, 2010 Wenqin Xu
//
//==============================================================================
//  pmainHF2e.cpp
//    
//  This is an example program to study c -> e and b -> e decays
//  in 200 GeV pp collisions with Pythia8.
//
//  The decays are stored in a ROOT tree and written to file.
//
//  Once written most things can be controlled through the runcard,
//  so there's no need to recompile.
//
//  Usage: pmainHF2e  runcard  rootfile
// 
//  Author: Thomas Ullrich  
//  Last update: September 9, 2008
//==============================================================================
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "Pythia.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#define PI 3.141592654	
#define Phibin 80
#define PR(x) std::cout << #x << " = " << (x) << std::endl;




using namespace Pythia8; 
/*
struct hfedecay_t
{
    int event_no;	

    int hf_index;   
 
    int hf_id;         // mother (c/b hadron)
    int hf_status;
   
    float hf_pt;
    float hf_pz;
    float hf_phi;
    float hf_eta;
    float hf_y;
    
    int e_index;

    int e_id;            // electron
    int e_status;
   
    float e_pt;
    float e_pz;
    float e_phi;
    float e_eta;
    float e_y;
	
    int   q1_id;
    float q1_x;
    int   q2_id;
    float q2_x;
    float Q2fac;
   
    float alphas;
    float ptHat;
    int   nFinal;
    float pdf1;
    float pdf2;
    int   code;
    float sigmaGen;
    float weight;   // useful for normalization/x-section
};
*/
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

//hfedecay_t hfe_event;

particle_t asso_particle;
particle_t upsilon_daughter;

// stores data in tuple)
int Flavor(int id)
{   
	id=fabs(id);
	if(id>9900000) id-=9900000; //deal with directly gluon pair
	if(id<10000)return static_cast<int>(id/pow(10.,static_cast<int>(log10(id))));
        else	
	{// to deal with the rear excited states, e.g 10521
		int n= static_cast<int>(id/pow(10.,static_cast<int>(log10(id)))); 
		id=id-n*pow(10.,static_cast<int>(log10(id)));
		//id=id-n*10000; 
		return static_cast<int>(id/pow(10.,static_cast<int>(log10(id))));
	}
}	






int main(int argc, char* argv[]) {
    
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " runcard  rootfilejobnumber mPrint" << endl;
        return 2;
    }
    char* runcard  = argv[1];
    char rootfile[256];
    sprintf(rootfile,"%s.root",argv[2]);
    cout<<"INFO: runcard is "<<runcard<<" rootfile is "<<rootfile<<endl;
    const char* xmlDB    = "pythia8145/xmldoc";
    bool mPrint=false;
    if(strcmp(argv[3],"true")==0 )mPrint=true;
    //--------------------------------------------------------------
    //  Initialization
    //--------------------------------------------------------------
    
    //
    //  ROOT
    //
    TFile *hfile  = new TFile(rootfile,"RECREATE");
	TTree tree("tree", "c -> e decays");
	
	tree.Branch("asso_particle", &asso_particle.event_no,
				"event_no/I:"
				"id/I:pt/F:pz/F:phi/F:eta/F:y/F:"
				"mother_id/I");
    TH1F *hfid = new TH1F("hfid","id for any particles",12000,-6000,6000);
    TH1F *htally  = new TH1F("htally","tally",100,-0.5,99.5);
    TH1F *hdPhihh[6][2];
    TH1F *hdPhigh[6][2];
    TH1F *hdPhiuh[6][2];
    TH1F *hdPhidh[6][2];
    for(int trigpT=0;trigpT<6;trigpT++){			
	for(int asspT=0;asspT<2;asspT++){
    	hdPhihh[trigpT][asspT] = new TH1F(Form("hdPhihh%d_%d",trigpT,asspT),Form("hdPhihh%d_%d",trigpT,asspT),Phibin, -0.5*PI,1.5*PI);
    	hdPhigh[trigpT][asspT] = new TH1F(Form("hdPhigh%d_%d",trigpT,asspT),Form("hdPhigh%d_%d",trigpT,asspT),Phibin, -0.5*PI,1.5*PI);
    	hdPhiuh[trigpT][asspT] = new TH1F(Form("hdPhiuh%d_%d",trigpT,asspT),Form("hdPhiuh%d_%d",trigpT,asspT),Phibin, -0.5*PI,1.5*PI);
    	hdPhidh[trigpT][asspT] = new TH1F(Form("hdPhidh%d_%d",trigpT,asspT),Form("hdPhidh%d_%d",trigpT,asspT),Phibin, -0.5*PI,1.5*PI);
	}
    }
    TH1F *hPt1 = new TH1F("hPt1","Pt of trigger hadrons",110,0,11);
    TH1F *hPt = new TH1F("hPt","Pt of associated hadrons",100,0,10);

	TH1F *hHFePt = new TH1F("hHFePt", "Pt of electron from HF decay", 100, 0.0, 20.0);
	
    TH1F *hNpartons=new TH1F("hNpartons","number of partons in one event",100,-0.5,99.5);
	TH1F *hNHFpartons = new TH1F("hNHFpartons", "Number of HF partons per event", 10, -.5, 9.5);
	TH1F *hNHFelectrons = new TH1F("hNHFelectrons", "Number of electrons from HF decay per event", 10, -.5, 9.5);
    TH2F *hPartonParent=new TH2F("hPartonParent","parton index v.s. parent index;parton index;parent index",100,-0.5,99.5,100,-0.5,99.5);
	
	TH2F *hHFePtCorr = new TH2F("hHFePtCorr", "Pt correlation for hf ee events", 100, 0.0, 20.0, 100, 0.0, 20.0);
	TH2F *hHFePtCorr_acc = new TH2F("hHFePtCorr_acc", "Pt correlation for hf ee events in acceptance", 100, 0.0, 20.0, 100, 0.0, 20.0);
	TH1F *hPhiee = new TH1F("hPhiee", "Phi angle between hf ee", 100, 0.0, 2*PI);
	TH1F *hThetaee = new TH1F("hThetaee", "Theta angle between hf ee", 100, 0.0, 2*PI);
	TH1F *hHFeEta = new TH1F("hHFeEta", "Eta of HF electron", 50, -3.0, 3.0);
    //  Create instance of Pythia 
    //
    Pythia pythia(xmlDB); // the init parameters are read from xml files
	// stored in the xmldoc directory. This includes
	// particle data and decay definitions.
    
    //
    // Shorthand for (static) settings
    //
    Settings& settings = pythia.settings;
    
    //
    //  Read in runcard
    //
    pythia.readFile(runcard);  
    cout << "Runcard '" << runcard << "' loaded." << endl;
   // pythia.readString("431:mayDecay = off"); 
   //   // pythia.readString("421:oneChannel = 1 1 0 -321 211"); // to focus on one decay channel only, may affect correlations 
   //      // pythia.readString("431:33:products = 333 211"); 
   //        // pythia.particleData.readString("431:33:products = 333 211");
   //
 
    //
    //  Retrieve number of events and other parameters from the runcard.
    //  We need to deal with those settings ourself. Getting
    //  them through the runcard just avoids recompiling.
    //
    int  maxNumberOfEvents = settings.mode("Main:numberOfEvents");
    int  nList     = settings.mode("Main:numberToList");
    int  nShow     = settings.mode("Main:timesToShow");
    int  maxErrors = settings.mode("Main:timesAllowErrors");
    bool showCS    = settings.flag("Main:showChangedSettings");
    bool showAS    = settings.flag("Main:showAllSettings");
    int  pace = maxNumberOfEvents/nShow;
    	
    //
    //  Remark: in this example we do NOT alter the
    //  BRs since they are different for the various charm
    //  hadrons making the normalization at the end rather
    //  cumbersome. In a production version this is what
    //  one probably would implement to save processing time.
    //
    
    //
    //  Initialize Pythia, ready to go
    //
    pythia.init();
    
    //
    // List changed or all data
    //
    if (showCS) settings.listChanged();
    if (showAS) settings.listAll();
    
    //--------------------------------------------------------------
    //  Event loop
    //--------------------------------------------------------------
    int ievent = 0;
    int numberOfElectrons = 0;
    int iErrors = 0;
    const float pt_trig_up[6] = {6,7,8,9,10,11};           
    const float pt_trig_lo[6] = {5,6,7,8,9,10};           
    const float pt_asso_up[2] = {0.5,1};
    const float pt_asso_lo[2] = {0.15,0.5};
    while (ievent < maxNumberOfEvents) 
    {
        
        if (!pythia.next()) {
            if (++iErrors < maxErrors) continue;
            cout << "Error: too many errors in event generation - check your settings & code" << endl;
            break;
        }

	// ===============
	//  Event analysis, define eventwise parameters firstly.
	// ===============

	Event &event = pythia.event;
	int flavor =0;
	bool b_hfevent=false;
	bool is_hfe_event=false;
	int nhfe = 0;
	int npartons=0;
	int nhfpartons=0;
	htally->Fill(0);
	if(mPrint)cout<<"to start the 1st loop"<<endl;
	for (int i = 0; i < event.size(); i++) 
	{
		htally->Fill(1);
		hfid->Fill(event[i].id());	
		if(event[i].isFinal()==true) continue;//not trying to find final particles now
		if(event[event[i].mother1()].id()==event[i].id())continue; //recoil
		if(Flavor(event[i].id()) == 4 || Flavor(event[i].id())==5 )	//
		{
			htally->Fill(2);
			if(fabs(event[i].id())==5 )htally->Fill(3);
			b_hfevent=true;
			nhfpartons++;
			for(int idaughter = event[i].daughter1(); idaughter <= event[i].daughter2(); idaughter++) {
				if(fabs(event[idaughter].id()) == 11) {
					if(mPrint) cout << "HF e found" << endl;
					hHFePt->Fill(event[idaughter].pT());
					nhfe++;
					is_hfe_event = true;
				}
			}
			hPartonParent->Fill(i,event[i].mother1());
		
		}
		if(fabs(event[i].id()) == 1 || fabs(event[i].id())==2 || fabs(event[i].id())==3 ||fabs(event[i].id())==21)
		{
			npartons++;
			hPartonParent->Fill(i,event[i].mother1());	
		}
	}//the main purpose of first loop is to decide whether this is a heavy flavor event or not, 

	hNHFpartons->Fill(nhfpartons);	
	hNHFelectrons->Fill(nhfe);	
		
	//loop over final state particles and record
	if(is_hfe_event) {
		for (int i=0; i<event.size(); i++) {
			if (event[i].isFinal() && event[i].isCharged()) {
				asso_particle.event_no = ievent;
				
				asso_particle.id = event[i].id();
				asso_particle.pt = event[i].pT();
				asso_particle.pz = event[i].pz();
				asso_particle.phi = event[i].phi();
				asso_particle.eta = event[i].eta();
				asso_particle.y = event[i].y();
				
				asso_particle.mother_id = event[event[i].mother1()].id();
				
				tree.Fill();
			}
		}
	}
		
	//loop over final state electrons
	
	//hfedecay_t hfedecay[2];	
	int index_e = 0;	
	/*	
	for(int i=0; i<event.size(); i++) {
		
		if (!(event[i].isFinal()&&fabs(event[i].id()) == 11)) continue;
		
		int imother = event[i].mother1();
		
		if(Flavor(event[imother].id())==4||Flavor(event[imother].id())==5) { //is HF
		
			//record electron
			hfe_event.event_no = ievent;
		
			hfe_event.hf_index  = imother;
	
			hfe_event.hf_id     = event[imother].id();  
			hfe_event.hf_status = event[imother].status();
			hfe_event.hf_pt     = event[imother].pT();
			hfe_event.hf_pz     = event[imother].pz();
			hfe_event.hf_phi    = event[imother].phi();
			hfe_event.hf_eta    = event[imother].eta();     
			hfe_event.hf_y      = event[imother].y();

			hfe_event.e_index   = i;			

			hfe_event.e_id       = event[i].id();  	   
			hfe_event.e_status   = event[i].status();	
			hfe_event.e_pt 	   = event[i].pT();	
			hfe_event.e_pz 	   = event[i].pz();	
			hfe_event.e_phi 	   = event[i].phi();	
			hfe_event.e_eta      = event[i].eta();     
			hfe_event.e_y 	   = event[i].y();
			
			hfe_event.q1_id      = pythia.info.id1();
			hfe_event.q1_x       = pythia.info.x1();
			hfe_event.q2_id      = pythia.info.id2();
			hfe_event.q2_x       = pythia.info.x2();
			hfe_event.Q2fac      = pythia.info.Q2Fac();
			hfe_event.alphas     = pythia.info.alphaS();
			hfe_event.ptHat      = pythia.info.pTHat();
			hfe_event.nFinal     = pythia.info.nFinal();
			hfe_event.pdf1       = pythia.info.pdf1();
			hfe_event.pdf2       = pythia.info.pdf2();
			hfe_event.code       = pythia.info.code();
			hfe_event.sigmaGen   = pythia.info.sigmaGen();
			hfe_event.weight     = 0.0;	
			
			tree.Fill();
			
		if(index_e >= 2) {
			cout << "Too many electrons to handle" << endl;
			break;
		}

		//record electron information
		hfedecay[index_e].event_no = ievent;
	
	        hfedecay[index_e].hf_index = imother;
	
		hfedecay[index_e].hf_id     = event[imother].id();  
		hfedecay[index_e].hf_status = event[imother].status();
		hfedecay[index_e].hf_pt     = event[imother].pT();
		hfedecay[index_e].hf_pz     = event[imother].pz();
		hfedecay[index_e].hf_phi    = event[imother].phi();
		hfedecay[index_e].hf_eta    = event[imother].eta();     
		hfedecay[index_e].hf_y      = event[imother].y();
		
                hfedecay[index_e].e_index   = i;

		hfedecay[index_e].e_id       = event[i].id();  	   
		hfedecay[index_e].e_status   = event[i].status();	
		hfedecay[index_e].e_pt 	   = event[i].pT();	
		hfedecay[index_e].e_pz 	   = event[i].pz();	
		hfedecay[index_e].e_phi 	   = event[i].phi();	
		hfedecay[index_e].e_eta      = event[i].eta();     
		hfedecay[index_e].e_y 	   = event[i].y();
		
		hfedecay[index_e].q1_id      = pythia.info.id1();
		hfedecay[index_e].q1_x       = pythia.info.x1();
		hfedecay[index_e].q2_id      = pythia.info.id2();
		hfedecay[index_e].q2_x       = pythia.info.x2();
		hfedecay[index_e].Q2fac      = pythia.info.Q2Fac();
		hfedecay[index_e].alphas     = pythia.info.alphaS();
		hfedecay[index_e].ptHat      = pythia.info.pTHat();
		hfedecay[index_e].nFinal     = pythia.info.nFinal();
		hfedecay[index_e].pdf1       = pythia.info.pdf1();
		hfedecay[index_e].pdf2       = pythia.info.pdf2();
		hfedecay[index_e].code       = pythia.info.code();
		hfedecay[index_e].sigmaGen   = pythia.info.sigmaGen();
		hfedecay[index_e].weight     = 0.0;//pythia.info.sigmaGen()/nMaxEvent; //not using zeroed out to avoid errors
		
		index_e++;
		
		if(index_e==2) {
			cout << "2 hfe event found" << endl;
			hHFePtCorr->Fill(hfedecay[0].e_pt, hfedecay[1].e_pt);
			float p1mag, p2mag, p1dotp2;
			p1mag = sqrt(event[hfedecay[0].e_index].px()*event[hfedecay[0].e_index].px() +
				     event[hfedecay[0].e_index].py()*event[hfedecay[0].e_index].py() +
				     event[hfedecay[0].e_index].pz()*event[hfedecay[0].e_index].pz());
			p2mag = sqrt(event[hfedecay[1].e_index].px()*event[hfedecay[1].e_index].px() +
                                     event[hfedecay[1].e_index].py()*event[hfedecay[1].e_index].py() + 
				     event[hfedecay[1].e_index].pz()*event[hfedecay[1].e_index].pz());
			p1dotp2 = event[hfedecay[0].e_index].px()*event[hfedecay[1].e_index].px() +
                                  event[hfedecay[0].e_index].py()*event[hfedecay[1].e_index].py() + 
				  event[hfedecay[0].e_index].pz()*event[hfedecay[1].e_index].pz();
			float cosTheta = p1dotp2/(p1mag*p2mag);
			float Theta = acos(cosTheta);
			hThetaee->Fill(Theta);
			float deltaPhi = hfedecay[0].e_phi - hfedecay[1].e_phi;
			if(deltaPhi>=0.0&&deltaPhi<PI) {
				hPhiee->Fill(deltaPhi);
			}
			else if(deltaPhi<0.0&&fabs(deltaPhi)<PI) {
				deltaPhi *= -1.0;
				hPhiee->Fill(deltaPhi);
			}
			else if(deltaPhi>0.0&&deltaPhi>PI) {
				deltaPhi = 2*PI - deltaPhi;
				hPhiee->Fill(deltaPhi);
			}
			else {
				deltaPhi += 2*PI;
				hPhiee->Fill(deltaPhi);
			}
			if(fabs(hfedecay[0].e_eta)<=1.0&&fabs(hfedecay[1].e_eta)<=1.0) {
				hHFePtCorr_acc->Fill(hfedecay[0].e_pt, hfedecay[1].e_pt);
			}
		}
		if(index_e>2) {
			cout << "More than 2 HFe ignoring for now" << endl;	
		}
			
		}
		
	}
	*/
		
        ievent++;
		/*
        if (ievent%pace == 0) {
            cout << "# of events generated = " << ievent 
            << ", # of electrons from c or b hadron decays generated so far = " << numberOfElectrons << endl;
        }
        */
        // List first few events.
        if (ievent < nList) {
            pythia.info.list(); 
            pythia.process.list();
            pythia.event.list();
        }
	if(mPrint)cout<<"event "<<ievent-1<<" is done"<<endl;
    }//end of all events, i.e. end of the program
    
    //--------------------------------------------------------------
    //  Finish up
    //--------------------------------------------------------------
    pythia.statistics();
    cout << "End of the program. Writing File" << endl;
    hfile->Write();
    
    return 0;
}


