  /* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
*                   --<--<--  A fast simulator --<--<--     *
*                 / --<--<--     of particle   --<--<--     *
*  ----HECTOR----<                                          *
*                 \ -->-->-- transport through -->-->--     *
*                   -->-->-- generic beamlines -->-->--     *
*                                                           *
* JINST 2:P09005 (2007)                                     *
*      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *
*       http://www.fynu.ucl.ac.be/hector.html               *
*                                                           *
* Center for Cosmology, Particle Physics and Phenomenology  *
*              Universite catholique de Louvain             *
*                 Louvain-la-Neuve, Belgium                 *
 *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/// \file H_Acceptance.cpp
/// \brief Collection of routines about physics output from Hector, excepting explicit variable reconstruction (see also H_Reconstruction.cpp)

// C++ includes
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>

// local includes
#include "H_AbstractBeamLine.h"
#include "H_Beam.h"
#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_OpticalElement.h"
#include "H_RomanPot.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraph2D.h"
using namespace std;

static int hN =0; // histogram counter -- for a unique TH1 histo name

void acceptance_basiccheck(double emin, double emax, int steps) {
	extern bool relative_energy;
	relative_energy = true;
	extern int kickers_on;
	kickers_on = 1;
	H_BeamLine* beam = new H_BeamLine(1,500.);
	beam->fill("data/LHCB1IR5_v6.500.tfs");
	for(int i = 0; i < steps+1; i++) {
		cout<<"------------------------------------------"<<endl;
		cout<<"energy : "<<emin + i*(emax-emin)/((double)steps)<<endl;
		H_BeamParticle p;
		//p.setPosition(0,-500,0,-142.5,0); / atlas conditions
                p.setPosition(-500,0,-142.5,0,0);	
		p.emitGamma(emin + i*(emax-emin)/((double)steps),0);
		p.computePath(beam);
		if(p.stopped(beam)) {
			p.getStoppingElement()->printProperties();
		}
	}
	delete beam;
}

/// Draws the 1D-acceptance for any optical element, for a given energy range and a given virtuality
TH1F* acceptance_element(H_AbstractBeamLine * beamline, const string element_name, const unsigned int number_of_particles=10, const float E_min=0., const float E_max=150., const unsigned int Ne = 10, const float Q2 =-0.1, const bool NonLinear=true, const int crang_sign  = -1, bool drawlines = false) {
/// @param beamline : beamline to test
/// @param number_of_particles : in the beam
/// @param E_min : min value for the beam particle energy loss
/// @param E_max : max value for  the beam particle energy loss
/// @param Ne : number of steps in energy
/// @param Q2 : virtuality value (t) of the emitted gamma
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
	const bool show_path = false;

	vector<TVectorD> stop_pos;
	vector<TVectorD>::iterator spi_i;

	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {	
		cout << "ERROR : You should be in relative energy" << endl;
		return new TH1F();
	}
	if(!beamline) {
		cout << "ERROR : No beamline ! Empty pointer.\n";
		return new TH1F();
	}
	if(!element_name.size()) {
		cout << "ERROR : Empty element name !\n";
		return new TH1F();
	}
	if(!NonLinear) cout << "WARNING : Linear mode (NonLinear off).\n";

        float energies[Ne]; //, * log10x = new float[Ne];
	for(unsigned int e_i=0; e_i<Ne; e_i++) {energies[e_i] = E_min + (E_max - E_min)/(float)Ne * (e_i+0.5);}

	char can_name[100], title[20];
	sprintf(can_name,"Acceptance of %s",element_name.c_str());
	sprintf(title,"h1D_%d",hN++);
	TH1F * h1D = new TH1F(title,can_name,Ne,E_min,E_max); // 1D-histo at fixed Q^2 = 0
	TMultiGraph* hitx_m = new TMultiGraph("hitx_m","hitx_m");
	TMultiGraph* hity_m = new TMultiGraph("hity_m","hity_m");

	cout << "Generating particles " << endl; 
	{
		vector<H_OpticalElement> stopping_elements;
		vector<H_OpticalElement>::iterator el_i;
		vector<int> number_of_stopped_particles;

		extern int kickers_on;
		kickers_on =1;

		for (unsigned int e_i=0; e_i<Ne; e_i++) {
			H_Beam mybeam;
			mybeam.setPosition(PX,PY,crang_sign*CRANG,0.,0.);
			mybeam.setE(BE_DEF);
			if(drawlines) {
				mybeam.setDE(0);
				mybeam.setDX(0);
				mybeam.setDY(0);
				mybeam.setDTX(0);
				mybeam.setDTY(0);
			}
			mybeam.createBeamParticles(number_of_particles); // smearings

			for(unsigned int particle_i=0; particle_i<number_of_particles; particle_i++) {
				mybeam.getBeamParticle(particle_i)->emitGamma(energies[e_i],Q2); 
			}
			mybeam.computePath(beamline,NonLinear);
			if(drawlines) {
//				hitx_m->Add(mybeam.drawBeamX(1),"l");
//				hity_m->Add(mybeam.drawBeamY(1),"l");
			}

			vector<TVectorD> stop_pos_i = mybeam.getStoppingElements(beamline,stopping_elements,number_of_stopped_particles);
			for(spi_i = stop_pos_i.begin(); spi_i<stop_pos_i.end(); spi_i++) {
	//			if((*spi_i)[4]>=beamline->getElement(element_name)->getS() && (*spi_i)[4]<=beamline->getElement(element_name)->getS()+beamline->getElement(element_name)->getLength()) {
					stop_pos.push_back(*spi_i);
	//			}
			}

			for(el_i = stopping_elements.begin(); el_i < stopping_elements.end(); el_i++) {
				if((el_i->getName().c_str() == element_name) &&  (mybeam.getBeamParticle(1)->isPhysical())) 
					h1D->Fill(energies[e_i],number_of_stopped_particles[el_i-stopping_elements.begin()]/(float)(number_of_particles));
			}
		} // for e_i
	
        }// kills stopping_elements and other vector<...>'s
	cout << "\t done"<<endl;

	TGraph* svse = new TGraph();
	TGraph2D* hits = new TGraph2D();
	TGraph* hitx = new TGraph();
	if(drawlines) hitx->SetMarkerColor(2);
	TGraph* hity = new TGraph();
	if(drawlines) hity->SetMarkerColor(2);
	TGraph* hitp = new TGraph();
	for(spi_i = stop_pos.begin(); spi_i<stop_pos.end(); spi_i++) {
		svse->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[4],energies[(spi_i - stop_pos.begin())/(int)number_of_particles]);
		hits->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[0],(*spi_i)[4],(*spi_i)[2]);
		hitx->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[4],(*spi_i)[0]);
		hity->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[4],(*spi_i)[2]);
		hitp->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[0],(*spi_i)[2]);
	}
	if(drawlines) {
		hitx->SetMarkerStyle(7);
		hitx->SetMarkerColor(2);
		hity->SetMarkerStyle(7);
		hity->SetMarkerColor(2);
	}
	hitx_m->Add(hitx,"p");
	hity_m->Add(hity,"p");
	TCanvas* can = new TCanvas("can","can",1);
	can->Divide(2,2);
	can->cd(1);
	hits->Draw("P");
	can->cd(2);
	hitx_m->Draw("ap");
	can->cd(3);
	hity_m->Draw("ap");
	can->cd(4);
	hitp->Draw("ap");
	TCanvas* ca2 = new TCanvas("ca2","ca2",1);
	ca2->SetGrid();
	svse->Draw("ap+");
	double height = (E_max - E_min)/10.;
	relative_energy = false;
	beamline->offsetElements(0,5*height/URAD);
	beamline->drawX(1*height,1*height,1);
	beamline->draw();
	svse->SetMarkerStyle(6);
	svse->GetXaxis()->SetTitle("s (m)");
	svse->GetYaxis()->SetTitle("Eloss (GeV)");
	svse->Draw("p+");

	// shows the beam path, just in case
	if (show_path) {
		H_BeamParticle p1;
		p1.setE(BE_DEF);
		p1.setPosition(p1.getX()+PX,p1.getY()+PY,p1.getTX()+crang_sign*CRANG,p1.getTY(),0.);
		p1.computePath(beamline,NonLinear);
		TGraph * ppath = p1.getPath(0,kGreen);
		TCanvas * path = new TCanvas("path","path");
		path->cd();
		ppath->Draw("ALP");
	}
	return h1D;
} // acceptance_element

/// Draws the lego plot of irradiation levels, from \f$ pp \rightarrow pX \f$
void acceptance_fluence(float rpos = 220., float deltax = 2000., string file = "data/LHCB1IR5_v6.500.tfs", int side = 1, const int N = 1000, const bool save=false) {
	/// @param rpos is the s coordinate where to look at the irradiation, in [m] from IP
	/// @param deltax is the x coordinate where a fake roman pot is set, in [microm] from the beam center
	/// @param file is the optics table
	/// @param side is the direction of propagation (1 is downstream, -1 is upstream)
	/// @param N is the number of particles to propagate (ideal one : 10^6)
	/// @param save : boolean

#if !defined(_include_pythia_) && !defined(_include_pythia8_)
	cout << "You are not running Pythia in Hector (_include_pythia(8)_ is not defined)!\n";
	return;
#endif

	gStyle->SetPalette(8);
	extern bool relative_energy;
	relative_energy = true;
	extern int kickers_on;
        kickers_on = 1;
	if(!relative_energy) {
		cout << "You should be in relative energy" << endl;
	}
	
        gROOT->SetStyle("Plain");
        // cross-section of the process [mb]
        const float sigma = 7.15;
        // LHC lumi in mbarn-1 s-1
        const float lumi = 2E6; // = 2x10^33 cm^-2 s^-1
	// Bunch crossing rate
	//const float bxing = 31536000.; // number of seconds in 1 year
	const float bxing = 1E7; // to achieve 20 fb-1
        // weight of the events : sigma * lumi / BC_rate / N
        float weight = sigma*lumi*bxing/(float)N; // per year
        weight = weight*100*2; // converting bins into cm^2

	char title[100];
	sprintf(title,"Hits in VFD at %dm (L=%.0f fb^{-1}",(int)rpos,lumi*bxing*1E-12);
	int binx = (rpos==220.) ? 20 : 20;
	int biny = (rpos==220.) ? 20 : 20;
	float ymin = (rpos==220.) ? -3.5 : -1.50; 
	float ymax = (rpos==220.) ? 3.5 : 1.5;
	float xmin = (rpos==220.) ? 0 : -30;
	float xmax = (rpos==220.) ? 20 : 0;
        TH2F* rp_hits = new TH2F("rp_hits",title,binx,xmin,xmax,biny,ymin,ymax);
        int stopped_number = 0;

        gROOT->SetStyle("Plain");
        gRandom->SetSeed(0);
        H_BeamLine* beam = new H_BeamLine(side,rpos+5);

        beam->fill(file);
        H_RomanPot * rp = new H_RomanPot("RP",rpos,deltax,BE_DEF);
        beam->add(rp);
        beam->offsetElements(140,-0.08);
        beam->offsetElements(146,-0.097);

	for(int i=0;i<N;i++) {
		H_BeamParticle p1;
		p1.smearAng();
		p1.smearPos();
		p1.setPosition(p1.getX()-500.,p1.getY(),p1.getTX()-1*kickers_on*CRANG,p1.getTY(),0);
		p1.doInelastic();
		p1.computePath(beam);

		if(p1.stopped(beam)) 
			if(p1.getStoppingElement()->getName()>="RP") {
                                                stopped_number++;
                                                p1.propagate(rpos);
                                                rp_hits->Fill(p1.getX()/1000.,p1.getY()/1000.,weight);
                        }
        }

	TCanvas * mycan = new TCanvas();
        rp_hits->Draw("SURF5");
	rp_hits->GetXaxis()->SetTitle("x (mm)");
	rp_hits->GetYaxis()->SetTitle("y (mm)");
	rp_hits->GetZaxis()->SetTitle("Hits per year per cm^{2}");
	char name[100];
	sprintf(name,"hits_%d_%d_%d.root",(int)rpos,(int)deltax,N);
	if(save) {
		mycan->Print(name);
		delete rp_hits;
		delete mycan;
	}
	cout << "Assumptions : \n\t* luminosity = " << lumi << " mb-1 s-1\n";
	cout << "\t* integrated luminosity = " << lumi*bxing*1E-12 << " mb-1 s-1\n";
	cout << "\t* Cross section : " << sigma << " mb\n";
	cout << "\t* Bunch crossing rate : " << bxing << endl;
	
}

/// Returns the usual beamline with a RP
H_AbstractBeamLine * createBeamline(float rp_pos, float rp_x, string filename, int side, string& rp_name) {
	/// @param rp_pos is the RP position in [m], from the IP
	/// @param rp_x is the RP position in [\f$ \mu m \f$], when the nominal beam is assumed to be in 0.
	/// @param filename : optics source file
	/// @param side is the direction of propagation (1 : forward, -1 : backward)
	/// @param rp_name : returns the RP name, for subsequent use of this specific optical element
        int direction =  (side<0)?-1:1;
        H_BeamLine* beamline = new H_BeamLine(side,rp_pos+0.1);
        beamline->fill(filename);
	char title[30];
	sprintf(title,"RP %dm %dum",(int)rp_pos,(int)rp_x);
        H_RomanPot * rp = new H_RomanPot(title,rp_pos,rp_x,BE_DEF);
        rp_name = title;
        beamline->add(rp);
        beamline->offsetElements(140,-1*(direction)*0.080);
        beamline->offsetElements(146,-1*(direction)*0.097);
        beamline->offsetElements(120,-1*(direction)*0.097);
        return beamline;
}

/// Draws the acceptance for each element of the beamline hit by some particles, from the IP until the RP.
void acceptance_beamline(float rp_pos, float rp_x, string filename="data/LHCB1IR5_v6.500.tfs", int side=1, int crang_sign=-1, float phimin = 0, float phimax = 2*3.141592653589, const float log10x_min=-3.5, const float log10x_max=-0.3, const unsigned int Ne = 20, const float log10t_min=-3, const float log10t_max=0, const unsigned int Nq = 33, unsigned int number_of_particles=30, bool NonLinear = true, bool draw_RP_only=false) {
/// @param rp_pos : distance of the roman pots, in meters from the IP
/// @param rp_x   : roman pot transverse position from the center of the beam, in micrometers [\f$ \mu m \f$]
/// @param filename : name of the optics source file
/// @param side : direction of propagation (+1 : forward, -1: backward)
/// @param crang_sign : sign of the (half) crossign angle (should be -1 when forward)
/// @param phimin : lower bound for the random phi angle, in H_Particle::emitGamma()
/// @param phimax : higher bound for the random phi angle, in H_Particle::emitGamma()
/// @param Ne : number of steps in energy
/// @param log10x_min : min value for x (\f$ x = \xi \f$)
/// @param log10x_max : max value for x (\f$ x = \xi \f$)
/// @param Nq : number of steps in t
/// @param log10t_min : min of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param log10t_max : max of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param number_of_particles : in the beam
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
/// @param draw_RP_only : draw option

	extern bool relative_energy;
	relative_energy=true;
	if(!relative_energy) {
		cout << "You should be in relative energy" << endl;
		return;
	}

        gROOT->SetStyle("Plain");
        H_AbstractBeamLine * beamline;
        string rp_name;
        beamline = createBeamline(rp_pos+0.1,rp_x,filename,side,rp_name);

        vector<string> full_el_list;
        vector<float> full_n_list;
        vector<float> full_en_list;
        vector<float> full_q2_list;

        float * energies = new float[Ne], * log10x = new float[Ne];
        float * q2s = new float[Nq], *log10t = new float[Nq];
        float * inefficiency = new float[Ne*Nq];

	for(unsigned int log10x_i=0; log10x_i<Ne; log10x_i++) log10x[log10x_i]=log10x_min + (log10x_max - log10x_min)/(float)Ne*(log10x_i+0.5);
	for(unsigned int log10t_i=0; log10t_i<Nq; log10t_i++) log10t[log10t_i]=log10t_min + (log10t_max - log10t_min)/(float)Nq*(log10t_i+0.5);
        for(unsigned int e_i=0; e_i<Ne; e_i++) { energies[e_i]= BE_DEF * pow(10,log10x[e_i]); }
        for(unsigned int q2_i=0; q2_i<Nq; q2_i++) {q2s[q2_i]=  -1 * pow(10,log10t[q2_i]); }
	TMultiGraph * tm =0;

	cout << "Generating particles : " << endl; 
        {
                vector<H_OpticalElement> stopping_elements;
                vector<H_OpticalElement>::iterator el_i = stopping_elements.begin();
                vector<int> number_of_stopped_particles;

		extern int kickers_on;
		kickers_on =1;

                for (unsigned int q2_i=0; q2_i<Nq; q2_i++) {
                        for (unsigned int e_i=0; e_i<Ne; e_i++) {
                                H_Beam mybeam;
				mybeam.setPosition(-500.,0.,crang_sign*CRANG,0.,0.);
                                mybeam.createBeamParticles(number_of_particles);
				mybeam.emitGamma(energies[e_i],q2s[q2_i],phimin,phimax);

				mybeam.computePath(beamline,NonLinear);
				mybeam.getStoppingElements(beamline,stopping_elements,number_of_stopped_particles);
                                inefficiency[q2_i*Ne+e_i]=1-(mybeam.getStoppedNumber(beamline)/(float)number_of_particles);

                                for(el_i = stopping_elements.begin(); el_i < stopping_elements.end(); el_i++) {
                                        full_el_list.push_back(el_i->getName());
					if(mybeam.getBeamParticle(1)->isPhysical()) {
  	                                    	full_n_list.push_back(number_of_stopped_particles[el_i-stopping_elements.begin()]/(float)(number_of_particles));
					} else {
						full_n_list.push_back(0);
					}
                                        full_en_list.push_back(log10x[e_i]);
                                        full_q2_list.push_back(log10t[q2_i]);
                                        // creates 4 vectors: list of stopping elements and the corresponding number of stopped particles, for a given E and Q2
                                        // eg : full_el_list  = RP, TAN.4R5, RP
                                        //      full_n_list   = 10,  2     , 8
                                        //      full_en_list  = 500, 600   , 600
                                        //      full_q2_list  = -1, -1     , -1
                                } // for el_i
                        } // for e_i
                }//for q2_i
        }// kills stopping_elements and other vector<...>'s
	cout << "\t done"<<endl;

        vector<string> list_of_names;
        vector<string>::iterator li = list_of_names.begin();
        vector<string>::iterator el_l = list_of_names.begin();

        {
                vector<string> temp = full_el_list;
                sort(temp.begin(),temp.end());
                list_of_names.push_back(*(temp.begin()));
                for (el_l=temp.begin(); el_l<temp.end(); el_l++)
                        if(*(list_of_names.end()-1) != *el_l) list_of_names.push_back(*el_l);
        } // kills the temp vector

        cout << endl << "The elements that are hit by some particles are : " <<endl;
        for(el_l=list_of_names.begin(); el_l<list_of_names.end(); el_l++) cout << "\t" << *el_l << endl;
        cout << endl;

	TCanvas * path = new TCanvas("path","Beam without energy loss, for monitoring purpose");
	path->cd();
		H_Beam mybeam;
		mybeam.setPosition(-500.,0.,crang_sign*CRANG,0.,0.);
		mybeam.createBeamParticles(number_of_particles);
		mybeam.computePath(beamline,NonLinear);
		tm = mybeam.drawBeamX(kRed);
	tm->Draw("ALP");
	gPad->SetGrid();

        TCanvas * main_can = new TCanvas("main_can","");
        int cur_pad=1;
		if(draw_RP_only) main_can->Divide(2);
        else main_can->Divide(list_of_names.size()+1);

        TH2F * h2D;
        TCanvas * eff_can;
        char can_name[20];
        for(li = list_of_names.begin(); li<list_of_names.end(); li++) {
                strcpy(can_name,(*li).c_str());
				h2D = new TH2F(can_name,can_name,Nq,log10t_min,log10t_max,Ne,log10x_min,log10x_max);

                for(el_l=full_el_list.begin(); el_l<full_el_list.end(); el_l++) {
                        int i = el_l-full_el_list.begin();
                        if(*el_l == *li) {
				// fill (log Q, log x, weight)
				h2D->Fill(full_q2_list[i],full_en_list[i],full_n_list[i]);
                }
        }

		if(draw_RP_only && ((*li).c_str()!=rp_name)) continue; // draws only RP graph by default
	        eff_can = new TCanvas(can_name,can_name);
                h2D->SetTitle(can_name);
               	h2D->GetYaxis()->SetTitle("log_{10}(x)");
	        h2D->GetXaxis()->SetTitle("log_{10}(-Q^{2})");
        	h2D->GetZaxis()->SetTitle("Acceptance");

	        eff_can->cd();
        	h2D->Draw("surf1");
                main_can->cd(cur_pad);
	        cur_pad++;
		h2D->SetContour(10);
                h2D->Draw("Cont4");
	
		// save RP data
		if(li->c_str() == rp_name) {
			char outfilename[40];
			sprintf(outfilename,"RP_acceptance_%d_%d_%dx%d_%d_phi_%.2f_%.2f.txt",(int)rp_pos,(int)rp_x,Ne,Nq,number_of_particles,phimin,phimax);
			ofstream rp_acceptancefile(outfilename);
			for(unsigned int log10x_i=Ne; log10x_i>0; log10x_i--) {
				for(unsigned int log10t_i=0; log10t_i<Nq; log10t_i++) 
					rp_acceptancefile << fixed << setprecision(4) << h2D->GetBinContent(log10t_i+1,log10x_i) << " ";
				rp_acceptancefile << endl;
			}
			rp_acceptancefile.close(); 
		}

        } // for li

        strcpy(can_name,"particle not stopped");
		h2D = new TH2F(can_name,can_name,Nq,log10t_min,log10t_max,Ne,log10x_min,log10x_max);
	
        for (unsigned int q2_i=0; q2_i<Nq; q2_i++) {
                for (unsigned int e_i=0; e_i<Ne; e_i++) {
			h2D->Fill(log10t[q2_i],log10x[e_i],inefficiency[q2_i*Ne+e_i]);
                }
        }
        h2D->SetTitle(can_name);
        h2D->GetYaxis()->SetTitle("log_{10}(x)");
        h2D->GetXaxis()->SetTitle("log_{10}(t)");
        h2D->GetZaxis()->SetTitle("Acceptance");
	if(!draw_RP_only) {
		eff_can = new TCanvas(can_name,can_name);
		eff_can->cd();
        	h2D->Draw("surf1");
	}
	main_can->cd(cur_pad);
        cur_pad++;
        h2D->Draw("Cont4");

} // acceptance 

/// Same routine as acceptance_beamline, but dedicated to RP only. Should be slightly faster.
void acceptance_rp(float rp_pos, float rp_x, string filename="data/LHCB1IR5_v6.500.tfs", const int side=1, const int crang_sign = -1, unsigned int number_of_particles=10, const float E_min=0, const float E_max=3000, const unsigned int Ne = 10, const float log10t_min=-3, const float log10t_max=0, const unsigned int Nq = 10, bool NonLinear=true, const char * add2title="", string outfilename="") {
/// @param rp_pos : distance of the roman pots, in meters from the IP
/// @param rp_x   : roman pot transverse position from the center of the beam, in micrometers [\f$\mu m\f$]
/// @param filename : name of the optics source file
/// @param side : direction of propagation (+1 : forward, -1: backward)
/// @param crang_sign : sign of the (half) crossign angle (should be -1 when forward)
/// @param number_of_particles : in the beam
/// @param E_min : min value for the beam particle energy loss
/// @param E_max : max value for  the beam particle energy loss
/// @param Ne : number of steps in energy
/// @param log10t_min : min of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param log10t_max : max of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param Nq : number of steps in t
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
/// @param add2title : extra string added to the histogram title
/// @param outfilename : name of file to save to

	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {	
		cout << "You should be in relative energy" << endl;
		return;
	}

        gROOT->SetStyle("Plain");
	gRandom->SetSeed(0);

        H_AbstractBeamLine * beamline;
        string rp_name;
        beamline = createBeamline(rp_pos+0.1,rp_x,filename,side,rp_name);

        float * energies = new float[Ne]; //, * log10x = new float[Ne];
        float * q2s = new float[Nq], *log10t = new float[Nq];

	for(unsigned int log10t_i=0; log10t_i<Nq; log10t_i++) { log10t[log10t_i]=log10t_min + (log10t_max - log10t_min)/(float)Nq*(log10t_i+0.5);}
	for(unsigned int e_i=0; e_i<Ne; e_i++) {energies[e_i] = E_min + (E_max - E_min)/(float)Ne * (e_i+0.5); cout << energies[e_i] << " \t"; } cout << endl;

        for(unsigned int q2_i=0; q2_i<Nq; q2_i++) {q2s[q2_i]=  -1 * pow(10,log10t[q2_i]); }


	char * can_name = new char[100];
	sprintf(can_name,"Acceptance at %dm (%d #mum) for %s",(int)rp_pos,(int)rp_x,add2title);
	TH2F * myRP= new TH2F("myRP",can_name,Nq,log10t_min,log10t_max,Ne,E_min,E_max);
//	TMultiGraph * tm = 0;
//	TGraph * ppath = 0;

	cout << "Generating particles : " << endl; 
	{
		vector<H_OpticalElement> stopping_elements;
		vector<H_OpticalElement>::iterator el_i;
		vector<int> number_of_stopped_particles;

		extern int kickers_on;
		kickers_on =1;

		for (unsigned int q2_i=0; q2_i<Nq; q2_i++) {
			if (q2_i ==0) cout << "E_beam = " << endl;
			for (unsigned int e_i=0; e_i<Ne; e_i++) {

				H_Beam mybeam;
				mybeam.setPosition(-500.,0.,crang_sign*CRANG,0.,0.);
				mybeam.setE(7000.);
				mybeam.createBeamParticles(number_of_particles); // all smearings
				mybeam.emitGamma(energies[e_i],q2s[q2_i]);
				mybeam.computePath(beamline,NonLinear);

				mybeam.getStoppingElements(beamline,stopping_elements,number_of_stopped_particles);
//					if(q2_i==1 && e_i==1) tm = mybeam.drawBeamX(1);

				for(el_i = stopping_elements.begin(); el_i < stopping_elements.end(); el_i++) {
					if(el_i->getName().c_str() == rp_name) {
						if(mybeam.getBeamParticle(1)->isPhysical()) myRP->Fill(1,1);
						if(mybeam.getBeamParticle(1)->isPhysical())
							 myRP->Fill(log10t[q2_i],energies[e_i],number_of_stopped_particles[el_i-stopping_elements.begin()]/(float)(number_of_particles));
					} 
				} 
			} // for e_i

			if (q2_i == 0) cout << endl;
			cout << 100*(q2_i+1)/(float)Nq << "% completed" << endl;
	
		}//for q2_i
		cout<<"deleting pointers"<<endl;
        }// kills stopping_elements and other vector<...>'s
	cout << "\t done"<<endl;

	// shows the beam path, just in case
//	TCanvas * path = new TCanvas("path","path");
//		path->cd();
//	tm->Add(ppath);
//	tm->Draw("ALP");

	// draws the contour graph
	TCanvas * eff_can = new TCanvas("cont2",can_name);
	myRP->Draw();
	myRP->GetYaxis()->SetTitle("E_{loss} (GeV)");
	myRP->GetXaxis()->SetTitle("log_{10}(t)");
	myRP->GetXaxis()->SetTitleOffset(1.18);
	myRP->GetYaxis()->SetTitleOffset(1.2);
	myRP->GetZaxis()->SetTitleOffset(1.2);
	myRP->GetZaxis()->SetTitle("Acceptance");
	double * levels = new double[5];
	for (int i=0; i<5; i++) levels[i] = i*0.25;
	myRP->SetContour(5,levels);
	myRP->SetStats(0); // no statistics drawn
	myRP->Draw("Cont2");

	if(outfilename.size()) {  // i.e. : if(save_into_a_file)
		eff_can->Print(outfilename.c_str());
		delete eff_can;
	}

	// draws the lego plot
	eff_can = new TCanvas("lego",can_name);
	myRP->Draw("Lego");

	if (outfilename.size()) { // i.e. : if (save_into_a_file)
		delete beamline;
		delete [] energies;
		delete [] q2s;
		delete [] can_name;
		delete myRP;
//		delete tm;
		delete eff_can;
//		delete path;
	}	
} // acceptance_rp

void acceptance_aperture(bool drawlines = false, string element_name = "\"TAN.4R5\"", string filename="data/LHCB1IR5_v6.500.tfs", int side=1, int crang_sign = -1, unsigned int number_of_particles=10, float E_min=0., float E_max=150.,  unsigned int Ne = 10, float Q2_0 =-0.1, bool NonLinear=true) {
    extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {
		cout << "You should be in relative energy" << endl;
		return;
	}
	gROOT->SetStyle("Plain");
	gRandom->SetSeed(0);

	H_AbstractBeamLine * beamline;
	string rpname = "kp";
	beamline = createBeamline(420.,3000.,filename,side,rpname);

	acceptance_element(beamline,element_name,number_of_particles,E_min,E_max,Ne,Q2_0,NonLinear,crang_sign,drawlines);
	return;
}

/// Draws the 1D-acceptance for a RP detector, for two given virtuality values
void acceptance_rp_1D(float rp_pos=420., float rp_x=4000., string filename="data/LHCB1IR5_v6.500.tfs", const int side=1, const int crang_sign = -1, unsigned int number_of_particles=10, const float E_min=0., const float E_max=150., const unsigned int Ne = 10, const float Q2_0 =-0.1, const float Q2_1=-1, bool NonLinear=true, const char * add2title="", string outfilename="") {
/// @param rp_pos : distance of the roman pots, in meters from the IP
/// @param rp_x   : roman pot transverse position from the center of the beam, in micrometers [\f$\mu m\f$]
/// @param filename : name of the optics source file
/// @param side : direction of propagation (+1 : forward, -1: backward)
/// @param crang_sign : sign of the (half) crossign angle (should be -1 when forward)
/// @param number_of_particles : in the beam
/// @param E_min : min value for the beam particle energy loss
/// @param E_max : max value for  the beam particle energy loss
/// @param Ne : number of steps in energy
/// @param Q2_0 : first virtuality value
/// @param Q2_1 : second virtuality value
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
/// @param add2title : extra string added to the histogram title
/// @param outfilename : name of file to save to


	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {	
		cout << "You should be in relative energy" << endl;
		return;
	}
        gROOT->SetStyle("Plain");
	gRandom->SetSeed(0);

        H_AbstractBeamLine * beamline;
        string rp_name;
        beamline = createBeamline(rp_pos+0.1,rp_x,filename,side,rp_name);

	char can_name[100];
	sprintf(can_name,"Acceptance at %dm (%d #mum) %s",(int)rp_pos,(int)rp_x,add2title);
	TH1F * h1D_0 = acceptance_element(beamline,rp_name,number_of_particles,E_min,E_max,Ne,Q2_0,NonLinear,crang_sign);
	TH1F * h1D_1 = acceptance_element(beamline,rp_name,number_of_particles,E_min,E_max,Ne,Q2_1,NonLinear,crang_sign);


	// draws the contour graph
	TCanvas * eff_can = new TCanvas("cont2",can_name);
	THStack * histos = new THStack();
	histos->SetTitle(can_name);
	h1D_0->SetLineWidth(2);
	histos->Add(h1D_1);
	histos->Add(h1D_0);
	h1D_1->SetLineWidth(0);
	h1D_1->SetFillColor(39);
	histos->Draw("nostack");	
        histos->GetXaxis()->SetTitle("E_{loss} (GeV)");
        histos->GetYaxis()->SetTitle("Relative number");

	TLegend * leg = new TLegend(0.75,0.55,0.95,0.78);
	char legitem[100];
	sprintf(legitem,"Q^{2} = %.1f  GeV^{2}",-Q2_0);
	leg->AddEntry(h1D_0,legitem);
	sprintf(legitem,"Q^{2} = %.1f  GeV^{2}",-Q2_1);
	leg->AddEntry(h1D_1,legitem);
	leg->SetBorderSize(1);
	leg->SetFillColor(0);
	leg->Draw();

	if(outfilename.size()) { 
		eff_can->Print(outfilename.c_str(),"eps");
		delete eff_can;
		delete beamline;
		delete histos;
		delete h1D_1;
		delete h1D_0;
		delete leg;
	}	
} // acceptance_rp_1D


/// draws the aperture shape of a given optical element, as well as protons "sensing" this aperture
void acceptance_profile(string element_name = "\"MB.B9R5.B1\"", const float energy = 110, const float xmin = -26, const float xmax = -17, const float ymin = -2, const float ymax = 2, const int n = 100) {
	/// @param element_name is the id of the optical element whose aperture is being drawn

	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {	
		cout << "You should be in relative energy" << endl;
		return;
	}

	char title[200];
	sprintf(title,"Aperture effect of %s on %d GeV energy loss protons",element_name.c_str(),(int)energy);
	TH2F* pnotstopped = new TH2F("ns",title,200,xmin,xmax,200,ymin,ymax);
	TH2F* pstopped = new TH2F("ps",title,200,xmin,xmax,200,ymin,ymax);
	pstopped->SetMarkerColor(kRed);

	extern int kickers_on;
	kickers_on = 1;

        gROOT->SetStyle("Plain");
	gRandom->SetSeed(0);

	H_BeamLine* beamline = new H_BeamLine(1,500.);
        beamline->fill("data/LHCB1IR5_v6.500.tfs");

	for(int i=0;i<n;i++) {
		H_BeamParticle p1;
	        p1.setE(7000);
		p1.emitGamma(energy,0);
		p1.smearPos();
		p1.smearAng();
		p1.smearE();
		p1.setPosition(p1.getX()-500.,p1.getY(),p1.getTX()-CRANG,p1.getTY(),0.);
		p1.computePath(beamline,1);
		p1.propagate(beamline,element_name);
		if((p1.stopped(beamline))&&(p1.getStoppingElement()->getName()==element_name)) {
				pstopped->Fill(p1.getX()/1000.,p1.getY()/1000.);
				if (!i) cout << p1.getX()/1000. << "\t" << p1.getY()/1000. << endl;
			} else {
				pnotstopped->Fill(p1.getX()/1000.,p1.getY()/1000.);
				if (!i) cout << p1.getX()/1000. << "\t" << p1.getY()/1000. << endl;
			}
		}
	TCanvas* can = new TCanvas("can","can",1);
	can->SetGrid();
	pstopped->GetXaxis()->SetTitle("x [mm]");
	pstopped->GetYaxis()->SetTitle("y [mm]");
	pstopped->Draw();
	pstopped->SetStats(0);
	beamline->getElement(element_name)->getAperture()->draw(1E-3);
	pnotstopped->Draw("same");
	pnotstopped->SetStats(0);
	return;
}

int main() {
	cout << "Ca commence" << endl;
//	acceptance_rp(420,3550,"data/LHCB1IR5_v6.500.tfs",1,-1,10,0,150,10,-3,0,10,1,"robert","robo.eps");
//	acceptance_rp_1D(220.,2000.,"data/LHCB1IR5_v6.500.tfs", 1,-1,1000,0.,1400.,30,-0.1,-1.,1," for beam 1","roboo.eps");
	acceptance_fluence(220., 2000., "data/LHCB1IR5_v6.500.tfs", 1, 1000, true);
//	acceptance_basiccheck(0,1000,10000);
	cout << "C'est fini" << endl;
	return 0;
}
