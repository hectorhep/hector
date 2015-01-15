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

/// \file H_Display.cpp
/// \brief Some nice plots to show the potential of Hector

// c++ #includes
#include <iostream>
#include <string>
#include <cstdlib> //srand

// ROOT #includes
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TRandom.h"

// local #includes
#include "H_AbstractBeamLine.h"
#include "H_Beam.h"
#include "H_BeamLine.h"
#include "H_BeamParticle.h"

using namespace std;

/// \brief Allows the multiplication of a TGraph by a float; Not implemented in ROOT.
TGraph * rescale_tgraph(const TGraph * input, const float a) {
	const int N = input->GetN();
	float x[N], y[N];
	for (int i=0; i<N; i++) {
		x[i] = (input->GetX())[i];
		y[i] = a*(input->GetY())[i];
	}
	TGraph * output = new TGraph(N,x,y);
	output->SetLineColor(input->GetLineColor());
	return output;
}


TMultiGraph * beam_profile_x=0, * beam_profile_y=0;
const unsigned int Nparticles=100;
TCanvas * can=0;
TCanvas * canprof=0;
TH2F * pprof=0, * nprof=0;
extern int kickers_on;
extern bool relative_energy;
H_BeamLine* beamline=0;
const string beam1file = "data/LHCB1IR5_7TeV.tfs";
const string beam2file = "data/LHCB2IR5_7TeV.tfs";
const string ipname = "IP5";
const float length = 500;

void run(const bool relative, const bool beam1) {
	
	unsigned int stopped_number_b1=0, color=0;
	string can1name, can2name, can1title, can2title, mg1name, mg2name, th2d1name, th2d2name;
	if(relative) {
		can1name="can_relative"; can1title="laboratory frame";
		mg1name="beam_profile_x_relative"; mg2name="beam_profile_y_relative";
		th2d1name="pprof_relative"; th2d2name="nprof_relative";
	} else {
		can1name="can_absolute"; can1title="beam frame";
		mg1name="beam_profile_x_absolute"; mg2name="beam_profile_y_absolute";
		th2d1name="pprof_absolute"; th2d2name="nprof_absolute";
	}

	if(gROOT->FindObject("can_absolute")) {

		can = (TCanvas*) gROOT->FindObject("can_absolute");
		delete can;
	}
	if(relative) can = new TCanvas(can1name.c_str(),can1title.c_str(),420,10,400,400);
	else can = new TCanvas(can1name.c_str(),can1title.c_str(),10,10,400,400);
	can->cd();

	relative_energy = relative;
	string beamfilename;
	int sign;

	if(beam1) {
	// beam 1 forward
	  sign = -1;
	  beamfilename = beam1file;
	  beamline = new H_BeamLine(1,length);  // side 1: right of IP
	  beamline->fill(beamfilename,1,ipname);   // 1 is forward direction
	  beamline->offsetElements(120,-0.097); // -0.097 because forward
	}
	else {
	// beam 2 forward
	  sign = -1;
	  beamfilename = beam2file;
          beamline = new H_BeamLine(-1,length); // side -1: left of IP
          beamline->fill(beamfilename,1,ipname);   // 1 is forward direction 
          beamline->offsetElements(120,+0.097); // -0.097*-1 because forward
	}
	if(gROOT->FindObject(mg1name.c_str())) {

		beam_profile_x = (TMultiGraph*) gROOT->FindObject(mg1name.c_str());
		delete beam_profile_x;
	}
	beam_profile_x = new TMultiGraph(mg1name.c_str(),"");

	if(gROOT->FindObject(mg2name.c_str())) {

		beam_profile_y = (TMultiGraph*) gROOT->FindObject(mg2name.c_str());
	}
	beam_profile_y = new TMultiGraph(mg2name.c_str(),"");

	stopped_number_b1 =0;
	if(gROOT->FindObject(th2d1name.c_str())) {

		pprof = (TH2F*) gROOT->FindObject(th2d1name.c_str());
		delete pprof;
	}
	if (relative) pprof = new TH2F(th2d1name.c_str(),"protons",100,3,6,100,-5,5);
	else pprof = new TH2F(th2d1name.c_str(),"protons",100,-93,-89,100,-5,5);

	if(gROOT->FindObject(th2d2name.c_str())) {

		nprof = (TH2F*) gROOT->FindObject(th2d2name.c_str());
		delete nprof;
	}
	float nnn=20.;
	if (relative) nprof = new TH2F(th2d2name.c_str(),"neutrons",100,-nnn,nnn,100,-nnn,nnn);
	else nprof = new TH2F(th2d2name.c_str(),"neutrons",100,-nnn,nnn,100,-nnn,nnn);

	float q=0;

	for(unsigned int i=0;i<2*Nparticles;i++) {
		color = 1+i%2; // protons in black; neutrons in red
		q = i%2; // protons or neutrons

		H_BeamParticle p1(MP,q);

		p1.smearPos();	p1.smearAng();
		if(ipname=="IP1") 
			p1.setPosition(p1.getX(),p1.getY()+500.,p1.getTX(),p1.getTY()+sign*CRANG*(i%2),0);
		else p1.setPosition(p1.getX()-500.,p1.getY(),p1.getTX()+sign*CRANG*(i%2),p1.getTY(),0);

		p1.computePath(beamline,1);

		TGraph * ppath_x = p1.getPath(0,color);
		TGraph * ppath_y = p1.getPath(1,color);
		beam_profile_x->Add(rescale_tgraph(ppath_x,1E-3));
		beam_profile_y->Add(rescale_tgraph(ppath_y,1E-3)); // put the graphs in mm and not in um
		if(p1.stopped(beamline)) stopped_number_b1++;

	}

	for(unsigned int i=0;i<2*Nparticles;i++) {

		color = 1+i%2; // protons in black; neutrons in red
		q = i%2; // protons or neutrons

		H_BeamParticle p1(MP,q);

		p1.smearPos();	p1.smearAng();
		if(ipname=="IP1") 
			p1.setPosition(p1.getX(),p1.getY()+500.,p1.getTX(),p1.getTY()+sign*CRANG*(i%2),0);
		else p1.setPosition(p1.getX()-500.,p1.getY(),p1.getTX()+sign*CRANG*(i%2),p1.getTY(),0);

		p1.computePath(beamline,1);
		p1.propagate(length);
		if(i==1) cout << p1.getX() << " " << p1.getY() << endl;
		if(i%2) { pprof->Fill(p1.getX()/1e3,p1.getY()/1e3);}
		else    { nprof->Fill(p1.getX()/1e3,p1.getY()/1e3);}
	}
	

	can->Divide(1,2);
	can->cd(1); 
	gPad->SetGrid();
	beam_profile_x->Draw("AL+");
	beam_profile_x->GetXaxis()->SetTitle("s [m]");
	beam_profile_x->GetXaxis()->SetTitleOffset(0.5);
	beam_profile_x->GetXaxis()->SetTitleSize(0.07);
	beam_profile_x->GetYaxis()->SetTitle("x [mm]");
	beam_profile_x->GetYaxis()->SetTitleOffset(0.5);
	beam_profile_x->GetYaxis()->SetTitleSize(0.07);

	// for drawing purposes
        beamline->offsetElements(0,-0.01);
        beamline->offsetElements(120,-0.097);
	beamline->drawX(beam_profile_x->GetYaxis()->GetXmin()/5.,beam_profile_x->GetYaxis()->GetXmax()/5.,1E-3);
	beam_profile_x->Draw("L+");

	can->cd(2); beam_profile_y->Draw("AL+");
	gPad->SetGrid();
	beam_profile_y->GetXaxis()->SetTitle("s [m]");
	beam_profile_y->GetXaxis()->SetTitleOffset(0.5);
	beam_profile_y->GetXaxis()->SetTitleSize(0.07);
	beam_profile_y->GetYaxis()->SetTitle("y [mm]");
	beam_profile_y->GetYaxis()->SetTitleOffset(0.5);
	beam_profile_y->GetYaxis()->SetTitleSize(0.07);
	beam_profile_y->Draw("L+");
//	cout << endl << stopped_number_b1 << " particles have been stopped in beam 1" << endl;

	if(relative) { canprof->cd(2); }
	else { canprof->cd(1);}
	pprof->SetMarkerColor(kRed);
	pprof->Draw("P");
	if(relative) { canprof->cd(4); }
	else { canprof->cd(3);}
	nprof->SetMarkerColor(kBlack);
	nprof->Draw("P");
}


/// \brief Shows what is happening around the ZDC region
void zdc_lhcbeams(const bool beam1=true) {
	/// @param length : maximum of the s axis
	/// @param beam1file : source file for beam 1 optics
	/// @param ipname : string identifying the IP position ("IP1" or "IP5")

	kickers_on = 1;
	gROOT->SetStyle("Plain");

	if(gROOT->FindObject("can_profiles")) {

		can = (TCanvas*) gROOT->FindObject("can_profiles");
		delete can;
	}
	canprof = new TCanvas("can_profiles","",50,50,800,800);
	canprof->Divide(2,2);


	run(true,beam1);
	run(false,beam1);

} // display_lhcbeams
