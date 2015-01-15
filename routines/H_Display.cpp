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

// ROOT #includes
#include "TCanvas.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TEllipse.h"
#include "TROOT.h"

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

void display_nonip(double via = 100., float length=500., string beam1file="data/LHCB1IR5_v6.500.tfs") {

	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = true;

	H_BeamLine* beamline = new H_BeamLine(1,length);
	beamline->fill(beam1file,1,"IP5");

	TMultiGraph * beam_profile_x1 = new TMultiGraph("beam_profile_x1","");

	const int N = 10;

	for(int i=0; i < N; i++) {
		H_BeamParticle p;
		p.smearPos();
		p.smearAng();
		p.setPosition(p.getX()-500.,p.getY(),p.getTX()-CRANG,p.getTY(),0);
		p.computePath(beamline,1);
		p.propagate(via);
		TGraph * ppath_x1 = p.getPath(0,1);
		beam_profile_x1->Add(ppath_x1);

		H_BeamParticle pp;
		pp.setPosition(p.getX(),p.getY(),p.getTX(),p.getTY(),via);
		pp.computePath(beamline,1);
		TGraph * ppath_x2 = pp.getPath(0,2);
		beam_profile_x1->Add(ppath_x2);
	}

	TCanvas* can = new TCanvas("can","can",1);
	can->cd();
	gPad->SetGrid();
	beam_profile_x1->Draw("al");

	return;
}

int main() {
	display_nonip(100);
	return 0;
}

/// \brief Shows both LHC beams aside, from the top and from the side
/// <BR> Usage : display_lhcbeams(500.,"data/LHCB1IR1.tfs"       , "data/LHCB2IR1.tfs"       , "IP1")
/// <BR> Usage : display_lhcbeams(500.,"data/LHCB1IR5_v6.500.tfs", "data/LHCB2IR5_v6.500.tfs", "IP5")
void display_lhcbeams(float length=500., string beam1file="data/LHCB1IR5_v6.500.tfs", string beam2file="data/LHCB2IR5_v6.500.tfs", char * ipname = "IP5") {
	/// @param length : maximum of the s axis
	/// @param beam1file : source file for beam 1 optics
	/// @param beam2file : source file for beam 2 optics
	/// @param ipname : string identifying the IP position ("IP1" or "IP5")
	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = false;
	if(relative_energy) {
		cout << "You should better run it in absolute energy... " << endl;
	}

	gRandom->SetSeed(0);
	gROOT->SetStyle("Plain");
	H_BeamLine* beamline1 = new H_BeamLine(1,length);
	H_BeamLine* beamline2 = new H_BeamLine(1,length);
	H_AbstractBeamLine * beamline;

	beamline1->fill(beam1file,1,ipname);
	beamline2->fill(beam2file,-1,ipname);

	beamline1->offsetElements(120,-0.097);
	beamline2->offsetElements(120,+0.097);

	cout << "IP : " << beamline1->getIP() << " " << beamline2->getIP() << endl;
	cout << "============BEAM 1==============" << endl;
	beamline1->showElements();
//	beamline1->calcMatrix();
	cout << "============BEAM 2==============" << endl;
	beamline2->showElements();
//	beamline2->calcMatrix();

	TMultiGraph * beam_profile_x = new TMultiGraph("beam_profile_x","");
	TMultiGraph * beam_profile_y = new TMultiGraph("beam_profile_y","");

	int color =1;
	int stopped_number_b1 =0, stopped_number_b2 =0;

	for(int j=1;j<=2; j++)	// BEAM 1 or 2
	for(int i=0;i<100;i++) {
		if(j==1) beamline = beamline1;
		else beamline = beamline2;
		int sign = (j==1) ? 1 : -1;
		color = 1;
		
		H_BeamParticle p1;

		p1.smearPos();
		p1.smearAng();
		if(strstr(ipname,"IP1")) 
			p1.setPosition(p1.getX(),p1.getY()+500.,p1.getTX(),p1.getTY()-sign*CRANG,0);
		else p1.setPosition(p1.getX()-500.,p1.getY(),p1.getTX()-sign*CRANG,p1.getTY(),0);

		p1.computePath(beamline,1);

		TGraph * ppath_x = p1.getPath(0,color+j-1);
		TGraph * ppath_y = p1.getPath(1,color+j-1);
		

		beam_profile_x->Add(rescale_tgraph(ppath_x,1E-3));
		beam_profile_y->Add(rescale_tgraph(ppath_y,1E-3)); // put the graphs in mm and not in um

		if(p1.stopped(beamline)) {
			if (j==1) stopped_number_b1++;
			else stopped_number_b2++;
		}
	}

   	TCanvas* ccc = new TCanvas("ccc","the particle test canvas",1);
	ccc->Divide(1,2);
	ccc->cd(1); 
	gPad->SetGrid();
	beam_profile_x->Draw("AL+");
	beam_profile_x->GetXaxis()->SetTitle("s [m]");
	beam_profile_x->GetYaxis()->SetTitle("beam top view (x [mm])");
	beam_profile_x->GetYaxis()->SetTitleOffset(0.75);

	// for drawing purposes
        beamline1->offsetElements(0,-0.01);
        beamline1->offsetElements(120,-0.097);
        beamline2->offsetElements(0,+0.01);
        beamline2->offsetElements(120,+0.097);

	beamline1->drawX(beam_profile_x->GetYaxis()->GetXmin()/5.,beam_profile_x->GetYaxis()->GetXmax()/5.,1E-3);
	beamline2->drawX(beam_profile_x->GetYaxis()->GetXmin()/5.,beam_profile_x->GetYaxis()->GetXmax()/5.,1E-3);
	beamline2->draw(0.8,0.25,0.95,0.75); // legend
	beam_profile_x->Draw("L+");

	ccc->cd(2); beam_profile_y->Draw("AL+");
	gPad->SetGrid();
	beam_profile_y->GetXaxis()->SetTitle("s [m]");
	beam_profile_y->GetYaxis()->SetTitle("beam side view (y [mm])");
	beam_profile_y->GetYaxis()->SetTitleOffset(0.75);
   	beamline2->drawY(beam_profile_y->GetYaxis()->GetXmin(),beam_profile_y->GetYaxis()->GetXmax());
	beam_profile_y->Draw("L+");

	cout << endl << stopped_number_b1 << " particles have been stopped in beam 1" << endl;
	cout << endl << stopped_number_b2 << " particles have been stopped in beam 2" << endl;
	
} // display_lhcbeams

/// Creates a animation (set of gif files) for the beam shape evolution through the beamline
void display_animatedprofile(const int Nstep = 20, float length = 500., string filename="data/LHCB1IR5_v6.500.tfs", int side = 1, bool draw=true, float gkk=100, float gq2=-0.1, const int Nparticles=500, const int xmin = -100000, const int xmax = 10000, const int ymin = -3000, const int ymax = 3000, const bool save=true) {
	/// @param Nstep : number of images in the animation
	/// @param length : beamline length
	/// @param filename : optics source file
	/// @param side : direction of propagation (forward 1, backward -1)
	/// @param draw : boolean
	/// @param gkk : Energy loss of the particle
	/// @param gq2 : Virtuality (or angular shift) in \f$ GeV^2 \f$
	/// @param Nparticles : beam line content
	/// @param xmin : min value for the horizontal axis
	/// @param xmax : max value for the horizontal axis
	/// @param ymin : min value for the vertical axis
	/// @param ymax : max value for the vertical axis
	/// @param save : boolean

	gROOT->SetStyle("Plain");
	TCanvas * c1 = new TCanvas("c1","",200,10,600,400);

        // beamline concerns
        H_BeamLine* beamline = new H_BeamLine(side,length+0.1);
        beamline->fill(filename);
        int direction = (side<0)?-1:1;
        beamline->offsetElements(120,-1*direction*0.097);
	extern int kickers_on;
	kickers_on = 1;

	// particle beam concerns
        H_Beam pbeam,gbeam;
	pbeam.setPosition(-500.,0.,-CRANG,0.,0.);
	gbeam.setPosition(-500.,0.,-CRANG,0.,0.);
	pbeam.createBeamParticles(Nparticles); // no gamma emission
	gbeam.createBeamParticles(Nparticles); // gamma emission
	for( int particle_i=0; particle_i<Nparticles; particle_i++) {
		gbeam.getBeamParticle(particle_i)->emitGamma(gkk,gq2);
	}
	pbeam.computePath(beamline);
	gbeam.computePath(beamline);

	// histograms concerns
	const int xbin = 400;
	const int ybin = 50;

	// loop on all positions	
	float s=0; // position from IP in [m]
	char gif_name[50];
	for (int i=0; i<Nstep; i++) {

		s += length/Nstep;
                cout << "Position s=" << s << endl;
		sprintf(gif_name,"%.2f",s);
		
	       	TH2F * phist = new TH2F("h1",gif_name,xbin,xmin,xmax,ybin,ymin,ymax);
        	phist->SetMarkerColor(kBlack); 
	        TH2F * ghist = new TH2F("h2",gif_name,xbin,xmin,xmax,ybin,ymin,ymax);
        	ghist->SetMarkerColor(kRed);

		pbeam.propagate(s);
		gbeam.propagate(s);

		for( int particle_i=0; particle_i<Nparticles; particle_i++) {
			phist->Fill(pbeam.getBeamParticle(particle_i)->getX(),pbeam.getBeamParticle(particle_i)->getY());
			ghist->Fill(gbeam.getBeamParticle(particle_i)->getX(),gbeam.getBeamParticle(particle_i)->getY());
                }

		if(draw) { ghist->Draw(); phist->Draw("same"); }

		// saving the file
		if (save) {
			if(i<=9) sprintf(gif_name,"p00%d.gif",i);
			else if(i<=99) sprintf(gif_name,"p0%d.gif",i);
			else sprintf(gif_name,"p%d.gif",i);
			c1->SaveAs(gif_name);
			cout << "Saved into " << gif_name << endl;
		}
		delete phist;
		delete ghist;
	}
	delete c1;
} // display_animatedprofile

/// used in display_beamprofile() : not working !
void getEllipseParameters(const float * x_data, const float * y_data, const unsigned int N, float& x_width, float& y_width, float& angle) {
	// In order to fit a good ellipse on the scattered plot :
	// 1) The TH2 is copied into a TGraph, to fit it with y(x) = ax => to retrieve the angle
	// 2) Rotation of the Graph to get the RMS in X and Y
	// 3) Creation of the final ellipse, with the good widths and angle

	TCanvas * ca0 = new TCanvas;
	ca0->Divide(2,1);
	ca0->cd(1);

	TGraph * draft = new TGraph(N,x_data,y_data);
	draft->Draw("AP");
	draft->Fit("pol1","Q");
	TF1 * pol1 = draft->GetFunction("pol1");
	pol1->Draw("same");

	// gets the angle [rad]
	angle = asin(1.) - atan(pol1->GetParameter(1));
	
	double x_datarot[N], y_datarot[N];
	for (unsigned int i=0; i<N; i++) {
		x_datarot[i]=  x_data[i]*cos(angle) - y_data[i]*sin(angle);
		y_datarot[i]=  x_data[i]*sin(angle) + y_data[i]*cos(angle);
	}

	ca0->cd(2);
	TGraph * draft2 = new TGraph(N,x_datarot,y_datarot);
	draft2->Draw("AP");
	x_width =  draft2->GetRMS(1); 
	y_width =  draft2->GetRMS(2);
	angle = 180-90*angle/asin(1.);
//	draft->Draw("AP");
	ca0->cd(1);
	TEllipse * ell = new TEllipse(draft->GetMean(1),draft2->GetMean(2),x_width*3,y_width*3);
	ell->SetTheta(angle);
	ell->Draw("same");

	//cout << "x = " << x_width << "\t y = " << y_width << "\t angle = " << angle << endl;

//	delete draft2;
//	delete draft;
//	delete ca0;
	return;
}


/// plots the beam profiles in (x1,x2) and (x,x') planes
void display_beamprofile(float s, string filename="data/LHCB1IR5_v6.500.tfs", char * ipname = "IP5", int side = 1, char * title ="", unsigned int NParticle=1000, const int crang_sign=-1, const bool save=false, char * outfilename="") {
	/// @param s : distance from IP [m]
	/// @param filename : optics source file
	/// @param ipname : string identifier for the IP position
	/// @param side : direction of propagation (forward 1, backward -1)
	/// @param title : for the graph
	/// @param NParticle : beam content
	/// @param crang_sign : direction for the (half) crossing angle at IP
	/// @param save : boolean
	/// @param outfilename: file to be written


// note : beam 1 forward  : side = 1  crang_sign =-1
// note : beam 1 backward : side = -1 crang_sign = 1
// note : beam 2 forward  : side = -1 crang_sign =-1
// note : beam 2 backward : side = 1  crang_sign = 1

	extern bool relative_energy;
	relative_energy = false;
	if(relative_energy) {
		cout << "You should be in absolute energy" << endl;
		return;
	}

        extern int kickers_on;
        kickers_on = 1;

	
	int max = (crang_sign>0)?100:-95;
	int min = (crang_sign<0)?-100:95;
	TH2F * hp  = new TH2F("Positions","",100,min,max,100,-2.5,2.5);
	TH2F * ha  = new TH2F("Angles","",100,-50,50,100,-50,50);
	TH2F * hax  = new TH2F("Phase_x","",100,min,max,100,-50,50);
	TH2F * hay  = new TH2F("Phase_y","",100,-2.5,2.5,100,-50,50);
	float draftx[NParticle], drafty[NParticle], drafttx[NParticle], draftty[NParticle];
//	float rmsx=0, rmsy=0, angle=0;
	TMultiGraph * profile = new TMultiGraph("prof","");

        H_BeamLine* beamline = new H_BeamLine(side,s+0.1);
        beamline->fill(filename,-1*side*crang_sign,ipname);
	beamline->offsetElements(120,0.097*crang_sign);

//	extern int kickers_on;
//	kickers_on = 1;

	for (unsigned int i=0; i<NParticle ; i++) {
		H_BeamParticle p1;
		p1.smearPos();
		p1.smearAng();
		p1.setPosition(p1.getX()-500.,p1.getY(),p1.getTX()+crang_sign*CRANG,p1.getTY(),0);
		p1.computePath(beamline);
		p1.propagate(beamline);
		p1.propagate(s);
		hp->Fill(p1.getX()/1000.,p1.getY()/1000.);
		ha->Fill(p1.getTX(),p1.getTY());
		hax->Fill(p1.getX()/1000.,p1.getTX());
		hay->Fill(p1.getY()/1000.,p1.getTY());
		draftx[i]=p1.getX()/1000.;
		drafty[i]=p1.getY()/1000.;
		drafttx[i]=p1.getTX();
		draftty[i]=p1.getTY();
		TGraph * path = p1.getPath(0,1);
		profile->Add(path);
	}

	TCanvas * can = new TCanvas;
	can->cd();
	hp->SetTitle(title);
	hp->Draw();	
	hp->GetXaxis()->SetTitle("x (mm)");
	hp->GetYaxis()->SetTitleOffset(1.2);
	hp->GetYaxis()->SetTitle("y (mm)");
	TEllipse * ellipse = new TEllipse(hp->GetMean(1),hp->GetMean(2),3*(hp->GetRMS(1)),3*(hp->GetRMS(2)));
	cout << "mean = " << hp->GetMean(1) << " " << hp->GetMean(2) << endl;
	ellipse->SetLineColor(kRed);
	ellipse->Draw();
	
	TCanvas * ca2 = new TCanvas;
	ca2->cd();
	profile->Draw("ACP");

	TCanvas *ca3 = new TCanvas;
	ca3->cd();
	ha->SetTitle(title);
	ha->Draw();
	ha->GetXaxis()->SetTitle("#theta_{x} (#murad)");
	ha->GetYaxis()->SetTitle("#theta_{y} (#murad)");
	TEllipse * ellips2 = new TEllipse(ha->GetMean(1),ha->GetMean(2),3*(ha->GetRMS(1)),3*(ha->GetRMS(2)));
	ellips2->SetLineColor(kRed);
	ellips2->Draw();

	TCanvas *ca4 = new TCanvas;
	ca4->cd();
	hax->SetTitle(title);
	hax->Draw();
	hax->SetStats(0);
	hax->GetXaxis()->SetTitle("x (mm)");
	hax->GetYaxis()->SetTitle("#theta_{x} (#murad)");
//	getEllipseParameters(draftx,drafttx,NParticle,rmsx,rmsy,angle);
//	ca4->cd();
//	cout << rmsx << " " << rmsy << " " << angle << endl;
//	TEllipse * ellips3 = new TEllipse(hp->GetMean(1),ha->GetMean(1),3*rmsx,3*rmsy);
//	ellips3->SetTheta(angle);
//	ellips3->SetLineColor(kRed);
//	ellips3->Draw();

	TCanvas *ca5 = new TCanvas;
	ca5->cd();
	hay->SetTitle(title);
	hay->Draw();
	hay->SetStats(0);
	hay->GetXaxis()->SetTitle("y (mm)");
	hay->GetYaxis()->SetTitle("#theta_{y} (#murad)");
//	getEllipseParameters(drafty,draftty,NParticle,rmsx,rmsy,angle);
//	ca5->cd();
//	cout << rmsx << " " << rmsy << " " << angle << endl;
//	TEllipse * ellips4 = new TEllipse(hp->GetMean(2),ha->GetMean(2),3*rmsx,3*rmsy);
//	ellips4->SetTheta(angle);
//	ellips4->SetLineColor(kRed);
//	ellips4->Draw();

	if(save) {
		char filetitle_pos[50], filetitle_ang[50], filetitle_phasex[50], filetitle_phasey[50];
		sprintf(filetitle_pos,"%s_pos.eps",outfilename);
		cout << filetitle_pos << endl;
		can->Print(filetitle_pos,"eps");
		sprintf(filetitle_ang,"%s_ang.eps",outfilename);
		cout << filetitle_ang << endl;
		ca3->Print(filetitle_ang,"eps");
		sprintf(filetitle_phasex,"%s_px.eps",outfilename);
		cout << filetitle_phasex << endl;
		ca4->Print(filetitle_phasex,"eps");
		sprintf(filetitle_phasey,"%s_py.eps",outfilename);
		cout << filetitle_phasey << endl;
		ca5->Print(filetitle_phasey,"eps");
		delete can;
		delete ca2;
		delete ca3;
		delete ca4;
		delete ca5;
		delete hp;
		delete ha;
		delete hax;
		delete hay;
		delete profile;
		delete beamline;
		delete ellipse;
		delete ellips2;
//		delete ellips3;
//		delete ellips4;
	}

}
