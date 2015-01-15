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

/// \file H_Reconstruction.cpp
/// \brief Collection of routines about explicit physical variable reconstruction (see also H_Acceptance.cpp)

// c++ #includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

// ROOT #includes
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMarker.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TRandom.h"

// local #includes
#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_RecRPObject.h"
using namespace std;

/// 2D- chromaticity grid in (x,x') and (x1,x2) planes
void reconstruction_drawgrid(double pos1 = 220., double pos2 = 224., int nee = 11,  double emin = 100., double emax = 1000., int nxp = 11, double xpmin = -500., double xpmax = 500., string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1, string outfilename=" ") {
	/// @param pos1 = 220. is the s position of RP1, in [m] from IP;
	/// @param pos2 = 224. is the s position of RP2, in [m] from IP;
	/// @param emin = 100. is the minimal photon energy in the grid
	/// @param emax = 1000. is the maximal photon energy in the grid
	/// @param xpmin = 0. is the minimal x-prime in the grid
	/// @param xpmax = 0. is the maximal x-prime in the grid
	/// @param filename is the default optics file
	/// @param side is the direction of the beam (1 is clockwise, -1 is counterclockwise)
	/// @param outfilename : file to be written -- do not save if empty
	/// Q2 is always minimal here => no phi angle dependence.
	// !!! Q2 is always minimal : p1.emitGamma(energy,0);
	
	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {
        cout << "You should be in relative energy" << endl;
        return;
	}
	
	H_BeamLine* beam = new H_BeamLine(side,pos2+5);
	beam->fill(filename);

	char outfile_name[40];
	sprintf(outfile_name,"chromaticity_grid_%d.txt",(int)pos1);
	ofstream outfile(outfile_name);
	outfile<<"#------------------------------------------------------------------------------"<<endl;
	outfile<<"# This is HECTOR chromaticity grid @ "<<(int)pos1<<" m"<<endl;
	outfile<<"# Energy runs (vertically) from "<<(int)emin<<" to "<<(int)emax<<" Gev with "<<nee<<" values."<<endl;
	outfile<<"# X-Angle runs (horizontally) from "<<(int)xpmin<<" to "<<(int)xpmax<<" µrad with "<<nxp<<" values."<<endl;
	outfile<<"# Format : x(µm),theta_x(µrad)"<<endl;
	outfile<<"#------------------------------------------------------------------------------"<<endl;
	outfile<<"# please check the optics version is correct : "<<filename<<endl; 
	outfile<<"#------------------------------------------------------------------------------"<<endl;
	outfile<<"# Contact persons  : Xavier Rouby, Jerome de Favereau, Krzysztof Piotrzkowski  "<<endl;
	outfile<<"# rouby@fynu.ucl.ac.be, favereau@fynu.ucl.ac.be, k.piotrzkowski@fynu.ucl.ac.be "<<endl;
	outfile<<"#------------------------------------------------------------------------------"<<endl;
	outfile<<"# variables : "<<endl;
	outfile.width(7);
	outfile.fill(' ');
	outfile<<left<<(int)pos1<<"\t # detector distance from IP"<<endl;
	outfile.width(7);
	outfile.fill(' ');
	outfile<<left<<nee<<"\t # energy bins"<<endl;
	outfile.width(7);
	outfile.fill(' ');
	outfile<<left<<(int)emin<<"\t # first energy"<<endl;
	outfile.width(7);
	outfile.fill(' ');
	outfile<<left<<(int)emax<<"\t # last energy "<<endl;
	outfile.width(7);
	outfile.fill(' ');
	outfile<<left<<nxp<<"\t # x-angle bins"<<endl;
	outfile.width(7);
	outfile.fill(' ');
	outfile<<left<<(int)xpmin<<"\t # first angle"<<endl;
	outfile.width(7);
	outfile.fill(' ');
	outfile<<left<<(int)xpmax<<"\t # last angle "<<endl;
	outfile<<"#------------------------------------------------------------------------------"<<endl;

	const double milli = 0.001;
	const double mega = 1000000.;

	double x1rp[nee*nxp], x2rp[nee*nxp], t1rp[nee*nxp];
	double theee[nee*nxp], thexp[nee*nxp];
	double energy, xp;

	for(int i = 0; i < nee; i++) {
		for(int j = 0; j < nxp; j++) {
			energy = emin + i*(emax-emin)/((double)nee-1);
			xp = xpmin + j*(xpmax-xpmin)/((double)nxp-1);
			theee[nee*i+j] = energy;
			thexp[nee*i+j] = xp;
			H_BeamParticle p1;
			p1.setPosition(PX,PY,xp-CRANG,0,0);
//			p1.emitGamma(energy,0);
			p1.setE(BE-energy);
			p1.computePath(beam,1);
			p1.propagate(pos1);
			x1rp[nee*i+j] = p1.getX()*milli;
			p1.propagate(pos2);
			x2rp[nee*i+j] = p1.getX()*milli;
			t1rp[nee*i+j] = atan(milli*(x2rp[nee*i+j] - x1rp[nee*i+j])/(double)(pos2-pos1))*mega;
			outfile.width(7);
			outfile.fill(' ');
			outfile<<right<<fixed<<setprecision(2)<<x2rp[nee*i+j]<<","<<t1rp[nee*i+j]<<"\t";
		}
		outfile<<endl;
	}


	// x1,x2 for iso-E
	double x1e[nee], x2e[nee], t1e[nee];

	TLegend* leg = new TLegend(0.65,0.15,0.88,0.35);
	leg->SetFillColor(0);
	leg->SetHeader("(Energy, #theta)");

	TLegend* le2 = new TLegend(0.65,0.15,0.88,0.35);
	le2->SetFillColor(0);
	le2->SetHeader("(Energy, #theta)");

	char multit[500];
	sprintf(multit,"Chromaticity grid at %d m",(int)pos1);

	TMultiGraph* multi = new TMultiGraph("multi",multit);
	TMultiGraph* mult2 = new TMultiGraph("mult2",multit);

	for(int i = 0; i < nee; i++) {
		for(int j = 0; j < nxp; j++) {
			x1e[j] = x1rp[nee*i+j];
			t1e[j] = t1rp[nee*i+j];
			x2e[j] = x2rp[nee*i+j];
		}
		TGraph* ioee = new TGraph(nxp,x1e,t1e);
		ioee->SetLineColor(2);

		TGraph* joee = new TGraph(nxp,x1e,x2e);
		joee->SetLineColor(2);
		if(!i) { 
			leg->AddEntry(ioee,"Iso-energy (GeV)", "pl"); 
			le2->AddEntry(joee,"Iso-energy (GeV)", "pl");
		}
		multi->Add(ioee);
		mult2->Add(joee);
	}

	// x1,x2 for iso-xp
	double x1x[nxp], x2x[nxp], t1x[nxp];
	         
	for(int i = 0; i < nxp; i++) {
		for(int j = 0; j < nee; j++) {
			x1x[j] = x1rp[i+nxp*j];
			t1x[j] = t1rp[i+nxp*j];
			x2x[j] = x2rp[i+nxp*j];
		}
		TGraph* ioxp = new TGraph(nee,x1x,t1x);
		ioxp->SetLineColor(4);
		ioxp->SetLineStyle(2);

		TGraph* joxp = new TGraph(nee,x1x,x2x);
		joxp->SetLineColor(4);
		joxp->SetLineStyle(2);
		if(!i) {
			leg->AddEntry(ioxp,"Iso-angle (#murad)", "pl");
			le2->AddEntry(joxp,"Iso-angle (#murad)", "pl");
		}
		multi->Add(ioxp);
		mult2->Add(joxp);
	}

	float x1,x2,y1,y2,width,height,xmin,xmax,ymin,ymax;

	TCanvas * c2 = new TCanvas();
	c2->cd();
	mult2->Draw("AL");
	mult2->GetXaxis()->SetTitle("x_{1} (mm)");
	mult2->GetYaxis()->SetTitle("x_{2} (mm)");
	le2->Draw();
	gPad->SetGrid();

 	xmin = mult2->GetXaxis()->GetXmin();
	xmax = mult2->GetXaxis()->GetXmax();
        ymin = mult2->GetYaxis()->GetXmin();
	ymax = mult2->GetYaxis()->GetXmax();
	width = (xmax-xmin)/8.;
	height = (ymax-ymin)/15.;
	
	int ind2[4] = {0,nee-1,nee*(nxp-1),nee*nxp-1};
	for(int i = 0; i<4; i++) {
		x1 = x1rp[ind2[i]]-width/2.;
		x2 = x1rp[ind2[i]]+width/2.;
		y1 = x2rp[ind2[i]]-height;
		y2 = x2rp[ind2[i]]-height/5.;
		if(x1<xmin) { x1 = xmin + width/10.; x2 = x1 + width; }
		if(x2>xmax) { x2 = xmax - width/10.; x1 = x2 - width; }
		if(y1<ymin) { y2 = x2rp[ind2[i]] + height; y1 = x2rp[ind2[i]] + height/5.; }
		TPaveText* lowlo2 = new TPaveText(x1,y1,x2,y2);
		char lltext[50];
		sprintf(lltext,"(%d,%d)",(int)theee[ind2[i]],(int)thexp[ind2[i]]);
		lowlo2->AddText(lltext);
		lowlo2->SetBorderSize(0);
		lowlo2->SetFillColor(0);
		lowlo2->SetFillStyle(1001);
		lowlo2->Draw();
		TMarker* mark2 = new TMarker(x1rp[ind2[i]],x2rp[ind2[i]],8);
		mark2->Draw();
	}
	
	TCanvas * c1 = new TCanvas();
	c1->cd();
	multi->Draw("AL");
	multi->GetXaxis()->SetTitle("x_{1} (mm)");
	multi->GetYaxis()->SetTitle("#theta_{x1} (#murad)");
	leg->Draw();
	gPad->SetGrid();
	
 	xmin = multi->GetXaxis()->GetXmin();
	xmax = multi->GetXaxis()->GetXmax();
        ymin = multi->GetYaxis()->GetXmin();
	ymax = multi->GetYaxis()->GetXmax();
	width = (xmax-xmin)/8.;
	height = (ymax-ymin)/15.;
	
	int ind[4] = {0,nee-1,nee*(nxp-1),nee*nxp-1};
	for(int i = 0; i<4; i++) {
		x1 = x1rp[ind[i]]-width/2.;
		x2 = x1rp[ind[i]]+width/2.;
		y1 = t1rp[ind[i]]-height;
		y2 = t1rp[ind[i]]-height/5.;
		if(x1<xmin) { x1 = xmin + width/10.; x2 = x1 + width; }
		if(x2>xmax) { x2 = xmax - width/10.; x1 = x2 - width; }
		if(y1<ymin) { y2 = t1rp[ind[i]] + height; y1 = t1rp[ind[i]] + height/5.; }
		TPaveText* lowlow = new TPaveText(x1,y1,x2,y2);
		char lltext[50];
		sprintf(lltext,"(%d,%d)",(int)theee[ind[i]],(int)thexp[ind[i]]);
		lowlow->AddText(lltext);
		lowlow->SetBorderSize(0);
		lowlow->SetFillColor(0);
		lowlow->SetFillStyle(1001);
		lowlow->Draw();
		TMarker* mark = new TMarker(x1rp[ind[i]],t1rp[ind[i]],8);
		mark->Draw();
	}

	if (! outfilename.size() ) {
		char output[100];
		sprintf(output,"%s_x1VSx2.eps",outfilename.c_str());
		c2->Print(output);
		sprintf(output,"%s_x1VSxp1.eps",outfilename.c_str());
		c1->Print(output);

		delete c1;
		delete c2;
		delete leg;
		delete le2;
		delete beam;
		delete multi;
		delete mult2;
	}
	
	return;
}

/// Chromaticity plot for vertical variable
void reconstruction_drawgrid_y(double pos1 = 220., double pos2 = 224., double ypmin = 0, double ypmax = 500., char * filename = "data/LHCB1IR5_v6.500.tfs", int side = 1) {
	/// @param pos1 = 220. is the s position of RP1, in [m] from IP;
	/// @param pos2 = 224. is the s position of RP2, in [m] from IP;
	/// @param ypmin = 0. is the minimal y-prime in the grid
	/// @param ypmax = 0. is the maximal y-prime in the grid
	/// @param filename is the default optics file
	/// @param side is the direction of the beam (1 is clockwise, -1 is counterclockwise)
	/// Q2 is always minimal here => no phi angle dependence.
	// !!! Q2 is always minimal : p1.emitGamma(energy,0);

        extern int kickers_on;
        kickers_on = 1;

	extern bool relative_energy;
	if(!relative_energy) {	
		cout << "You should be in relative energy" << endl;
		return;
	}
	
	H_BeamLine* beam = new H_BeamLine(side,pos2+5);
	beam->fill(filename);

	const int nthy = 10;
	double thy[nthy], y1[nthy], y2[nthy];

	for(int i = 0; i < nthy; i++) {
		thy[i] = ypmin + i*(ypmax - ypmin)/((double)nthy-1);
		H_BeamParticle p1;
		//p1.setPosition(0,0,0,thy[i],0);
		p1.setPosition(PX,PY,TX-CRANG,thy[i],0);
		p1.emitGamma(100,0);
		p1.computePath(beam,1);
		p1.propagate(pos1);
		y1[i] = p1.getY()/1000.;
		p1.propagate(pos2);
		y2[i] = p1.getY()/1000.;
	}
	
	TGraph* gr = new TGraph(nthy,y1,y2);
	char tit[100];
	sprintf(tit,"Effect of the vertical angle @ IP");
	gr->SetTitle(tit);
	gr->Draw("ALP*");
	sprintf(tit,"y_{1} (mm)");
	gr->GetXaxis()->SetTitle(tit);
	sprintf(tit,"y_{2} (mm)");
	gr->GetYaxis()->SetTitle(tit);
	gPad->SetGrid();

	TMarker* mark1 = new TMarker(y1[0],y2[0],8);
	mark1->Draw();
	TMarker* mark2 = new TMarker(y1[nthy-1],y2[nthy-1],8);
	mark2->Draw();

	TPaveText* low = new TPaveText(y1[0]+0.2,y2[0]+0.1, y1[0]+0.8,y2[0]+0.4); 
	sprintf(tit,"#theta_{y}^{IP} = %d #murad",(int)thy[0]);
	low->AddText(tit);
	low->SetBorderSize(0);
	low->SetFillColor(0);
	low->SetFillStyle(1001);
	low->SetTextSize(0.03);
	low->Draw();

	TPaveText* lo2 = new TPaveText(y1[nthy-1]+0.2,y2[nthy-1]+0.1, y1[nthy-1]+0.8,y2[nthy-1]+0.4); 
	sprintf(tit,"#theta_{y}^{IP} = %d #murad",(int)thy[nthy-1]);
	lo2->AddText(tit);
	lo2->SetBorderSize(0);
	lo2->SetFillColor(0);
	lo2->SetFillStyle(1001);
	lo2->SetTextSize(0.03);
	lo2->Draw();

	return;
}	

/// Energy reconstruction
void reconstruction_energy(unsigned int method = 1, float q2 = 0, float emin =0, float emax = 200, const int e_n = 5, const int p_n = 10, const float pos1=420., const float pos2=428., string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1, const bool save=false) {
	/// @param method : used for reconstruction (1=TM, 2=ACM)
	/// @param q2 : virtuality or momentum transfert (GeV^2) 
	/// @param emin : minimal energy loss (GeV)
	/// @param emax : maximal energy loss (GeV)
	/// @param e_n  : number of energy steps
	/// @param p_n  : number of particles to propagate, per step
	/// @param pos1 : RP position 1
	/// @param pos2 : RP position 2 (pos1 < pos2);
	/// @param filename : beam optics file
	/// @param side : direction of propagation

	float ee, er; // real and reconstructed energies

	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = true;

	H_BeamLine* beam = new H_BeamLine(side,pos2+5);
	beam->fill(filename);
	H_RecRPObject r(pos1,pos2,beam);
	float x1s, x2s, y1s, y2s;

	emin +=0.00000001; // if emin=0 => problem with the profile drawing

	TProfile * he = new TProfile("he","Energy reconstruction",e_n,emin,emax*1.1,emin,emax*1.1);
	TProfile * hde = new TProfile("hde","Energy resolution",e_n,emin,emax,-5.,5.);
	for(int e_i = 0; e_i < e_n; e_i++) {
		for (int p_i = 0; p_i < p_n; p_i++) {		
			ee = emin + (e_i+0.5)*(emax-emin)/((float)e_n);
			H_BeamParticle p1;
			p1.setPosition(PX,PY,-CRANG,0,0);
			p1.smearPos();
			p1.smearAng();
			//p1.smearE();
			p1.emitGamma(ee,q2);
			// should be non linear when using r.getE
			p1.computePath(beam,1);
			p1.propagate(pos1);
			x1s = p1.getX();
			y1s = p1.getY();
			p1.propagate(pos2);
			x2s = p1.getX();
			y2s = p1.getY();
			r.setPositions(x1s,y1s,x2s,y2s);
			er = r.getE(method);
			he->Fill(ee,er);
			hde->Fill(ee,er-ee);
			//if(!(e_i+p_i)) p1.getPath(0,kBlack)->Draw("ALP*");
			if (p_i == p_n-1) cout << "E = " << ee << " \tE_rec = " << er << " \t dE = " << er-ee <<  "\t dE/E = " << (er-ee)/ee << endl;
		}
	}

	char gtitle[100];
	if (method==2) sprintf(gtitle,"ACM");
	else if(method==1) sprintf(gtitle,"TM");
	else sprintf(gtitle," ");

	TCanvas* ene_can = new TCanvas("ene_can",gtitle,1);
	ene_can->cd();
	he->Draw();
	gPad->SetGrid();

	char hetitle[150];
	sprintf(hetitle,"Energy reconstruction - %s",gtitle);
	he->SetTitle(hetitle);
	he->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	he->GetYaxis()->SetTitle("Reconstructed energy loss E_{rec} (GeV)");
	he->GetXaxis()->SetRangeUser(emin,emax);
	he->GetYaxis()->SetRangeUser(emin*0.9,1.1*emax);
	he->GetYaxis()->SetLabelSize(0.02);
	he->GetXaxis()->SetLabelSize(0.02);
	he->SetStats(0);
	
	TCanvas* ene_ca2 = new TCanvas("ene_ca2",gtitle,1);
	ene_ca2->cd();
	hde->Draw();
	gPad->SetGrid();

	sprintf(hetitle,"Energy resolution @ %.0fm - %s",pos1,gtitle);
	hde->SetTitle(hetitle);
	hde->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	hde->GetYaxis()->SetTitle("E_{rec}-E_{rel} (GeV)");
	hde->GetXaxis()->SetRangeUser(emin,emax);
	hde->GetXaxis()->SetLabelSize(0.02);
	hde->GetYaxis()->SetLabelSize(0.02);
	hde->SetStats(0);

	if (save) {
		char outfilename[50],tmethod[5]="",linearity[5]="";
		if(method==TM) sprintf(tmethod,"_TM");
		else if(method==AM) sprintf(tmethod,"_AM");

		sprintf(outfilename,"energy_reconstruction%.0f%s%s.eps",pos1,tmethod,linearity);
		ene_ca2->Print(outfilename);
		
//		delete he;
		delete hde;
		delete ene_can;
		delete ene_ca2;
//		delete beam;
	}
	return;
}

/// Position reconstruction
void reconstruction_position(float pos = 100, float q2 = 0, float emin =20, float emax = 120, const int e_n = 12, const int p_n = 20, const float pos1=420., const float pos2=428., string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1, const bool save=false) {
	/// @param pos : particle offset at IP, to be reconstructed 
	/// @param q2 : virtuality or momentum transfert (GeV^2) 
	/// @param emin : minimal energy loss (GeV)
	/// @param emax : maximal energy loss (GeV)
	/// @param e_n  : number of energy steps
	/// @param p_n  : number of particles to propagate, per step
	/// @param pos1 : RP position1
	/// @param pos2 : RP position2 (pos1 < pos2)
	/// @param filename : beam optics file
	/// @param side : direction of propagation

	float ee, ar, ay; // real energies and reconstructed angles
//	pos1 , pos2  : roman pot positions
//      pos : lateral offset at IP

	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = true;
	
	H_BeamLine* beam = new H_BeamLine(side,pos2+5);
	beam->fill(filename);
	H_BeamLine* bea2 = new H_BeamLine(side,pos2+5);
	bea2->fill(filename);
	H_RecRPObject rx(pos1,pos2,beam), ry(pos1,pos2,bea2);	
	float x1s, x2s, y1s, y2s;
	float x1t, x2t, y1t, y2t;

	emin +=0.00000001; // if emin=0 => problem with the profile drawing

	
	//TProfile * he = new TProfile("he","Horizontal position reconstruction",e_n,emin,emax,-2*pos,2*pos);
	TProfile * hde = new TProfile("hde","Horizontal position resolution",e_n,emin,emax,-250.,250.);
	//TH2F * hde = new TH2F("hde","Horizontal position resolution",e_n,emin,emax,100,-120,0);
	//TProfile * hy = new TProfile("hy","Vertical position reconstruction",e_n,emin,emax,-2*pos,2*pos);
	TProfile * hdy = new TProfile("hdy","Vertical position resolution",e_n,emin,emax,-250.,250.);
	//TH2F * hdy = new TH2F("hdy","Vertical position resolution",e_n,emin,emax,100,-1,0);
	for(int e_i = 0; e_i < e_n; e_i++) {
		for (int p_i = 0; p_i < p_n; p_i++) {		
			ee = emin + (e_i+0.5)*(emax-emin)/((float)e_n);
			// horizontal position reconstruction
			H_BeamParticle px;
			px.setPosition(PX+pos,PY,-CRANG,0,0);
			// px.smearPos();
			px.smearAng();
			px.smearE();
			px.emitGamma(ee,q2);
			// should be non linear when using r.getE
			px.computePath(beam,1);
			px.propagate(pos1);
			x1s = px.getX();
			y1s = px.getY();
			px.propagate(pos2);
			x2s = px.getX();
			y2s = px.getY();
			//cout << "horizontal: " << x1s << " " << x2s << endl;
			rx.setPositions(x1s,y1s,x2s,y2s);
			rx.getE(TM);
			ar = rx.getX1();
//			he->Fill(ee,ar);
			hde->Fill(ee,ar-pos);
			if (p_i==p_n-1) cout << "dx = " << ar-pos << endl;


			// vertical position reconstruction
			H_BeamParticle py;
			py.setPosition(PX,PY+pos,-CRANG,0,0);
//			py.smearPos();
			py.smearAng();
			py.smearE();
			py.emitGamma(ee,q2);
			// should be non linear when using r.getE
			py.computePath(bea2,1);
			py.propagate(pos1);
			x1t = py.getX();
			y1t = py.getY();
			py.propagate(pos2);
			x2t = py.getX();
			y2t = py.getY();
			//cout << "vertical : " << x1t << " " << x2t << endl;
			ry.setPositions(x1t,y1t,x2t,y2t);
			ry.getE(TM);
			ay = ry.getY1(); 
//			hy->Fill(ee,ay);
			hdy->Fill(ee,ay-pos);
		}
	}

	char hetitle[150];
	TCanvas* xang = new TCanvas("xang","",1);
/*	xang->Divide(2,1);
	xang->cd(1);
	he->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Horizontal position reconstruction @ %.0fm (%.0f #mum)",pos1,pos);
	he->SetTitle(hetitle);
	he->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	he->GetYaxis()->SetTitle("Reconstructed  position #x_{rec} (#mum)");
	he->GetXaxis()->SetRangeUser(emin,emax);
	he->GetYaxis()->SetTitleOffset(1.3);
//	he->GetYaxis()->SetRangeUser(emin*0.9,1.1*emax);
//	he->GetYaxis()->SetLabelSize(0.02);
//	he->GetXaxis()->SetLabelSize(0.02);
	he->SetStats(0);
	
	xang->cd(2);
*/	hde->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Horizontal position reconstruction @ %.0fm (%.0f #mum)",pos1,pos);
	hde->SetTitle(hetitle);
	hde->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	hde->GetYaxis()->SetTitle("x_{rec}-x_{rel} (#mum)");
	hde->GetXaxis()->SetRangeUser(emin,emax);
	hde->GetYaxis()->SetTitleOffset(1.2);
//	hde->GetXaxis()->SetLabelSize(0.02);
//	hde->GetYaxis()->SetLabelSize(0.02);
	hde->SetStats(0);


	TCanvas* yang = new TCanvas("yang","",1);
/*	yang->Divide(2,1);
	yang->cd(1);
	hy->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Vertical position reconstruction @ %.0fm (%.0f #mum)",pos1,pos);
	hy->SetTitle(hetitle);
	hy->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	hy->GetYaxis()->SetTitle("Reconstructed position #y_{rec} (#mum)");
	hy->GetXaxis()->SetRangeUser(emin,emax);
	hy->GetYaxis()->SetTitleOffset(1.3);
//	hy->GetYaxis()->SetRangeUser(emin*0.9,1.1*emax);
//	hy->GetYaxis()->SetLabelSize(0.02);
//	hy->GetXaxis()->SetLabelSize(0.02);
	hy->SetStats(0);
	
	yang->cd(2);
*/	hdy->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Vertical position reconstruction @ %.0fm (%.0f #mum)",pos1,pos);
	hdy->SetTitle(hetitle);
	hdy->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	hdy->GetYaxis()->SetTitle("y_{rec}-y_{rel} (#mum)");
	hdy->GetXaxis()->SetRangeUser(emin,emax);
	hdy->GetYaxis()->SetTitleOffset(1.2);
//	hdy->GetXaxis()->SetLabelSize(0.02);
//	hdy->GetYaxis()->SetLabelSize(0.02);
	hdy->SetStats(0);

	if (save) {
		char outfilename[50];

		sprintf(outfilename,"pos_x_reconstruction_%.0f_%.0f.eps",pos,pos1);
		xang->Print(outfilename);
		sprintf(outfilename,"pos_y_reconstruction_%.0f_%.0f.eps",pos,pos1);
		yang->Print(outfilename);
		
//		delete he;
//		delete hy;
		delete hde;
		delete hdy;
		delete xang;
		delete yang;
//		delete beam;
//		delete bea2;
	}

	return;
}

/// Angle reconstruction
void reconstruction_angle(float angl = 100, float q2 = 0, float emin =20, float emax = 120, const int e_n = 12, const int p_n = 20, const float pos1=420., const float pos2=428., string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1, const bool save=false) {
	/// @param angl : particle angle at IP, to be reconstructed 
	/// @param q2 : virtuality or momentum transfert (GeV^2) 
	/// @param emin : minimal energy loss (GeV)
	/// @param emax : maximal energy loss (GeV)
	/// @param e_n  : number of energy steps
	/// @param p_n  : number of particles to propagate, per step
	/// @param pos1 : RP position1
	/// @param pos2 : RP position2 (pos1 < pos2)
	/// @param filename : beam optics file
	/// @param side : direction of propagation
	float ee, ar, ay; // real energies and reconstructed angles
//	pos1 , pos2  : roman pot positions

	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = true;
	
	H_BeamLine* beam = new H_BeamLine(side,pos2+5); // for horizontal reconstruction
	beam->fill(filename);
	H_BeamLine* bea2 = new H_BeamLine(side,pos2+5); // for vertical reconstruction
	bea2->fill(filename);
	H_RecRPObject rx(pos1,pos2,beam), ry(pos1,pos2,bea2);	
	float x1s, x2s, y1s, y2s;
	float x1t, x2t, y1t, y2t;

	emin +=0.00000001; // if emin=0 => problem with the profile drawing
	
//	TProfile * he = new TProfile("he","Horizontal angle reconstruction",e_n,emin,emax,-2*angl,2*angl);
	TProfile * hde = new TProfile("hde","Horizontal angle resolution",e_n,emin,emax,-250.,250.);
//	TProfile * hy = new TProfile("hy","Vertical angle reconstruction",e_n,emin,emax,-2*angl,2*angl);
	TProfile * hdy = new TProfile("hdy","Vertical angle resolution",e_n,emin,emax,-250.,250.);
	for(int e_i = 0; e_i < e_n; e_i++) {
		for (int p_i = 0; p_i < p_n; p_i++) {		
			ee = emin + (e_i+0.5)*(emax-emin)/((float)e_n);
			// horizontal angle reconstruction
			H_BeamParticle px;
			px.setPosition(PX,PY,angl-CRANG,0,0);
			px.smearPos();
			// px.smearAng();
			px.smearE();
			px.emitGamma(ee,q2);
			// should be non linear when using r.getE
			px.computePath(beam,1);
			px.propagate(pos1);
			x1s = px.getX();
			y1s = px.getY();
			px.propagate(pos2);
			x2s = px.getX();
			y2s = px.getY();
			rx.setPositions(x1s,y1s,x2s,y2s);
			rx.getE(TM);
			ar = rx.getTX();
//			he->Fill(ee,ar);
			hde->Fill(ee,ar-angl);
			if (p_i==p_n-1) cout << "dtheta = " << ar-angl << endl;


			// vertical angle reconstruction
			H_BeamParticle py;
			py.setPosition(PX,PY,-CRANG,angl,0);
			py.smearPos();
			// py.smearAng();
			py.smearE();
			py.emitGamma(ee,q2);
			// should be non linear when using r.getE
			py.computePath(bea2,1);
			py.propagate(pos1);
			x1t = py.getX();
			y1t = py.getY();
			py.propagate(pos2);
			x2t = py.getX();
			y2t = py.getY();
			ry.setPositions(x1t,y1t,x2t,y2t);
			ry.getE(TM);
			ay = ry.getTY();
//			hy->Fill(ee,ay);
			hdy->Fill(ee,ay-angl);
		}
	}

	char hetitle[150];
	TCanvas* xang = new TCanvas("xang","",1);
/*	xang->Divide(2,1);
	xang->cd(1);
	he->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Horizontal angle reconstruction @ %.0fm (%.0f #murad)",pos1,angl);
	he->SetTitle(hetitle);
	he->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	he->GetYaxis()->SetTitle("Reconstructed angle #theta^{x}_{rec} (#murad)");
	he->GetXaxis()->SetRangeUser(emin,emax);
	he->GetYaxis()->SetTitleOffset(1.3);
//	he->GetYaxis()->SetRangeUser(emin*0.9,1.1*emax);
//	he->GetYaxis()->SetLabelSize(0.02);
//	he->GetXaxis()->SetLabelSize(0.02);
	he->SetStats(0);
	
	xang->cd(2);
*/	hde->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Horizontal angle reconstruction @ %.0fm (%.0f #murad)",pos1,angl);
	hde->SetTitle(hetitle);
	hde->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	hde->GetYaxis()->SetTitle("#theta^{x}_{rec}-#theta^{x}_{rel} (#murad)");
	hde->GetXaxis()->SetRangeUser(emin,emax);
	hde->GetYaxis()->SetTitleOffset(1.2);
//	hde->GetXaxis()->SetLabelSize(0.02);
//	hde->GetYaxis()->SetLabelSize(0.02);
	hde->SetStats(0);


	TCanvas* yang = new TCanvas("yang","",1);
/*	yang->Divide(2,1);
	yang->cd(1);
	hy->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Vertical angle reconstruction @ %.0fm (%.0f #murad)",pos1,angl);
	hy->SetTitle(hetitle);
	hy->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	hy->GetYaxis()->SetTitle("Reconstructed angle #theta^{y}_{rec} (#murad)");
	hy->GetXaxis()->SetRangeUser(emin,emax);
	hy->GetYaxis()->SetTitleOffset(1.3);
//	hy->GetYaxis()->SetRangeUser(emin*0.9,1.1*emax);
//	hy->GetYaxis()->SetLabelSize(0.02);
//	hy->GetXaxis()->SetLabelSize(0.02);
	hy->SetStats(0);
	
	yang->cd(2);
*/	hdy->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	sprintf(hetitle,"Vertical angle reconstruction @ %.0fm (%.0f #murad)",pos1,angl);
	hdy->SetTitle(hetitle);
	hdy->GetXaxis()->SetTitle("Real energy loss E_{rel} (GeV)");
	hdy->GetYaxis()->SetTitle("#theta^{y}_{rec}-#theta^{y}_{rel} (#murad)");
	hdy->GetXaxis()->SetRangeUser(emin,emax);
	hdy->GetYaxis()->SetTitleOffset(1.2);
//	hdy->GetXaxis()->SetLabelSize(0.02);
//	hdy->GetYaxis()->SetLabelSize(0.02);
	hdy->SetStats(0);

	if (save) {
		char outfilename[50];

		sprintf(outfilename,"angle_x_reconstruction_%.0f_%.0f.eps",angl,pos1);
		xang->Print(outfilename);
		sprintf(outfilename,"angle_y_reconstruction_%.0f_%.0f.eps",angl,pos1);
		yang->Print(outfilename);
		
//		delete he;
//		delete hy;
		delete hde;
		delete hdy;
		delete xang;
		delete yang;
//		delete beam;
//		delete bea2;
	}

	return;
}

/// Returns the TGraph of the bias
TGraph * bias_graph(const float xmax, const float thmax, const float energyloss, const float pos , string filename, int side) {
	/// @param xmax : max X coordinate at IP, to be scanned (µm)
	/// @param thmax : max theta_X coordinate at IP, to be scanned (µm)
	/// @param energyloss : of the particle at IP (GeV)
	/// @param pos : RP position (m)
	/// @param filename : beam optics file
	/// @param side : direction of propagation
	const int n = 50;
	float x[n], th[n], rece[n], x1,y1,x2,y2;
	const float RP_dist = 8.; // distance between 2 planes of RP

	H_BeamLine* beam = new H_BeamLine(side,pos+RP_dist+1.);
	beam->fill(filename);
	H_RecRPObject r(pos,pos+RP_dist,beam);

	for(int i = 0; i < n; i++) {
		x[i] = i*xmax/((float)n-1);
		th[i] = i*thmax/((float)n-1);

		H_BeamParticle p1;
		p1.emitGamma(energyloss,0);

		p1.setPosition(PX+x[i],PY,-CRANG-th[i],0,0);
		p1.computePath(beam,1);

		p1.propagate(pos);
		x1 = p1.getX();
		y1 = p1.getY();

		p1.propagate(pos+RP_dist);
		x2 = p1.getX();
        	y2 = p1.getY();

		r.setPositions(x1,y1,x2,y2);
		rece[i] = r.getE(TM);
	}

	// removing energy bias due to nonlinearities :
	for(int i = n-1; i >= 0; i--) rece[i] -= rece[0];

	if(xmax )return new TGraph(n,x,rece);
	else return new TGraph(n,th,rece);
}

/// Reconstructed energy bias
void reconstruction_bias(float pos = 420, float xmax = 100., float thmax = 100., float energyloss = 100.,  float energyloss2 = 110., string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1) {
	/// @param pos : RP position1 from IP (m) 
	/// @param xmax : max X coordinate at IP, to be scanned (µm)
	/// @param thmax : max theta_X coordinate at IP, to be scanned (µm)
	/// @param energyloss : of the particle at IP (GeV), graph1
	/// @param energyloss2 : of the particle at IP (GeV), graph2
	/// @param filename : beam optics file
	/// @param side : direction of propagation
	char title[100];

	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = true;

	// x graphs
	sprintf(title,"Bias on E_{rec} at %.0fm from x(IP)",pos);
	TMultiGraph * mgx = new TMultiGraph("mgx",title);

	TGraph* gxx = bias_graph(xmax, 0, energyloss, pos, filename, side);
	gxx->SetMarkerStyle(22);
	gxx->SetMarkerColor(kRed);
	mgx->Add(gxx,"P");

	TGraph* gx1 = bias_graph(xmax, 0, energyloss2, pos, filename, side);
	gx1->SetMarkerStyle(23);
	gx1->SetMarkerColor(kBlue);
	mgx->Add(gx1,"P");

        TCanvas* tc = new TCanvas("tc1","bias",1);
        tc->cd();
        gPad->SetGrid();
        mgx->Draw("A");
        mgx->GetYaxis()->SetTitle("E_{rec} bias (GeV)");
        mgx->GetXaxis()->SetTitle("x at the IP (#mum)");

	TLegend * lef = new TLegend(0.14,0.63,0.42,0.78);
	sprintf(title,"E = %.0f GeV  ",energyloss);
	lef->AddEntry(gxx,title,"p");
	sprintf(title,"E = %.0f GeV  ",energyloss2);
	lef->AddEntry(gx1,title,"p");
	lef->Draw();
	lef->SetBorderSize(1);
	lef->SetFillColor(0);



	// theta graphs
	sprintf(title,"Bias on E_{rec} at %.0fm from #theta_{x}(IP)",pos);
	TMultiGraph * mgt = new TMultiGraph("mgt",title);

	TGraph* gth = bias_graph(0, thmax, energyloss, pos, filename, side);
	gth->SetMarkerStyle(22);
	gth->SetMarkerColor(kRed);
	mgt->Add(gth,"P");

        TGraph* gt2= bias_graph(0, thmax, energyloss2, pos, filename, side);
        gt2->SetMarkerStyle(23);
        gt2->SetMarkerColor(kBlue);
        mgt->Add(gt2,"P");


	TCanvas* tc2 = new TCanvas("tc2","bias",1);
	tc2->cd();
	gPad->SetGrid();
	mgt->Draw("A");
	mgt->GetYaxis()->SetTitle("E_{rec} bias (GeV)");
	mgt->GetXaxis()->SetTitle("#theta_{x} at the IP (#murad)");

	if(pos<400) mgx->GetYaxis()->SetRangeUser(-(mgt->GetYaxis()->GetXmax()),0);
	else mgx->GetYaxis()->SetRangeUser(0,mgt->GetYaxis()->GetXmax());

	mgt->Draw("A");
	lef->Draw();
	return;
}


///////////////////////////////////////////:

/// used in show_energy_res
float energy_res(float energy = 100, float q2 = 0., float detres = 10., int method = TM, int meanorres = 0, float pos = 420., const int N = 100, string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1) {
	/// \param energy : energy loss (\f$ GeV \f$ )
	/// \param q2 : momentum transfert (\f$ GeV^2 \f$ )
	/// \param detres : detector resolution in position (\f$ \mu m \f$ )
	/// \param method : for E reconstruction (TM or ACM)
	/// \param meanorres : choice of the return value : mean or RMS
	/// \param pos : RP distance from IP (m)
	/// \param filename : optics file
	/// \param side : direction of propagation

	using namespace TMath;

	extern int kickers_on;
	kickers_on = 1;

	extern bool relative_energy;
	relative_energy = true;
	
	float pos1 = pos;
	float pos2 = pos+8.;

	H_BeamLine* beam = new H_BeamLine(side,pos2+1);
	beam->fill(filename);

	H_RecRPObject r(pos1,pos2,beam);
	
	float x1s = 0;
	float x2s = 0;
	float y1s = 0;
	float y2s = 0;

	float rece[N];
	for(int i=0;i<N;i++) {

		H_BeamParticle p1;
		p1.setPosition(PX,PY,-CRANG,0,0);
		p1.smearPos();
		p1.smearAng();
		p1.emitGamma(energy,q2);
		p1.computePath(beam,1);

		p1.propagate(pos1);
		x1s = gRandom->Gaus(p1.getX(),detres);
		y1s = gRandom->Gaus(p1.getY(),detres);
		p1.propagate(pos2);
		x2s = gRandom->Gaus(p1.getX(),detres);
		y2s = gRandom->Gaus(p1.getY(),detres);

		r.setPositions(x1s,y1s,x2s,y2s);
		
		rece[i] = r.getE(method);
	}

//	delete beam;

	if(meanorres) return RMS(N,rece);
	else return Mean(N,rece);
}

void reconstruction_e_resolution(float pos = 420., float emin = 20., float emax = 120., int N = 1000) {

	TLegend* leg = new TLegend(0.7,0.3,0.99,0.5,"Detector resolution");
	float det_res1 = 5; // mum
	float det_res2 =30; // mum


	// first graphs : dE vs E
	// ------------------------
	char tit1[100];
	sprintf(tit1,"Energy resolution at %.0f m - trivial method",pos);
	TMultiGraph* meth1 = new TMultiGraph("met1",tit1);

	const int nepoints = 10;
	float ee[nepoints];  // e_loss : x-axis
	float de1[nepoints]; // e_res for 5 microns 
	//float de2[nepoints]; // e_res for 30 microns

	float ep1[3], /*ep2[3],*/ ep3[3]/*, ep4[3]*/;
	float q2[3] = {0.01,1.,10.};
	
	for(int i = 0; i < nepoints; i++) { 
		ee[i] = emin + i*(emax - emin)/(nepoints-1); 
	}

	// lines with detector resolution = 5;
	cout<<"dE vs E, detector resolution = " << det_res1 << " #mum"<<endl;
	for(int k = 0; k < 3; k++) {
		for(int i = 0; i < nepoints; i++) {
			de1[i] = energy_res(ee[i],-q2[k],det_res1,TM,1,pos,N);
			ep1[k] = de1[i];
		}
		TGraph* gr1 = new TGraph(nepoints,ee,de1);
		gr1->SetLineWidth(2);
		meth1->Add(gr1);
		sprintf(tit1,"%.0f #mum",det_res1);
		if(!k) leg->AddEntry(gr1,tit1,"l");
	}

	// lines with detector resolution = 30;
	cout<<"dE vs E, detector resolution = "<< det_res2 << " #mum"<<endl;
	for(int k = 0; k < 3; k++) {
		for(int i = 0; i < nepoints; i++) {
			de1[i] = energy_res(ee[i],-q2[k],det_res2,TM,1,pos,N);
		}
		TGraph* gr1 = new TGraph(nepoints,ee,de1);
		gr1->SetLineStyle(2);
		gr1->SetLineColor(kRed);
		gr1->SetLineWidth(2);
		meth1->Add(gr1);
		sprintf(tit1,"%.0f #mum",det_res2);
		if(!k) leg->AddEntry(gr1,tit1,"l");
	}


	// next graph : dE vs Q²
	// ------------------------
	sprintf(tit1,"Energy resolution at %.0f m - trivial method",pos);
	TMultiGraph* mqth1 = new TMultiGraph("mqt1",tit1);
	float qq[nepoints];
	float qmin = -5, qmax = 0.;
	float e2[3] = {20.,50,120.};

	for(int i = 0; i < nepoints; i++) {
		qq[i] = -(qmin + i*(qmax - qmin)/(nepoints-1));
	}

	// lines with detector resolution = 5;
	cout<<"dE vs Q², detector resolution = " << det_res1 << " #mum"<<endl;

	for(int k = 0; k < 3; k++) {
		for(int i = 0; i < nepoints; i++) {
			de1[i] = energy_res(e2[k],-qq[i],det_res1,TM,1,pos,N);
		}
		ep3[k] = de1[0];
		TGraph* gr1 = new TGraph(nepoints,qq,de1);
		gr1->SetLineWidth(2);
		mqth1->Add(gr1);
	}

	// lines with detector resolution = 30;
	cout<<"dE vs Q², detector resolution = " << det_res2 << " #mum"<<endl;

	for(int k = 0; k < 3; k++) {
		for(int i = 0; i < nepoints; i++) {
			de1[i] = energy_res(e2[k],-qq[i],det_res2,TM,1,pos,N);
		}
		TGraph* gr1 = new TGraph(nepoints,qq,de1);
		gr1->SetLineStyle(2);
		gr1->SetLineColor(kRed);
		gr1->SetLineWidth(2);
		mqth1->Add(gr1);
	}


	TCanvas* eresc = new TCanvas("eresc","Energy resolution",1);
	eresc->cd();
	meth1->Draw("AL");
	meth1->GetYaxis()->SetTitle("#deltaE (GeV)");
	meth1->GetXaxis()->SetTitle("E_{loss} (GeV)");


	TCanvas* eres2 = new TCanvas("eres2","Energy resolution",1);
	eres2->cd();
	mqth1->Draw("AL");
	mqth1->GetYaxis()->SetTitle("#deltaE (GeV)");
	mqth1->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
	leg->SetFillColor(0);
	leg->Draw();

	
	float max1,max3;
	max1 = meth1->GetYaxis()->GetXmax();
	max3 = mqth1->GetYaxis()->GetXmax();

	char cad[50];
	for(int i=0;i<3;i++) {
		eresc->cd();
		TPaveText* tp1 = new TPaveText(0.9*emax,ep1[i],1.1*emax,ep1[i]+max1/20.);
		tp1->SetBorderSize(0);
		tp1->SetFillColor(0);
		tp1->SetFillStyle(1);
		sprintf(cad,"Q^{2} = %.3f GeV^{2}",q2[i]);
		tp1->AddText(cad);
		tp1->Draw();

		eres2->cd();
		TPaveText* tp3 = new TPaveText(-0.9*qmin,ep3[i],-1.1*qmin,ep3[i]+max3/20.);
		tp3->SetBorderSize(0);
		tp3->SetFillColor(0);
		tp3->SetFillStyle(1);
		sprintf(cad,"E = %d GeV",(int)e2[i]);
		tp3->AddText(cad);
		tp3->Draw();

	}

	return;
}

int main() {
reconstruction_drawgrid();
return 0;
};
