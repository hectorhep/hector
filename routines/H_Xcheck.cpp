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

/// \file H_Xcheck.cpp
/// \brief Routines for mad-X validation.

// C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

// ROOT #includes
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TROOT.h"

// local #includes
#include "H_AbstractBeamLine.h"
#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_Beam.h"
#include "H_RectangularDipole.h"
using namespace std;

/// Draws the beta functions
void xcheck_beta(double length = 500., int number_of_points = 100, char * filename = "data/LHCB1IR5_v6.500.tfs", int side = 1, char * add2title = " ", const int crang_sign=-1, const bool save=false, char * outfilename="") {
	/// @param length is the length of the s axis
	/// @param number_of_points is the content of the graphs from Hector
	/// @param filename is the beamline used
	/// @param side gives the direction of the beam, in the beamline (forward 1, backward -1)
	/// @param add2title is an extra string one can add to the title of the graph
	/// @param crang_sign is the direction of the (half) crossing angle at IP
	/// @param save : boolean
	/// @param outfilename : to be written if "save" is true
	/// usage : testBeam_Betafunctions(546,200,"twiss_ip5_b1_v6.5.txt",1,"(beam 1)");

	extern int relative_energy;
	relative_energy  = false;
	if(relative_energy) {
		cout << "You should be in absolute energy" << endl;
		return;
	}

	gROOT->SetStyle("Plain");
        H_Beam mine;
	TRandom3 *r = new TRandom3();
        mine.createBeamParticles(10000,MP,QP,r);
        cout << mine.getNumberOfBeamParticles() << " particles have been created." << endl;

	H_BeamLine* beamline = new H_BeamLine(side,length);
        beamline->fill(filename);
	beamline->showElements();
	beamline->offsetElements(120,0.097*crang_sign);
        mine.computePath(beamline);

	cout << "Emittance x : " << mine.getEmittanceX() << endl;
	cout << "Emittance y : " << mine.getEmittanceY() << endl;

	char title[50];
	sprintf(title,"#beta functions %s",add2title);
        TCanvas * beta;
        if(gROOT->FindObject("beta")) {

                beta = (TCanvas*) gROOT->FindObject("beta");
                delete beta;
        }
        beta = new TCanvas("beta","",0,0,700,500);
        TMultiGraph * betag = new TMultiGraph("betag",title);
	TGraph * beta_x, *beta_y;
//	TGraphErrors * beta_x, *beta_y;
	beta_x = mine.getBetaX(length-0.1,number_of_points);
        beta_y = mine.getBetaY(length-0.1,number_of_points);
        TGraph * beta_xth = beamline->getBetaX();
        TGraph * beta_yth = beamline->getBetaY();

        beta->cd();
	beta_x->SetMarkerStyle(25);
	beta_x->SetMarkerSize(0.5);
	beta_y->SetMarkerStyle(25);
	beta_y->SetMarkerSize(0.5);
	beta_y->SetMarkerColor(kRed);
	beta_xth->SetLineStyle(1);
	beta_yth->SetLineStyle(2);

        betag->Add(beta_x);
        betag->Add(beta_y);
        betag->Add(beta_xth);
        betag->Add(beta_yth);

        betag->Draw("A");
        betag->GetXaxis()->SetTitle("s [m]");
        betag->GetYaxis()->SetTitle("#beta [m]");
        betag->GetYaxis()->SetTitleOffset(1.31);
//        betag->Draw("C");

	beta_x->Draw("P");
	beta_y->Draw("P");
	beta_xth->Draw("C");
	beta_yth->Draw("C");

        TLegend* leg=new TLegend(0.75,0.55,0.95,0.78," ");
        leg->AddEntry(beta_x,"#beta_{x} - Hector","p");
        leg->AddEntry(beta_y,"#beta_{y} - Hector","p");
        leg->AddEntry(beta_xth,"#beta_{x} - MAD-X","l");
        leg->AddEntry(beta_yth,"#beta_{y} - MAD-X","l");
	leg->SetBorderSize(1);
	leg->SetFillColor(0);
        leg->Draw();
	if (save) {
		sprintf(outfilename,"%s.eps",outfilename);
		beta->Print(outfilename,"eps");
		delete beta_x ; delete beta_y ; delete beta_xth ; delete beta_yth ; delete betag; delete leg; delete beta; delete beamline;
	} 

} // xcheck_beta

/// Draws the relative positions
void xcheck_position(double length=500., int side=1, int crang_sign = -1, int ip = 5, int gr = 0, const string& add2title = "", int number_of_points = 100, const string& filename="data/LHCB1IR5_v6.500.tfs", const bool save=false, const string& outfilename="") {
	/// @param length : distance of propagation
	/// @param side : direction of propagation (forward 1, backward -1)
	/// @param crang_sign : direction of the (half) crossing angle at IP
	/// @param ip : 1 or 5
	/// @param gr : choose the output type (0 or non 0)
	/// @param add2title : string to be added to the title of the graph
	/// @param number_of_points : beam content
	/// @param filename : optics source file
	/// @param save : boolean
	/// @param outfilename : file to be written to

	gROOT->SetStyle("Plain");

	double refpos[number_of_points];
	double abspos[number_of_points];
	double relpos[number_of_points];
	double thes[number_of_points];
	
	const string ipname = (ip==1) ? "IP1" : "IP5";

	// 1) propagation without crang or kickers : reference path
	// crang_off ,x_ini = 0, kickers_off, absolute energy,  energy_sdip

	extern int kickers_on;
	kickers_on = 0;
	H_BeamLine* beamline = new H_BeamLine(side,length);
	beamline->fill(filename,-side*crang_sign,ipname);
	beamline->offsetElements(120,0.097*crang_sign);
	// hg <-> data from Mad-X
	TGraph* hg = (ip==1) ? beamline->getRelY() : beamline->getRelX();

        ofstream backup("m7000.txt");
        const int NPoints = hg->GetN();
//	double * ms = new double[NPoints];
//	double * mx = new double[NPoints];
	Double_t * ms = hg->GetX();
	Double_t * mx = hg->GetY();
	for(int i=0; i<NPoints; i++)	
	        backup << setw(15) << ms[i] << setw(15) << -mx[i] << endl;
	backup.close();
	cout << "The file 'm7000.txt' has been created\n";

	H_BeamParticle* p = new H_BeamParticle;
	p->setPosition(0.,0.,0.,0.,0.);
	p->computePath(beamline,1);
	p->propagate(10.);
	TGraph* nocrang = (ip==1) ? p->getPath(1,0) : p->getPath(0,0);
	for(int i=0; i<number_of_points; i++) {
		thes[i] = i*0.99*length/((double)(number_of_points-1));
		p->propagate(thes[i]);
		refpos[i] = (ip==1) ? p->getY() : p->getX();
	}
	
	beamline->showElements(HKICKER);	
	delete beamline;
	delete p;

	// 2) absolute propagation position
	// crang_on ,x_ini = -500, kickers_on, absolute energy,  energy_sdip

	kickers_on = 1;
	H_BeamLine* beamlin2 = new H_BeamLine(side,length);
	beamlin2->fill(filename,-side*crang_sign,ipname);
	beamlin2->offsetElements(120,0.097*crang_sign);

	H_BeamParticle* p2 = new H_BeamParticle;
	if(ip==1) p2->setPosition(0.,500.,0.,CRANG*crang_sign,0.);
	else      p2->setPosition(-500.,0.,CRANG*crang_sign,0.,0.);
	p2->computePath(beamlin2,1);
	p2->propagate(10.);
	TGraph* wicrang = (ip==1) ? p2->getPath(1,1) : p2->getPath(0,1);
	
	for(int i=0; i<number_of_points; i++) {
		p2->propagate(thes[i]);
		abspos[i] = (ip==1) ? p2->getY() : p2->getX();
		relpos[i] = (ip==1) ? -(p2->getY() - refpos[i]) : -(p2->getX() - refpos[i]);
	}

	delete beamlin2;
	delete p2;

	// 3) drawing and praying...
	TCanvas * can;
        if(gROOT->FindObject("can")) {

                can = (TCanvas*) gROOT->FindObject("can");
                delete can;
        }
	can = new TCanvas("can","",0,0,700,500);
	can->cd();
	TGraph* tg = new TGraph(number_of_points,thes,abspos);
	TGraph* rg = new TGraph(number_of_points,thes,refpos);
	TGraph* cg = new TGraph(number_of_points,thes,relpos);
	if(gr) {
		nocrang->Draw("A*");
		wicrang->Draw("L");
		tg->Draw("A*");
		rg->Draw("L");
	} else {
		cg->SetMarkerStyle(25);
		cg->SetMarkerSize(0.5);
		cg->SetMarkerColor(kRed);
		cg->Draw("AP");
		char title[50];
		sprintf(title,"Relative position to ideal path %s",add2title.c_str());
		cg->SetTitle(title);
		cg->GetYaxis()->SetTitle("#deltax (#mum)");
		cg->GetXaxis()->SetTitle("s (m)");
		cg->GetYaxis()->SetTitleOffset(1.3);
		hg->SetLineStyle(0);
		hg->SetLineColor(kRed);
		hg->Draw("L");

		TLegend* leg = new TLegend(0.7,0.6,0.95,0.7);
		leg->AddEntry(cg,"#deltax -Hector","p");
		leg->AddEntry(hg,"#deltax -MadX","l");
		leg->SetFillColor(0);
		leg->Draw();
		if (save) {
			can->SetGridx();
			can->SetGridy();
			can->Print(outfilename.c_str());
			delete leg;
		}
	}

	if (save) {
		delete hg;
		delete nocrang;
		delete wicrang;
		delete cg;
		delete tg;
		delete rg;
		delete can;
	}
	return;
} // xcheck_position


/// Draws the beta functions
void xcheck_dispersion(double length = 500., int number_of_points = 100, const string& filename = "data/LHCB1IR5_v6.500.tfs", int side = 1, char * add2title = " ", const int crang_sign=-1, const bool save=false, string outfilename="") {
	/// @param length is the length of the s axis
	/// @param number_of_points is the content of the graphs from Hector
	/// @param filename is the beamline used
	/// @param side gives the direction of the beam, in the beamline (forward 1, backward -1)
	/// @param add2title is an extra string one can add to the title of the graph
	/// @param crang_sign is the direction of the (half) crossing angle at IP
	/// @param save : boolean
	/// @param outfilename : to be written if "save" is true


	extern int relative_energy;
	relative_energy=false;
	if(relative_energy) {
		cout << "You should be in absolute energy" << endl;
		return;
	}

	gROOT->SetStyle("Plain");
        H_Beam mine;
        mine.createBeamParticles(10000);
        cout << mine.getNumberOfBeamParticles() << " particles have been created." << endl;

	H_BeamLine* beamline = new H_BeamLine(side,length);
        beamline->fill(filename);
	beamline->showElements();
	beamline->offsetElements(120,0.097*crang_sign);
        mine.computePath(beamline);

	char title[50];
	sprintf(title,"#beta functions %s",add2title);
        TCanvas * disp = new TCanvas("disp","");
        TMultiGraph * dispg = new TMultiGraph("dispg",title);
	TGraph * disp_x, *disp_y;
//	TGraphErrors * disp_x, *disp_y;
	disp_x = mine.getBetaX(length-0.1,number_of_points);
        disp_y = mine.getBetaY(length-0.1,number_of_points);
        TGraph * disp_xth = beamline->getDX();
        TGraph * disp_yth = beamline->getDY();

        disp->cd();
	disp_x->SetMarkerStyle(25);
	disp_x->SetMarkerSize(0.5);
	disp_y->SetMarkerStyle(25);
	disp_y->SetMarkerSize(0.5);
	disp_y->SetMarkerColor(kRed);
	disp_xth->SetLineStyle(1);
	disp_yth->SetLineStyle(2);

        dispg->Add(disp_x);
        dispg->Add(disp_y);
        dispg->Add(disp_xth);
        dispg->Add(disp_yth);

        dispg->Draw("A");
        dispg->GetXaxis()->SetTitle("s [m]");
        dispg->GetYaxis()->SetTitle("#disp [m]");
        dispg->GetYaxis()->SetTitleOffset(1.31);
//        dispg->Draw("C");

//	disp_x->Draw("P");
//	disp_y->Draw("P");
	disp_xth->Draw("L");
	disp_yth->Draw("L");

        TLegend* leg=new TLegend(0.75,0.55,0.95,0.78," ");
        leg->AddEntry(disp_x,"#disp_{x} - Hector","p");
        leg->AddEntry(disp_y,"#disp_{y} - Hector","p");
        leg->AddEntry(disp_xth,"#disp_{x} - MAD-X","l");
        leg->AddEntry(disp_yth,"#disp_{y} - MAD-X","l");
	leg->SetBorderSize(1);
	leg->SetFillColor(0);
        leg->Draw();
	if (save) {
		outfilename +=".eps";
		disp->Print(outfilename.c_str(),"eps");
		delete disp_x ; delete disp_y ; delete disp_xth ; delete disp_yth ; delete dispg; delete leg; delete disp; delete beamline;
	} 

} // xcheck_dispersion

void xcheck_theory(const double lowe=0.1, const double highe=1.0, const double theta_x=0.0001, const double R=1000.) {

	double xp=theta_x; // in rad
	const double l=10.; // in m
	char title[50];
	TMultiGraph * tmg = new TMultiGraph("tmg","Deflection by a rectangular dipole : position");
	TMultiGraph * thf = new TMultiGraph("thf","Deflection by a rectangular dipole : angle");
	TMultiGraph * diff = new TMultiGraph("diff","Position error on the deflection by a rectangular dipole");
	TMultiGraph * diffa = new TMultiGraph("diffa","Angle error on the deflection by a rectangular dipole");
	sprintf(title,"R=%dm, L=%dm",(int)R,(int)l);
        TLegend * leg=new TLegend(0.75,0.55,0.95,0.78,title);
        TLegend * ldi=new TLegend(0.75,0.22,0.95,0.45,title);

	// Theory for the dx
	TF1 * eq2 = new TF1("eq2","-([1]*[1]/[0])/(1 + sqrt(1-([1]/[0])*([1]/[0]))) + [0]*x*sqrt(1-[3]*[3]) - [0]*x*sqrt(1 - (([1]/([0]*x) - [3])*([1]/([0]*x) - [3])))",lowe,highe);
	eq2->SetParameter(0,R); // R
	eq2->SetParameter(1,l); // l
	eq2->SetParameter(3,xp); // x'
	//eq2->Draw();

	// Theory for the dtheta
	TF1 * eq1 = new TF1("eq1","- asin([1]/(x*[0]) - [3]) + asin([1]/[0] - [3])",lowe,highe);
	eq1->SetParameter(0,R); // R
	eq1->SetParameter(1,l); // l
	eq1->SetParameter(3,xp); // x'


	// theory
	const int Nb = 50;

	for (int j=0; j<2; j++) {
		xp *= (9*j+1);
		double * et = new double[Nb]; // list of abscisses : percent of p0
		double * xt = new double[Nb]; // horizontal shift
		double * tt = new double[Nb]; // theta_f
		for (int i=0; i<Nb; i++) {
		// what is called here xp is in fact x'=dx/ds=sin(theta_x); 
			eq2->SetParameter(3,xp); 	
			et[i]=lowe + (highe - lowe)*i/(double)(Nb-1);
			xt[i] = eq2->Eval(et[i])*1000.; // in mm
			tt[i] = eq1->Eval(et[i])*1000.; // in mrad
		}
	
		TGraph *th1 = new TGraph(Nb,et,xt); // delta x
		TGraph *th2 = new TGraph(Nb,et,tt); // delta theta_x
		sprintf(title,"Theory : %.2f mrad",1000*xp);
	        leg->AddEntry(th1,title,"l");
		delete [] et;
		delete [] xt;
		th1->SetLineColor(kRed+j);
		th2->SetLineColor(kRed+j);
		th1->SetLineStyle(j+1);
		th2->SetLineStyle(j+1);
		tmg->Add(th1,"l");
		thf->Add(th2,"l");
	}

	// Hector
        extern bool relative_energy;
        relative_energy = true;

        H_AbstractBeamLine* beam = new H_AbstractBeamLine(10.+1.0e-7);
        H_RectangularDipole* rdip = new H_RectangularDipole(0,1/(double)R,10.);
        beam->add(rdip);
	beam->showElements();
	
	const int N = 9;

	xp=theta_x;
	for (int j=0; j<2; j++) {
		xp *= (9*j+1);
		cout << "xp = " << xp << endl;
		double * e = new double[N]; // percentage : from 0.9 to 1.0
		double * x = new double[N]; // hector : x 
		double * s = new double[N]; // hector : sin(theta)
		double * t = new double[N]; // theory pos
		double * ta = new double[N]; // theory angle
		double * d = new double[N]; // difference
		double * da = new double[N]; // difference angle
		for (int i=0; i<N; i++) {
	        	H_BeamParticle p1;
			e[i]=lowe + (highe - lowe) * i /(double)(N-1);

			p1.setPosition(0,0,1E6*xp,0,0);
        		p1.setE(7000*e[i]);
		        p1.computePath(beam,1);
			p1.propagate(10.+1.0e-8);
			x[i] = -p1.getX()/1000.; // x_hector = x[i] in mm
			s[i] = p1.getTX()/1000. -xp*1000.; // thetax_hector = s[i] in mrad
			eq2->SetParameter(3,xp);
			t[i] = eq2->Eval(e[i])*1000.;  //  x_theory = t[i] in mm
			d[i] = (x[i] - t[i])*1000.;    //  d[i] = x_theory - x_hector in mm
			ta[i] = eq1->Eval(e[i])*1000;  //  theta_theory = ta[i] in mrad
			da[i] = (s[i] - ta[i])*1000.;	       //  da[i] = thetax_hector - thetax_theory, in mrad

			cout<<da[i]/1000.<<" "<<s[i]<<" "<<ta[i]<<endl;
			cout<<s[i]*1000.-ta[i]*1000.<<" "<<s[i]*1000.<<" "<<ta[i]*1000.<<endl;

		}
		TGraph * gr = new TGraph(N,e,x);
		TGraph * h2 = new TGraph(N,e,s);
		TGraph * di = new TGraph(N,e,d);
		TGraph * dia = new TGraph(N,e,da);
		sprintf(title,"Hector : %.2f mrad",1000*xp);
	        leg->AddEntry(gr,title,"p");
		sprintf(title,"%.2f mrad",1000*xp);
		ldi->AddEntry(di,title,"p");
		delete [] e;
		delete [] x;
		delete [] t;
		delete [] d;
		delete [] da;
		gr->SetMarkerColor(kRed+j);
		gr->SetMarkerStyle(25+j);
                di->SetMarkerColor(kRed+j);
                di->SetMarkerStyle(25+j);
                h2->SetMarkerColor(kRed+j);
                h2->SetMarkerStyle(25+j);
                dia->SetMarkerColor(kRed+j);
                dia->SetMarkerStyle(25+j);

		tmg->Add(gr,"p");
		thf->Add(h2,"p");
		diff->Add(di,"p");
		diffa->Add(dia,"p");
	}

        leg->SetBorderSize(1);
        leg->SetFillColor(0);

	// Deflection in position
        TCanvas * ca;
        if(gROOT->FindObject("ca")) {

                ca = (TCanvas*) gROOT->FindObject("ca");
                delete ca;
        }
        ca = new TCanvas("ca","",0,0,700,500);
	ca->cd();
	tmg->Draw("a");	
	tmg->GetXaxis()->SetTitle("p/p_{0}");
	tmg->GetYaxis()->SetTitle("#Deltax (mm)");
        leg->Draw();
	gPad->SetGrid();

	// Deflection in angle
        TCanvas * cc;
        if(gROOT->FindObject("cc")) {

                cc = (TCanvas*) gROOT->FindObject("cc");
                delete cc;
        }
	cc = new TCanvas("cc");
	cc->cd();
	thf->Draw("a");
	thf->GetXaxis()->SetTitle("p/p_{0}");
	thf->GetYaxis()->SetTitle("#Delta #theta_{x} (mrad)");
        leg->Draw();
	gPad->SetGrid();

	// Error in angle
        TCanvas * cdia;
        if(gROOT->FindObject("cdia")) {

                cdia = (TCanvas*) gROOT->FindObject("cdia");
                delete cdia;
        }
        cdia = new TCanvas("cdia");
	cdia->cd();
	diffa->Draw("a");
	diffa->GetXaxis()->SetTitle("p/p_{0}");
	diffa->GetYaxis()->SetTitle("#Delta #theta^{x}_{Hector} - #Delta #theta^{x}_{th} (#murad)");
        ldi->Draw();
	gPad->SetGrid();

 	// Error in position
        TCanvas * cb;
        if(gROOT->FindObject("cb")) {

                cb = (TCanvas*) gROOT->FindObject("cb");
                delete cb;
        }
        cb = new TCanvas("cb");
        cb->cd();
        diff->Draw("a");
        diff->GetXaxis()->SetTitle("p/p_{0}");
        diff->GetYaxis()->SetTitle("#Delta x_{Hector} - #Delta x_{th} (#mum)");
        ldi->SetBorderSize(1);
        ldi->SetFillColor(0);
        ldi->Draw();
	gPad->SetGrid();

}

