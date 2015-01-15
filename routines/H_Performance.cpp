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

/// \file H_Performance.cpp
/// \brief Contains the routines for the computation of Hector performances (benchmarks)

// c++ #includes
#include <iostream>
#include <string>
using namespace std;

// ROOT #includes
#include "TCanvas.h" 
#include "TMultiGraph.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TROOT.h"

// local #includes
#include "H_Beam.h" 
#include "H_BeamLine.h"

void emitGamma_perf(const int N) {
	H_BeamParticle p;
	for (int i=0; i<N; i++) {
		p.emitGamma(1,-0.01);
	}
}

/// subroutine used by performance_compute()
void propagate_perf(int & n_el, float length, string filename, int side, const int Nparticles=100, bool nonlinear=false, bool draw=false) {

        gROOT->SetStyle("Plain");
        H_Beam mine;
        mine.createBeamParticles(Nparticles);
        cout << mine.getNumberOfBeamParticles() << " particles have been created." << endl;

        H_BeamLine * beamline;
        int direction =  (side<0)?-1:1;
	beamline = new H_BeamLine(side,length);
	beamline->fill(filename);
        beamline->offsetElements(140,-1*(direction)*0.080);
        beamline->offsetElements(146,-1*(direction)*0.097);
        mine.computePath(beamline,nonlinear);
	n_el = beamline->getNumberOfElements();

	if (draw) {
		TCanvas* cccp = new TCanvas("cccp","the particle test canvas",1);
		TMultiGraph * beamx = mine.drawBeamX(kBlack);
		TMultiGraph * beamy = mine.drawBeamY(kBlack);
		cccp->Divide(1,2);
		cccp->cd(1);
		beamx->Draw("AL");
		beamx->GetXaxis()->SetTitle("s [m]");
		beamx->GetYaxis()->SetTitle("beam width [m]");
		beamline->drawX(beamx->GetYaxis()->GetXmin(),beamx->GetYaxis()->GetXmax());
		beamx->Draw("L");

		cccp->cd(2);
		beamy->Draw("AL");
		beamy->GetXaxis()->SetTitle("s [m]");
		beamy->GetYaxis()->SetTitle("beam width [m]");
		beamline->drawY(beamy->GetYaxis()->GetXmin(),beamy->GetYaxis()->GetXmax());
		beamy->Draw("L");
	}
        delete beamline;
} 

/// Benchmark for brute propagation of protons by Hector
void performance_compute(const int Nparticles=10000, int steps = 10, string filename="data/LHCB1IR5_v6.500.tfs", int side=1) {
	/// @param Nparticles : beam content
	/// @param steps : number of points for the linear fit
	/// @param filename : optics source file
	/// @param side : direction of propagation (forward 1, backward -1)
	float length =10.;
	int nel = 0;
	float * x = new float[steps];
	float * y = new float[steps];
	float * y1 = new float[steps];
	float * n = new float[steps];
	TStopwatch mywatch;
	for (int i=0; i<steps; i++) {
		length = i==0?steps:50.*i;
		cout << Nparticles*(i+1)/1000 <<"E3 : ";
		mywatch.Start();
		propagate_perf(nel,length,filename,side,Nparticles,false);
		mywatch.Stop();
		cout << mywatch.CpuTime() << " " << mywatch.RealTime() << endl;
		x[i]=length;
		y[i]=mywatch.CpuTime();
		n[i]=(float)nel;
	}

        for (int i=0; i<steps; i++) {
                length = i==0?steps:50.*i;
                cout << Nparticles*(i+1)/1000 <<"E3 : ";
                mywatch.Start();
                propagate_perf(nel,length,filename,side,Nparticles,true);
                mywatch.Stop();
                cout << mywatch.CpuTime() << " " << mywatch.RealTime() << endl;
                x[i]=length;
                y1[i]=mywatch.CpuTime();
                n[i]=(float)nel;
        }

	char title[100];
	sprintf(title,"Computing time for %d particles",Nparticles);
	TCanvas * mycan = new TCanvas("mycan",title);
	mycan->Divide(2);
	mycan->cd(1);

	TMultiGraph *mgr1 = new TMultiGraph("mgr1",title);
	TGraph *pgr = new TGraph(steps,x,y);
	pgr->SetMarkerStyle(24);
	TGraph *pgr1 = new TGraph(steps,x,y1);
	pgr1->SetMarkerStyle(20);
	mgr1->Add(pgr,"P");
	mgr1->Add(pgr1,"P");
	mgr1->Draw("AP*");
	mgr1->GetXaxis()->SetTitle("Beam length (m)");
	mgr1->GetYaxis()->SetTitle("CPU time (s)");
	mgr1->GetYaxis()->SetTitleOffset(1.2);
	mgr1->SetTitle(title);
	gPad->SetGrid();

	mycan->cd(2);
	TGraph *pg2 = new TGraph(steps,n,y);
	TGraph *pg3 = new TGraph(steps,n,y1);
	TMultiGraph *mgr = new TMultiGraph("mgr","");

	TF1 *f1 = new TF1("f1","[0]+[1]*x",0.,length);
	f1->SetLineStyle(2);
	f1->SetLineWidth(1);
	pg2->Fit("f1","");

        TF1 *f2 = new TF1("f2","[0]+[1]*x",0.,length);
	f2->SetLineWidth(1);
        pg3->Fit("f2","");

	mgr->Add(pg2,"P");
	pg2->SetMarkerStyle(24);
	mgr->Add(pg3,"P");
	pg3->SetMarkerStyle(20);
	mgr->Draw("A");
        mgr->GetXaxis()->SetTitle("Number of optical elements");
        mgr->GetYaxis()->SetTitle("CPU time (s)");
        mgr->GetYaxis()->SetTitleOffset(1.2);
        mgr->SetTitle("");
	gPad->SetGrid();
	
}

int main() {
	performance_compute();
}
